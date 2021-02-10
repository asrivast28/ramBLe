#!/usr/bin/env python

##
# @file mnets_lemontree.py
# @brief Script for running lemon-tree using mnets arguments.
# @author Ankit Srivastava <asrivast@gatech.edu>
#
# Copyright 2021 Georgia Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from collections import OrderedDict
import os.path as path
from tempfile import NamedTemporaryFile
from time import time


def parse_args():
    '''
    Parse command line arguments.
    '''
    from argparse import ArgumentParser
    from os import makedirs
    from shutil import copyfile

    LEMONTREE_DIR = path.join(path.expanduser('~'), 'data', 'lemon-tree')
    LEMONTREE_RUN = path.join(LEMONTREE_DIR, 'LemonTree', 'run_task.sh')

    parser = ArgumentParser(description='Run lemon-tree using mnets configuration')
    parser.add_argument('executable', nargs='?', metavar='RUN', default=LEMONTREE_RUN, help='Path to the lemon-tree task executable')
    parser.add_argument('-n', '--nvars', metavar='N', type=int, help='Number of variables in the dataset')
    parser.add_argument('-m', '--nobs', metavar='M', type=int, help='Number of observations in the dataset')
    parser.add_argument('-f', '--file', metavar='FILE', type=str, required=True, help='Name of the file from which dataset is to be read')
    parser.add_argument('-c', '--colobs', action='store_true', help='The file contains observations in columns')
    parser.add_argument('-s', '--separator', metavar='CHAR', type=str, help='Delimiting character in the file')
    parser.add_argument('-v', '--varnames', action='store_true', help='The file contains variable names')
    parser.add_argument('-i', '--indices', action='store_true', help='The file contains observation indices')
    parser.add_argument('-a', '--algorithm', metavar='NAME', type=str, help='Name of the algorithm to be used')
    parser.add_argument('-o', '--outdir', metavar='DIR', type=str, default='.', help='Name of the directory to which the output files should be written')
    parser.add_argument('-g', '--config', metavar='FILE', type=str, help='JSON file with algorithm specific configurations')
    parser.add_argument('--prefix', type=str, help='Prefix for the generated output files')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()

    if not path.exists(args.file):
        raise RuntimeError('Couldn\'t find the data file')
    if not (args.colobs and args.varnames and args.indices):
        raise RuntimeError('Only files with "-c -v -i" options are currently supported')
    if not args.config:
        args.config = 'lemontree_configs.json'
    if not path.exists(args.config):
        raise RuntimeError('Couldn\'t find the algorithm configuration file')
    if not path.isdir(args.outdir):
        try:
            makedirs(args.outdir)
        except:
            raise RuntimeError('Output directory doesn\'t exist and could not be created')
    copyfile(args.config, path.join(args.outdir, 'configs.json'))
    return args


def read_task_configs(configfile, datafile):
    '''
    Read configurations for all the tasks from the JSON file.
    '''
    from json import load

    with open(configfile, 'r') as cf:
        all_configs = load(cf, object_pairs_hook=OrderedDict)
    task_configs = []
    for task, configs in all_configs.items():
        configs['data_file'] = datafile
        configs.move_to_end('data_file', last=False)
        configs['task'] = task
        configs.move_to_end('task', last=False)
        task_configs.append(configs)
    return task_configs


def read_dataset(name, sep, colobs, varnames, indices):
    '''
    Read dataset from the given CSV.
    '''
    from pandas import read_csv

    header = None
    index = False
    if colobs:
        if indices:
            header = 0
        if varnames:
            index = 0
    else:
        if varnames:
            header = 0
        if indices:
            index = 0
    dataset = read_csv(name, sep=sep, header=header, index_col=index)
    if not colobs:
        dataset = dataset.T
    return dataset


def write_variables(name, sep, colobs, varnames, indices):
    '''
    Read the variables from the given CSV file
    and write them to a temporary file.
    '''
    dataset = read_dataset(name, sep, colobs, varnames, indices)
    varfile = None
    with NamedTemporaryFile(suffix='.txt', delete=False) as vf:
        vf.write('\n'.join(dataset.index).encode())
        varfile = vf.name
    return varfile


def run_task(executable, configs, verbose, previous=None):
    '''
    Function for running any lemon-tree task.
    '''
    import subprocess
    from sys import stdout

    output_file = configs['output_file']
    if previous is not None:
        configs['cluster_file'] = previous
    arguments = ' '.join('-%s %s' % (name, val) for name, val in configs.items())
    arguments = executable + ' ' + arguments
    if verbose:
        print(arguments)
    process = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
    if verbose:
        for line in iter(process.stdout.readline, b''):
            line = line.decode('utf-8')
            print(line, end='')
    stdout.flush()
    process.communicate()
    if process.returncode != 0:
        raise RuntimeError('Error while executing the task %s' % configs['task'])
    return output_file


def run_ganesh(executable, configs, verbose):
    '''
    Function for running GaneSH with some
    task-specific changes.
    '''
    init_seed = configs['seed']
    num_runs = configs['num_runs']
    configs['num_runs'] = 1
    output_file = configs['output_file']
    with open(output_file, 'w') as of:
        for r in range(num_runs):
            run_file = output_file + '.' + str(r)
            configs['seed'] = init_seed + r
            configs['output_file'] = run_file
            configs['burn_in'] = configs['num_steps'] - 1
            run_task(executable, configs, verbose)
            of.write(run_file + '\n')
    return output_file


def run_lemontree(executable, configs, outdir, prefix, verbose):
    '''
    Function for running lemon-tree end-to-end.
    '''
    runtimes = OrderedDict()
    prev_output = None
    for task_configs in configs:
        output_file = task_configs['output_file']
        if not output_file:
            output_file = NamedTemporaryFile(suffix='.txt', delete=True).name
        if prefix:
            output_file = prefix + output_file
        output_file = path.join(outdir, output_file)
        task_configs['output_file'] = output_file
        t = time()
        if task_configs['task'] == 'ganesh':
            prev_output = run_ganesh(executable, task_configs, verbose)
        else:
            prev_output = run_task(executable, task_configs, verbose, prev_output)
        runtimes[task_configs['task']] = time() - t
        if verbose:
            print()
    return runtimes


def main():
    '''
    Main function.
    '''
    args = parse_args()
    task_configs = read_task_configs(args.config, args.file)
    t = time()
    for configs in task_configs:
        # Write a regulator file with all the variables if one does not exist
        if configs['task'] == 'regulators' and not configs['reg_file']:
            configs['reg_file'] = write_variables(args.file, args.separator, args.colobs, args.varnames, args.indices)
    runtimes = run_lemontree(args.executable, task_configs, args.outdir, args.prefix, args.verbose)
    print('Time taken in learning the network: %.6f' % (time() - t))
    message = {
        'ganesh'         : 'the GaneSH run',
        'tight_clusters' : 'consensus clustering',
        'regulators'     : 'learning the modules',
        }
    for task, rt in runtimes.items():
        print('Time taken in %s: %.6f' % (message[task], rt))


if __name__ == '__main__':
    main()
