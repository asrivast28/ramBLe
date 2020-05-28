#!/usr/bin/env python

##
# @file run_experiments.py
# @brief Script for running the experiments.
# @author Ankit Srivastava <asrivast@gatech.edu>
#
# Copyright 2020 Georgia Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from collections import OrderedDict
from itertools import product
import os
import os.path
from os.path import join
import sys
from tempfile import NamedTemporaryFile


small_datasets = OrderedDict([
    #(name        , (-f, -n, -m, -s, -c, -v, -i)),
    ('lizards'    , ('test/lizards.csv', 3, 409, ',', False, True, False)),
    ('coronary'   , ('test/coronary.csv', 6, 1841, ',', False, True, False)),
    ('asia'       , ('test/asia.csv', 8, 5000, ',', False, True, False)),
    ('child'      , ('test/child.csv', 20, 10000, ' ', True, False, False)),
    ('insurance'  , ('test/insurance.csv', 27, 10000, ' ', True, False, False)),
    ('mildew'     , ('test/mildew.csv', 30, 10000, ' ', True, False, False)),
    ('alarm'      , ('test/alarm.csv', 37, 10000, ' ', True, False, False)),
    ])

big_datasets = OrderedDict([
    #(name        , (-f, -n, -m, -s, -c, -v, -i)),
    ('yeast'      , ('data/yeast/yeast_data_discretized.tsv', 5716, 2577, ' ', True, True, True)),
    ('development', ('data/athaliana/athaliana_development_discretized.tsv', 18373, 5102, ' ', True, True, True)),
    ('complete'   , ('data/athaliana/athaliana_complete_discretized.tsv', 18380, 16838, ' ', True, True, True)),
    ])

all_datasets = OrderedDict(list(small_datasets.items()) + list(big_datasets.items()))

all_algorithms = [
    'gs',
    'iamb',
    'inter.iamb',
    ]

all_processes = [
    ]
for power in range(0, 11):
    all_processes.append(2 ** power)

ppn_mappings = OrderedDict([
    (16, '1:2:3:4:5:6:7:8:13:14:15:16:17:18:19:20'),
    (18, '1:2:3:4:5:6:7:8:9:13:14:15:16:17:18:19:20:21'),
    (20, '1:2:3:4:5:6:7:8:9:10:13:14:15:16:17:18:19:20:21:22'),
    (22, '1:2:3:4:5:6:7:8:9:10:11:13:14:15:16:17:18:19:20:21:22:23'),
    (24, '0:1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:23'),
    ])

NUM_REPEATS = 5


def parse_args():
    '''
    Parse command line arguments.
    '''
    import argparse
    from os.path import expanduser, realpath

    parser = argparse.ArgumentParser(description='Run scaling experiments')
    parser.add_argument('-b', '--basedir', metavar='DIR', type=str, default=realpath(join(expanduser('~'), 'ramBLe')), help='Base directory for running the experiments.')
    parser.add_argument('-s', '--scratch', metavar='DIR', type=str, default=realpath(join(expanduser('~'), 'scratch')), help='Scratch directory, visible to all the nodes.')
    parser.add_argument('-d', '--dataset', metavar='NAME', type=str, nargs='*', default=list(big_datasets.keys()), help='Datasets to be used.')
    parser.add_argument('-a', '--algorithm', metavar='NAME', type=str, nargs='*', default=all_algorithms, help='Algorithms to be used.')
    parser.add_argument('-p', '--process', metavar='P', type=int, nargs='*', default=all_processes, help='Processes to be used.')
    parser.add_argument('--ppn', metavar='PPN', type=int, nargs='*', default=[list(ppn_mappings.keys())[0]], help='Number of processes per node to be used.')
    parser.add_argument('-r', '--repeat', metavar='N', type=int, default=NUM_REPEATS, help='Number of times the experiments should be repeated.')
    parser.add_argument('--bnlearn', action='store_true', help='Flag for running bnlearn instead of our implementation.')
    parser.add_argument('--results', metavar = 'FILE', type=str, default='results_%s' % os.environ.get('PBS_JOBID', 0))
    args = parser.parse_args()
    return args


def get_executable_configurations(basedir, datasets, algorithms, bnlearn):
    boolean_args = ['-c', '-v', '-i']
    default_ramble_args = ['-r', '--warmup', '--hostnames']
    executable = join(basedir, 'ramble') if not bnlearn else join(basedir, 'scripts/ramble_bnlearn.R')
    configurations = []
    for name, algorithm in product(datasets, algorithms):
        dataset_args = all_datasets[name]
        script_args = [executable, '-d']
        script_args.append('-a %s -f %s' % (algorithm, join(basedir, dataset_args[0])))
        script_args.append('-n %d -m %d -s \'%s\'' % tuple(dataset_args[1:4]))
        script_args.extend(b for i, b in enumerate(boolean_args) if dataset_args[4 + i])
        if not bnlearn:
            script_args.extend(default_ramble_args)
        configurations.append((name, algorithm, ' '.join(script_args)))
    return configurations


def get_hostfile(scratch, ppn):
    nodefile = os.environ['PBS_NODEFILE']
    seen = set()
    hosts = []
    with open(nodefile, 'r') as nf:
        for n in nf.readlines():
            if n not in seen:
                hosts.append(n.strip() + ':%d' % ppn)
            seen.add(n)
    with NamedTemporaryFile(mode='w', suffix='.hosts', dir=scratch, delete=False) as hf:
        hf.write('\n'.join(hosts) + '\n')
    return hf.name


def get_mpi_configurations(scratch, processes, ppns):
    default_mpi_args = ['-env MV2_SHOW_CPU_BINDING 1']
    configurations = []
    ppn_hostfiles = dict((ppn, get_hostfile(scratch, ppn)) for ppn in ppns)
    for p, ppn in product(processes, ppns):
        mpi_args = ['mpirun -np %d -hostfile %s -env MV2_CPU_MAPPING %s' % (p, ppn_hostfiles[ppn], ppn_mappings[ppn])]
        mpi_args.extend(default_mpi_args)
        configurations.append((p, ' '.join(mpi_args)))
    return configurations


def get_runtime(action, output, required=True):
    import re

    float_pattern = r'((?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)'
    pattern = 'Time taken in %s: %s' % (action, float_pattern)
    match = re.search(pattern, output)
    if required:
        return float(match.group(1))
    else:
        return float(match.group(1) if match is not None else 0)

def parse_runtimes(output):
    # optional runtimes
    warmup = get_runtime('warming up MPI', output, required=False)
    redistributing = get_runtime('redistributing', output, required=False)
    blankets = get_runtime('getting the blankets', output, required=False)
    symmetry = get_runtime('symmetry correcting the blankets', output, required=False)
    sync = get_runtime('synchronizing the blankets', output, required=False)
    neighbors = get_runtime('getting the neighbors', output, required=False)
    direction = get_runtime('directing the edges', output, required=False)
    gsquare = get_runtime('G-square computations', output, required=False)
    # required runtimes
    reading = get_runtime('reading the file', output, required=True)
    network = get_runtime('getting the network', output, required=True)
    writing = get_runtime('writing the network', output, required=True)
    return [warmup, reading, redistributing, blankets, symmetry, sync, neighbors, direction, gsquare, network, writing]


def run_experiment(basedir, scratch, config, repeat, bnlearn):
    import subprocess

    runtimes = []
    dotfile = join(scratch, '%s_%s' % (config[0], config[1]))
    if bnlearn:
        dotfile += '_bnlearn'
    dotfile += '.dot'
    outfile = dotfile if not os.path.exists(dotfile) else NamedTemporaryFile(suffix='.dot', dir=scratch, delete=False).name
    for r in range(repeat):
        arguments = config[-1] + ' -o %s' % outfile
        print(arguments)
        sys.stdout.flush()
        output = subprocess.check_output(arguments, shell=True).decode('utf-8')
        print(output)
        sys.stdout.flush()
        if not bnlearn:
            print('Comparing generated file %s with %s' % (outfile, dotfile))
            sys.stdout.flush()
            subprocess.check_call(' '.join([join(basedir, 'scripts', 'compare_dot'), dotfile, outfile, '-d']), shell=True)
        runtimes.append(parse_runtimes(output))
    return runtimes


def main():
    '''
    Main function.
    '''
    args = parse_args()
    all_configs = get_executable_configurations(args.basedir, args.dataset, args.algorithm, args.bnlearn)
    if not args.bnlearn:
        mpi_configs = get_mpi_configurations(args.scratch, args.process, args.ppn)
        all_configs = list((executable[0], executable[1], mpi[0], mpi[-1] + ' ' + executable[-1]) for executable, mpi in product(all_configs, mpi_configs))
    with open(args.results, 'w') as results:
        results.write('# warmup,reading,redistributing,blankets,symmetry,sync,neighbors,direction,gsquare,network,writing\n')
        for config in all_configs:
            if not args.bnlearn:
                results.write('# our runtime for dataset=%s using algorithm=%s on processors=%d \n' % tuple(config[:3]))
            else:
                results.write('# bnlearn\'s runtime for dataset=%s using algorithm=%s \n' % tuple(config[:2]))
            for rt in run_experiment(args.basedir, args.scratch, config, args.repeat, args.bnlearn):
                results.write(','.join(str(t) for t in rt) + '\n')
                results.flush()


if __name__ == '__main__':
    main()
