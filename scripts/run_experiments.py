#!/usr/bin/env python

from collections import OrderedDict
from itertools import product
import os
import os.path
from os.path import join
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
    ])

NUM_REPEATS = 5


def parse_args():
    '''
    Parse command line arguments.
    '''
    import argparse
    from os.path import expanduser, realpath

    parser = argparse.ArgumentParser(description='Run scaling experiments')
    parser.add_argument('-b', '--basedir', metavar='DIR', type=str, default=realpath(join(expanduser('~'), 'discover-mb')), help='Base directory for running the experiments.')
    parser.add_argument('-s', '--scratch', metavar='DIR', type=str, default=realpath(join(expanduser('~'), 'scratch')), help='Scratch directory, visible to all the nodes.')
    parser.add_argument('-d', '--dataset', metavar='NAME', type=str, nargs='*', default=list(big_datasets.keys()), help='Datasets to be used.')
    parser.add_argument('-a', '--algorithm', metavar='NAME', type=str, nargs='*', default=all_algorithms, help='Algorithms to be used.')
    parser.add_argument('-p', '--process', metavar='P', type=int, nargs='*', default=all_processes, help='Processes to be used.')
    parser.add_argument('--ppn', metavar='PPN', type=int, nargs='*', default=[list(ppn_mappings.keys())[0]], help='Number of processes per node to be used.')
    parser.add_argument('-r', '--repeat', metavar='N', type=int, default=NUM_REPEATS, help='Number of times the experiments should be repeated.')
    parser.add_argument('--bnlearn', action='store_true', help='Flag for running bnlearn instead of our implementation.')
    parser.add_argument('--results', metavar = 'FILE', type=str, default='results_%s' % os.environ['PBS_JOBID'])
    args = parser.parse_args()
    return args


def get_executable_configurations(basedir, datasets, algorithms, bnlearn):
    boolean_args = ['-c', '-v', '-i']
    default_csl_args = ['-r', '-d', '--warmup', '--hostnames']
    executable = join(basedir, 'csl') if not bnlearn else join(basedir, 'scripts/csl_bnlearn.R')
    configurations = []
    for name, algorithm in product(datasets, algorithms):
        dataset_args = all_datasets[name]
        script_args = [executable]
        script_args.append('-a %s -f %s' % (algorithm, join(basedir, dataset_args[0])))
        script_args.append('-n %d -m %d -s \'%s\'' % tuple(dataset_args[1:4]))
        script_args.extend(b for i, b in enumerate(boolean_args) if dataset_args[4 + i])
        if not bnlearn:
            script_args.extend(default_csl_args)
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


def parse_runtimes(output):
    import re

    float_pattern = '((?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)'
    warmup = re.search('Time taken in warming up MPI: ' + float_pattern, output)
    warmup = float(0 if warmup is None else warmup.group(1))
    reading = float(re.search('Time taken in reading the file: ' + float_pattern, output).group(1))
    blankets = re.search('Time taken in getting the blankets: ' + float_pattern, output)
    blankets = float(0 if blankets is None else blankets.group(1))
    symmetry = re.search('Time taken in symmetry correcting the blankets: ' + float_pattern, output)
    symmetry = float(0 if symmetry is None else symmetry.group(1))
    neighbors = re.search('Time taken in getting the neighbors: ' + float_pattern, output)
    neighbors = float(0 if neighbors is None else neighbors.group(1))
    direction = float(re.search('Time taken in directing the edges: ' + float_pattern, output).group(1))
    gsquare = float(re.search('Time taken in G-square computations: ' + float_pattern, output).group(1))
    network = float(re.search('Time taken in getting the network: ' + float_pattern, output).group(1))
    writing = float(re.search('Time taken in writing the network: ' + float_pattern, output).group(1))
    return [warmup, reading, blankets, symmetry, neighbors, direction, gsquare, network, writing]


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
        output = subprocess.check_output(arguments, shell=True).decode('utf-8')
        print(output)
        if not bnlearn:
            print('Comparing generated file %s with %s' % (outfile, dotfile))
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
        results.write('# warmup,reading,blankets,symmetry,neighbors,direction,gsquare,network,writing\n')
        for config in all_configs:
            results.write('# Runtimes for dataset=%s using algorithm=%s on processors=%d \n' % tuple(config[:3]))
            for rt in run_experiment(args.basedir, args.scratch, config, args.repeat, args.bnlearn):
                results.write(','.join(str(t) for t in rt) + '\n')
                results.flush()


if __name__ == '__main__':
    main()
