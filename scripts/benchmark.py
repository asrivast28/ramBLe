#!/usr/bin/env python

##
# @file benchmark.py
# @brief Script for benchmarking structure learning.
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

datasets = [
    ('child',      20),
    ('insurance',  27),
    ('mildew',     35),
    ('alarm',      37),
    ('barley',     46),
    ('hailfinder', 56),
]

directories = [
    ('1K',   10 ** 3),
    ('10K',  10 ** 4),
    ('100K', 10 ** 5),
    ('1M',   10 ** 6),
    # ('10M',  10 ** 7),
]


def parse_args():
    '''
    Parse command line arguments.
    '''
    import argparse
    import os
    from os.path import abspath, dirname

    parser = argparse.ArgumentParser(description='Benchmark structure learning')
    parser.add_argument('-b', '--base', type=str, default=dirname(dirname(abspath(__file__))), metavar='DIR', help='Base directory for running the experiments.')
    parser.add_argument('-e', '--executable', type=str, metavar='NAME', required=True, help='Name of the executable to be used.')
    parser.add_argument('-r', '--repeat', type=int, default=3, metavar='N', help='Number of times the experiments should be repeated.')
    args = parser.parse_args()
    args.executable = abspath(args.executable)
    return args


def parse_results(output):
    '''
    Parse benchmark metrics from the given output.
    '''
    import re

    net = float(re.search('Time taken in getting the network: (\d+.\d+)', output).group(1))
    match = re.search('Time taken in G-square computations: (\d+.\d+)', output)
    gsq = float(match.group(1) if match is not None else 0)
    mem = int(re.search('Maximum resident set size \(kbytes\): (\d+)', output).group(1))
    return net, gsq, mem


def run_experiments(base, executable, repeat):
    '''
    Run experiments and print the benchmark metrics.
    '''
    import subprocess

    for d, m in directories:
        for s, n in datasets:
            results = []
            for r in range(repeat):
                command = '/usr/bin/time -v %s -n %d -m %d -f %s/data/%s/%s.csv -c -s \' \' -l' % (executable, n, m, base, d, s)
                output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True).decode('utf-8')
                results.append(parse_results(output))
            print('data/%s/%s.csv' % (d, s))
            print('Network =AVERAGE(%s)' % ','.join(str(r[0]) for r in results))
            print('G-square =AVERAGE(%s)' % ','.join(str(r[1]) for r in results))
            print('Memory =AVERAGE(%s) / 1024' % ','.join(str(r[2]) for r in results))
            print('\n')


def main():
    '''
    Main function.
    '''
    args = parse_args()
    run_experiments(args.base, args.executable, args.repeat)

if __name__ == '__main__':
    main()
