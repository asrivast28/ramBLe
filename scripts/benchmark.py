#!/usr/bin/env python

##
# @file benchmark.py
# @brief Script for benchmarking datasets 
# @author Ankit Srivastava <asrivast@gatech.edu>

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
    import argparse
    import os.path

    parser = argparse.ArgumentParser(description='Benchmark structure learning')
    parser.add_argument('-b', '--base', type=str, default=os.path.abspath(os.pardir), metavar='DIR', help='Base directory for running the experiments.')
    parser.add_argument('-e', '--executable', type=str, metavar='NAME', required=True, help='Name of the executable to be used.')
    parser.add_argument('-r', '--repeat', type=int, default=3, metavar='N', help='Number of times the experiments should be repeated.')
    args = parser.parse_args()
    args.executable = os.path.abspath(args.executable)
    return args

def read_results(output):
    import re

    net = float(re.search('Time taken in getting the network: (\d+.\d+) sec', output).group(1))
    match = re.search('Time taken in calls to SABNAtk: (\d+.\d+) sec', output)
    sab = float(match.group(1) if match is not None else 0)
    mem = int(re.search('Maximum resident set size \(kbytes\): (\d+)', output).group(1))
    return net, sab, mem

def run_experiments(base, executable, repeat):
    import os
    import subprocess
    import tempfile

    for d, m in directories:
        for s, n in datasets:
            results = []
            for r in range(repeat):
                outfile = tempfile.NamedTemporaryFile(suffix='.dot', delete=False).name
                command = '/usr/bin/time -v %s -n %d -m %d -f %s/data/%s/%s.csv -c -s \' \' -o %s -w' % (executable, n, m, base, d, s, outfile)
                output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True).decode('utf-8')
                os.remove(outfile)
                results.append(read_results(output))
            print('data/%s/%s.csv' % (d, s))
            print('Network =AVERAGE(%s)' % ','.join(str(r[0]) for r in results))
            print('SABNAtk =AVERAGE(%s)' % ','.join(str(r[1]) for r in results))
            print('Memory =AVERAGE(%s) / 1024' % ','.join(str(r[2]) for r in results))
            print('\n')

def main():
    args = parse_args()
    run_experiments(args.base, args.executable, args.repeat)

if __name__ == '__main__':
    main()
