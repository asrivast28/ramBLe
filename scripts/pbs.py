#!/usr/bin/env python

##
# @file pbs.py
# @brief Script for submitting PBS jobs.
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


def parse_args():
    '''
    Parse command line arguments.
    '''
    import argparse
    import multiprocessing
    import os
    from os.path import abspath, dirname

    parser = argparse.ArgumentParser(description='Submit custom script using PBS')
    parser.add_argument('-s', '--script', type=str, metavar='FILE', required=True, help='Name of the script to be submitted.')
    parser.add_argument('-n', '--name', type=str, default='benchmark-ramble', metavar='JOB', help='Name of the job.')
    parser.add_argument('-t', '--time', type=str, default='12:00:00', metavar='HH:MM:SS', help='Duration of the job.')
    parser.add_argument('-N', '--nodes', type=int, default=1, metavar='NODES', help='Number of required nodes.')
    parser.add_argument('-p', '--procs', type=int, default=multiprocessing.cpu_count(), metavar='PROCS', help='Number of processes per node.')
    parser.add_argument('-q', '--queue', type=str, default='hive', metavar='NAME', help='Name of the queue.')
    parser.add_argument('-o', '--output', type=str, metavar='FILE', help='Name of the output file.')
    parser.add_argument('-d', '--depend', type=str, metavar='JOBID', help='ID of the job on which this job depends.')
    parser.add_argument('-a', '--after', type=str, metavar='HHMM', help='Schedule the job after the given time.')
    parser.add_argument('-v', '--variables', type=str, nargs='*', default=[], metavar='VAR=VALUE', help='Variables to pass to the run script.')
    args = parser.parse_args()
    return args


def create_submission_script(args):
    '''
    Create preamble for the script.
    '''
    import os.path
    from shutil import copymode
    from tempfile import NamedTemporaryFile

    preamble_lines = [
        '#PBS -N %s\t# job name',
        '#PBS -l walltime=%s\t# duration of the job',
        '#PBS -l nodes=%d:ppn=%d\t# number of nodes and cores per node',
        '#PBS -q %s\t# queue name (where job is submitted)',
        '#PBS -j oe\t# combine output and error messages into a single file',
        '#PBS -o %s\t# output file name',
    ]
    output = args.output if args.output is not None else args.name + '.out'
    if os.path.exists(output):
        raise RuntimeError('Output file %s already exists' % output)
    if args.depend:
        preamble_lines.append('#PBS -W depend=afterok:%s\t# job ID of the job on which this job depends' % args.depend)
    if args.after:
        preamble_lines.append('#PBS -a %s\t# delay executing the job until the given time and date' % args.after)
    preamble = '\n'.join(preamble_lines) % (args.name, args.time, args.nodes, args.procs, args.queue, output) + '\n'

    with NamedTemporaryFile(mode='w', suffix=os.path.splitext(args.script)[1], delete=False) as pbs:
        with open(args.script, 'r') as orig:
            lines = orig.readlines()
            lines.insert(1 if lines[0].startswith(r'#!') else 0, preamble)
            pbs.write(''.join(lines))
    copymode(args.script, pbs.name)
    return pbs.name


def submit_job(script, variables):
    '''
    Submit a PBS job.
    '''
    import subprocess
    command = ['qsub', script]
    if variables:
        command.append('-v ' + ','.join(variables))
    subprocess.check_call(command)


def main():
    '''
    Main function.
    '''
    args = parse_args()
    script = create_submission_script(args)
    # print(script)
    submit_job(script, args.variables)

if __name__ == '__main__':
    main()
