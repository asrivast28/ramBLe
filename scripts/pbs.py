#!/usr/bin/env python

##
# @file pbs.py
# @brief Script for submitting PBS jobs.
# @author Ankit Srivastava <asrivast@gatech.edu>


def parse_args():
    '''
    Parse command line arguments.
    '''
    import argparse
    import os
    from os.path import abspath, dirname

    parser = argparse.ArgumentParser(description='Submit custom script using PBS')
    parser.add_argument('-s', '--script', type=str, metavar='FILE', required=True, help='Name of the script to be submitted.')
    parser.add_argument('-n', '--name', type=str, default='benchmark-csl', metavar='JOB', help='Name of the job.')
    parser.add_argument('-l', '--time', type=str, default='12:00:00', metavar='HH:MM:SS', help='Duration of the job.')
    parser.add_argument('-N', '--nodes', type=int, default=1, metavar='NODES', help='Number of required nodes.')
    parser.add_argument('-q', '--queue', type=str, default='hive', metavar='NAME', help='Name of the queue.')
    parser.add_argument('-o', '--output', type=str, metavar='FILE', help='Name of the output file.')
    parser.add_argument('-d', '--depend', type=str, metavar='JOBID', help='ID of the job on which this job depends.')
    args = parser.parse_args()
    return args


def create_submission_script(args):
    '''
    Create preamble for the script.
    '''
    from shutil import copymode
    from tempfile import NamedTemporaryFile

    preamble_format = \
'''\
#PBS -N %s                 # job name
#PBS -l walltime=%s        # duration of the job
#PBS -l nodes=%d:ppn=%d    # number of nodes and cores per node
#PBS -q %s                 # queue name (where job is submitted)
#PBS -j oe                 # combine output and error messages into a single file
#PBS -o %s                 # output file name
'''
    PROCS_PER_NODE = 24
    output = args.output if args.output is not None else args.name + '.out'
    preamble = preamble_format % (args.name, args.time, args.nodes, PROCS_PER_NODE, args.queue, output)
    if args.depend:
        preamble += \
'#PBS -W depend=afterok:%s  # job ID of the job on which this job depends' % args.depend
    preamble += '\n'

    with NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as pbs:
        pbs.write(preamble)
        with open(args.script, 'r') as orig:
            pbs.write(orig.read())
    copymode(args.script, pbs.name)
    return pbs.name


def submit_job(script):
    '''
    Submit a PBS job.
    '''
    import subprocess
    subprocess.check_call(['qsub', script])


def main():
    '''
    Main function.
    '''
    args = parse_args()
    script = create_submission_script(args)
    # print(script)
    submit_job(script)

if __name__ == '__main__':
    main()
