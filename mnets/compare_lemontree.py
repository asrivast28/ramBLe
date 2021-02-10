#!/usr/bin/env python

##
# @file compare_lemontree.py
# @brief Script for comparing output files for lemon-tree algorithm.
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
from os.path import exists, join

TOLERANCE = 1e-6

def parse_args():
    '''
    Parse command line arguments.
    '''
    global TOLERANCE
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Compare our output files with those from lemon-tree')
    parser.add_argument('first', metavar='DIR', help='Directory containing first set of output files')
    parser.add_argument('second', metavar='DIR', help='Directory containing second set of output files')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--tolerance', metavar='TOL', type=float, help='Tolerance for float comparison')
    args = parser.parse_args()

    if args.tolerance is not None:
        TOLERANCE = args.tolerance
    return args


def compare_configs(task, firstconfigs, secondconfigs):
    for name, firstval in firstconfigs.items():
        secondval = secondconfigs.pop(name, None)
        if secondval is None:
            print('WARNING: Second configs for task %s do not contain the key %s' % (task, name))
        elif firstval != secondval:
            print('WARNING: Mismatch in configs for task %s, key %s (%s != %s)' % (task, name, firstval, secondval))
    for name in secondconfigs.keys():
        print('WARNING: First configs for task %s do not contain key %s' % (task, name))


def read_task_outputs(firstdir, seconddir):
    '''
    Read configurations for all the tasks from the JSON file.
    '''
    from collections import OrderedDict
    from json import load

    with open(join(firstdir, 'configs.json'), 'r') as fc, open(join(seconddir, 'configs.json'), 'r') as sc:
        firstconfigs = load(fc, object_pairs_hook=OrderedDict)
        secondconfigs = load(sc, object_pairs_hook=OrderedDict)
    task_outputs = OrderedDict()
    for first, second in zip(firstconfigs.items(), secondconfigs.items()):
        if first[0] == second[0]:
            task = first[0]
            firstout = first[1].pop('output_file')
            firstout = join(firstdir, firstout) if firstout else None
            secondout = second[1].pop('output_file')
            secondout = join(seconddir, secondout) if secondout else None
            task_outputs[task] = tuple([firstout, secondout])
            compare_configs(task, first[1], second[1])
        else:
            print('WARNING: Mismatch in tasks in the config file (%s != %s)' % (first[0], second[0]))
    return task_outputs


def compare_clusters(firstfile, secondfile):
    match = True
    with open(firstfile, 'r') as ff, open(secondfile, 'r') as sf:
        firstclusters = dict(line.split() for line in ff)
        secondclusters = dict(line.split() for line in sf)
    for var, cluster in firstclusters.items():
        sc = secondclusters.pop(var, None)
        if sc is None:
            print('ERROR: %s not found in %s' % (var, secondfile))
        elif cluster != sc:
            print('ERROR: Mismatch in cluster of %s (%s != %s)' % (var, cluster, sc))
            match = False
    if len(secondclusters) > 0:
        match = False
        print('ERROR: %s is missing some clusters in %s' % (firstfile, secondfile))
    return match


def compare_ganesh(firstout, secondout, verbose):
    with open(firstout, 'r') as fo, open(secondout, 'r') as so:
        firstfiles = list(f.strip() for f in fo)
        secondfiles = list(f.strip() for f in so)
    if len(firstfiles) != len(secondfiles):
        print('ERROR: Mismatch in number of GaneSH runs (%d != %d)' % (len(firstfiles), len(secondfiles)))
        return False
    match = True
    for first, second in zip(firstfiles, secondfiles):
        if verbose:
            print('Comparing GaneSH clusters in %s and %s' % (first, second))
        match &= compare_clusters(first, second)
    return match


def compare_consensus(firstout, secondout, verbose):
    if verbose:
        print('Comparing consensus clusters in %s and %s' % (firstout, secondout))
    return compare_clusters(firstout, secondout)


def read_regulators(filename):
    regulators = {}
    with open(filename, 'r') as rf:
        for line in rf:
            reg, module, score = line.split()
            regulators[tuple([reg, module])] = float(score)
    return regulators


def compare_regulators(firstfile, secondfile):
    match = True
    firstregulators = read_regulators(firstfile)
    secondregulators = read_regulators(secondfile)
    for vc, score in firstregulators.items():
        ls = secondregulators.pop(vc, None)
        if ls is None:
            print('ERROR: %s not a regulator of module %s in %s' % (vc[0], vc[1], secondfile))
        elif abs(score - ls) >= TOLERANCE:
            print('ERROR: Mismatch in score of regulator %s of module %s (%f != %f)' % (vc[0], vc[1], score, ls))
            match = False
    if len(secondregulators) > 0:
        match = False
        print('ERROR: %s is missing some regulator assignments in %s' % (firstfile, secondfile))
    return match


def read_xml_tree(filename):
    import gzip
    import xml.etree.ElementTree as ET

    decompressed = gzip.open(filename, 'r')
    return ET.parse(decompressed)


def compare_xml_trees(e1, e2, path='*'):
    if e1.tag != e2.tag:
        print('ERROR: Mismatch in tag at path %s (%s != %s)' % (path, e1.tag, e2.tag))
        return False
    if len(e1.attrib) != len(e2.attrib):
        print('ERROR: Mismatch in number of attributes at path %s (%d != %d)' % (path, len(e1.attrib), len(e2.attrib)))
        return False
    match = True
    for a1, a2 in zip(e1.attrib.items(), e2.attrib.items()):
        if a1[0] == a2[0]:
            try:
                if abs(float(a1[1]) - float(a2[1])) >= TOLERANCE:
                    print('ERROR: Mismatch in value of %s at path %s (%s != %s)' % (a1[0], path, a1[1], a2[1]))
                    match = False
            except ValueError:
                if a1[1] != a2[1]:
                    print('ERROR: Mismatch in value of %s at path %s (%s != %s)' % (a1[0], path, a1[1], a2[1]))
                    match = False
        else:
            print('ERROR: Mismatch in attribute keys at path %s (%s != %s)' % (path, a1[0], a1[1]))
            match = False
    if not match:
        return False
    if len(e1) != len(e2):
        print('ERROR: Mismatch in number of children at path %s (%d != %d)' % (path, len(e1), len(e2)))
        return False
    return all([compare_xml_trees(c1, c2, '%s -> %s' % (path, e1.tag)) for c1, c2 in zip(e1, e2)])


def compare_modules(firstout, secondout, verbose):
    match = True
    for ext in ('all', 'top', 'random'):
        firstregs = '%s.%sreg.txt' % (firstout, ext)
        secondregs = '%s.%sreg.txt' % (secondout, ext)
        if verbose:
            print('Comparing regulators in %s and %s' % (firstregs, secondregs))
        match &= compare_regulators(firstregs, secondregs)
    firstxml = '%s.xml.gz' % firstout
    secondxml = '%s.xml.gz' % secondout
    if verbose:
        print('Comparing module networks in %s and %s' % (firstxml, secondxml))
    firstmodules = read_xml_tree(firstxml)
    secondmodules = read_xml_tree(secondxml)
    match &= compare_xml_trees(firstmodules.getroot(), secondmodules.getroot())
    return match


def main():
    '''
    Main function.
    '''
    args = parse_args()
    task_outputs = read_task_outputs(args.first, args.second)
    match = True
    for task, outputs in task_outputs.items():
        print('Comparisons for task %s' % task)
        first, second = outputs
        if first is None or second is None:
            print('WARNING: Skipping comparison for task %s because output files are missing\n' % task)
            continue
        if task == 'ganesh':
            match &= compare_ganesh(first, second, args.verbose)
        elif task == 'tight_clusters':
            match &= compare_consensus(first, second, args.verbose)
        elif task == 'regulators':
            match &= compare_modules(first, second, args.verbose)
        print('Done comparisons for task %s\n' % task)
    assert match, 'Mismatches found in the output files'


if __name__ == '__main__':
    main()
