#!/usr/bin/env python

##
# @file compare.py
# @brief Script for comparing two graphs.
# @author Ankit Srivastava <asrivast@gatech.edu>

import pydot
import warnings


def parse_args():
    '''
    Parse command line arguments.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='Discretize continuous data')
    parser.add_argument('-f', '--first', type=str, required=True, metavar='NAME', help='Name of the file which contains the first network.')
    parser.add_argument('--first_format', type=str, default='dot', metavar='FORMAT', help='Format of the file which contains the first network.')
    parser.add_argument('-s', '--second', type=str, required=True, metavar='NAME', help='Name of the file which contains the second network.')
    parser.add_argument('--second_format', type=str, default='dot', metavar='FORMAT', help='Format of the file which contains the second network.')
    args = parser.parse_args()
    return args


def read(name, file_format):
    '''
    Read the network from the given file.
    '''
    graph = None
    if file_format == 'dot':
        with open(name, 'rt') as f:
            s = f.read()
        s = s.replace('"', '')
        graph = pydot.graph_from_dot_data(s)
        if graph is not None:
            graph = graph[0]
    elif file_format == 'el':
        with open(name, 'r') as f:
            edges = [l.split() for l in f.readlines()]
        graph = pydot.graph_from_edges(edges)
    else:
        raise RuntimeError('Unknown file format', file_format)
    print('Type of %s is %s' % (name, graph.get_graph_type()))
    return graph


def compare(first, second):
    '''
    Compare the first network with the second network.
    '''
    tp = 0
    fn = 0
    for fe in first.get_edges():
        edge = tuple([fe.get_source(), fe.get_destination()])
        se = second.get_edge(edge)
        if len(se) == 1:
            tp += 1
        else:
            fn += 1
    fp = len(second.get_edges()) - tp
    return (tp, fp, fn)


def main():
    '''
    Main function.
    '''
    args = parse_args()
    first = read(args.first, args.first_format)
    second = read(args.second, args.second_format)
    if first.get_graph_type() != second.get_graph_type():
        warnings.warn('Comparing graphs of different types', RuntimeWarning, stacklevel=2)
    tp, fp, fn = compare(first, second)
    print('\nComparison results')
    print('# of edges found only in %s: %d' % (args.first, fn))
    print('# of edges found only in %s: %d' % (args.second, fp))
    print('# of edges common to both: %d' % tp)


if __name__ == '__main__':
    main()
