#!/usr/bin/env python

##
# @file discretize.py
# @brief Script for discretizing continuous datasets.
# @author Ankit Srivastava <asrivast@gatech.edu>

import pandas


def parse_args():
    '''
    Parse command line arguments.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='Discretize continuous data')
    parser.add_argument('-f', '--file', type=str, metavar='NAME', help='Name of the input file.')
    parser.add_argument('-s', '--separator', type=str, metavar='DELIM', help='Delimiting character in the file.')
    parser.add_argument('-c', '--colobs', action='store_true', help='The file contains observations in columns')
    parser.add_argument('-v', '--varnames', action='store_true', help='The file contains variable names.')
    parser.add_argument('-i', '--indices', action='store_true', help='The file contains observation indices.')
    parser.add_argument('-o', '--out', type=str, metavar='NAME', help='Name of the output file.')
    args = parser.parse_args()
    return args


def read(name, sep, colobs, varnames, indices):
    '''
    Read the dataset from the given CSV file.
    '''
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
    dataset = pandas.read_csv(name, sep=sep, header=header, index_col=index)
    if colobs:
        dataset = dataset.T
    return dataset


def discretize_column(column, intervals):
    '''
    Discretize a column using the given interval edges.
    '''
    # Lower bound the given intervals
    intervals.insert(0, column.min())
    # Upper bound the given intervals
    intervals.insert(len(intervals), column.max())
    # Remove duplicates and sort the interval edges
    bins = sorted(set(intervals))
    cut = pandas.cut(column, bins, labels=False, include_lowest=True)
    # if (cut.isnull().sum() > 0):
        # raise RuntimeError('Found null values in the discretized observations for %s' % column.name)
    return cut


def create_intervals(control, multiples):
    '''
    Create interval edges by multiplying a control
    with the given multiples.
    '''
    return [control*m for m in multiples]


def discretize(dataset):
    '''
    Discretize a given dataset as per Friedman et al. (2000)
    '''
    control = dataset.mean(axis=0)
    upper = 2 ** 0.5
    lower = 1 / upper
    multiples = [lower, upper]
    discretized = dataset.apply(lambda c: discretize_column(c, create_intervals(c.mean(), multiples)), axis=0)
    return discretized


def write(dataset, name, sep, colobs):
    '''
    Write the dataset as a CSV file.
    '''
    if colobs:
        dataset = dataset.T
    dataset.to_csv(name, sep=sep)


def main():
    '''
    Main function.
    '''
    args = parse_args()

    dataset = read(args.file, args.separator, args.colobs, args.varnames, args.indices)
    discretized = discretize(dataset)
    write(discretized, args.out, args.separator, args.colobs)


if __name__ == '__main__':
    main()
