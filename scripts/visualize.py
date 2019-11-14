#!/usr/bin/env python

import pydot

from compare import read


def parse_args():
    '''
    Parse command line arguments.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='Visualize a network')
    parser.add_argument('file', type=str, metavar='NAME', help='Name of the file.')
    parser.add_argument('format', type=str, nargs='?', default='dot', metavar='FORMAT', help='Format of the file.')
    args = parser.parse_args()
    return args


def visualize(graph):
    '''
    Visualize a given pydot graph.
    '''
    import subprocess
    import tempfile

    with tempfile.NamedTemporaryFile(suffix='.png') as image:
        graph.write_png(image.name)
        subprocess.check_call(['display', image.name])


def main():
    '''
    Main function.
    '''
    args = parse_args()
    graph = read(args.file, args.format)
    visualize(graph)


if __name__ == '__main__':
    main()
