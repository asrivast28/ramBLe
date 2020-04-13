#!/usr/bin/env python

##
# @file visualize.py
# @brief Script for visualizing graphs written as dot files.
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
