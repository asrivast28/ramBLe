#!/usr/bin/env python

import subprocess
import sys
import tempfile


def visualize(dotFile):
    with tempfile.NamedTemporaryFile(suffix='.png') as image:
        subprocess.check_call(['dot', '-Tpng', '-Gdpi=150', dotFile], stdout=image)
        subprocess.check_call(['display', image.name])

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: ./visualize_dot.py <dot file name>')
    else:
        visualize(sys.argv[1])
