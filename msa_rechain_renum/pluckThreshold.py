#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import argparse
import itertools
import textwrap

import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')

from res_naming import *

def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)


def get_parser():
    parser = argparse.ArgumentParser(description='JSON input file'
            ' to select sequences from.')
    parser.add_argument('--threshold', '-t', type=float, required=True,
            help='Select sequences with'
            ' alignment scores of this value or more.')
    parser.add_argument('input', type=str,
            help='Name of settings file with a sequence archetype')
    return parser


AAs = "ACDEFGHIKLMNPQRSTVY"


def main(argv=None):
    aparser = get_parser()
    args = aparser.parse_args(argv)

    with open(args.input) as input:
        data = json.load(input)

    for fname, chains in data['sequences'].items():
        for chain in chains:
            if chain['align_score'] >= args.threshold:
                print('>' + fname + '|' + chain['chain'])
                for l in textwrap.wrap(chain['sequence'], 80):
                    print(l)

    return 0

if __name__ == '__main__':
    exit(main(sys.argv[1:]))
