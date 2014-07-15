#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import argparse
import itertools
import textwrap
import sqlite3

import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')

from res_naming import *

def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)


def get_parser():
    parser = argparse.ArgumentParser(description='JSON input file'
            ' to select sequences from.')
    parser.add_argument('--threshold', '-t', type=float,
            help='Select sequences with alignment scores of this value or more.')
    parser.add_argument('input', type=str,
            help='Name of settings file with a sequence archetype')
    parser.add_argument('--dbase_name', '-d', type=str, required=True,
            help='Name of sqlite3 database.')
    parser.add_argument('--table_name', '-a',  type=str, default='ubiquitins',
            help='Name of table to insert ubiquitins into.')
    return parser


AAs = "ACDEFGHIKLMNPQRSTVY"


def main(argv=None):

    aparser = get_parser()
    args = aparser.parse_args(argv)

    with open(args.input) as input:
        data = json.load(input)

    ubqs_pdbs_chains = []
    for fname, chains in data['sequences'].items():
        for chain in chains:
            threshold = args.threshold
            if threshold and 'align_score' in chain and chain['align_score'] < args.threshold:
                continue
            pdbname = chain['protein_name']
            ubqs_pdbs_chains.append( (pdbname, chain['chain']) )

    # print(ubqs_pdbs_chains)
    # exit(1)

    conn = sqlite3.connect(args.dbase_name)
    c = conn.cursor()


    c.execute('CREATE TABLE IF NOT EXISTS '  + args.table_name + '(pdb_code TEXT NOT NULL, chain_id TEXT NOT NULL)')
    for pdbname, chainname in ubqs_pdbs_chains:
        stmt = 'INSERT INTO ' + args.table_name + ' VALUES ("' + pdbname + '", "' + chainname + '")'
        c.execute(stmt)
    conn.commit()
    conn.close()


    return 0

if __name__ == '__main__':
    exit(main(sys.argv[1:]))
