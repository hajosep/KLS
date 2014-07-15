#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import string
import argparse
import gzip
from collections import OrderedDict
import itertools
# from Bio import pairwise2
import pairwise2
import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')

from res_naming import *

def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)


def proteinSequencesFromFasta(fasta):
    seqs = {}
    it = iter(fasta)
    l = next(it)
    working_seq = []
    chain = l[6]
    try:
        while True:
            if not l.startswith('>') or string.count(l, '|') != 3 or \
                    len(l.partition('|')[0]) != 7:
                raise IOError('Malformed FASTA file line: ' + l)
            chain = l[6]
            if chain in seqs:
                raise IOError('Wtf? Redefinion of chain!')
            working_seq = []
            l = next(it)
            while not l.startswith('>'):
                working_seq += l.strip()
                l = next(it)
            seqs[chain] = working_seq
    except StopIteration:
        seqs[chain] = working_seq
        for k, v in seqs.items():
            seqs[k] = ''.join(v)
        return seqs


def threeLetterToOne(threeLetter):
    if threeLetter in longer_names:
        return longer_names[threeLetter]
    elif threeLetter in modres:
        unmod = longer_names[modres[threeLetter]]
        warning('Interpretting ' + threeLetter + ' as a modified ' + unmod)
        return unmod
    raise IOError('Bad three letter res record:', threeLetter)


def parseSEQRESline(l):
    tokens = l.strip().split()
    assert(tokens[0] == 'SEQRES')
    # Test if the
    try:
        int(tokens[1])
        int(tokens[3])
    except ValueError:
        print('Bad SEQRES record:', l.strip(), file=sys.stderr)
        raise
    chain = tokens[2]
    seq = [threeLetterToOne(t) for t in tokens[4:]]
    return chain, seq


def seqFromSEQRESRecs(it):
    line = next(it)
    if not line.startswith('SEQRES'):
        return {}
    ch, seq = parseSEQRESline(line)
    seqs = seqFromSEQRESRecs(it)
    if ch in seqs:
        seqs[ch] = seq + seqs[ch]
    else:
        seqs[ch] = seq
    return seqs


def parseATOMline(l):
    tokens = l.strip().split()
    assert(tokens[0] == 'ATOM')
    # Test if the
    try:
        int(tokens[1])
    except ValueError:
        raise IOError('Bad ATOM record:', l.strip())
    chain = tokens[4]
    if chain not in seqs:
        seqs[chain] = []
    to
    resName3 = l[17:20]
    seqs[chain].append(threeLetterToOne(resName3))

def parseATOMline(l):
    thisChain = l[21]
    thisRes = threeLetterToOne(l[17:20])
    return thisChain, thisRes

def seqFromATOMRecs(it):
    seqs = {}
    line = next(it)
    while line.startswith('ATOM  ') or line.startswith('TER') or \
            line.startswith('ANISOU') or line.startswith('HETATM'):
        if line.startswith('ATOM  ') and line[13:15] == 'CA':
            ch, letter = parseATOMline(line)
            if ch in seqs:
                seqs[ch] += letter
            else:
                seqs[ch] = letter
        try:
            line = next(it)
        except StopIteration:
            pass
    return seqs


def proteinSequencesFromPDB(pdb):
    lead_it, = itertools.tee(pdb, 1)
    try:
        while True:
            trail_it = lead_it.__copy__()
            l = next(lead_it)
            if l.startswith('SEQRES'):
                return seqFromSEQRESRecs(trail_it)
            elif l.startswith('ATOM  '):
                warning('No SEQRES records detected in ' + pdb.name +
                        ', falling back to ATOM records. Missing density'
                        ' cannot be accounted for.')
                return seqFromATOMRecs(trail_it)
    except StopIteration:
        raise ValueError('No sequence found in PDB ' + pdb.name)


def proteinSequencesFromFilename(fname):
    try:
        seqs = {}
        if fname.endswith('.pdb'):
            with open(fname) as pdb:
                print('Searching PDB file ' + fname + ' for sequence')
                seqs = proteinSequencesFromPDB(pdb)
        elif fname.endswith('.pdb.gz'):
            with gzip.open(fname) as pdb:
                print('Searching PDB file ' + fname + ' for sequence')
                seqs = proteinSequencesFromPDB(pdb)
        elif fname.endswith('.fasta') or fname.endswith('.fas') or \
                fname.endswith('.faa'):
            with open(fname) as fasta:
                print('Searching FASTA file ' + fname + ' for sequence')
                seqs = proteinSequencesFromFasta(fasta)
        elif fname.endswith('.fasta.gz') or fname.endswith('.fas.gz') or \
                fname.endswith('.faa.gz'):
            with gzip.open(fname) as fasta:
                print('Searching FASTA file ' + fname + ' for sequence')
                seqs = proteinSequencesFromFasta(fasta)
        else:
            raise ValueError('Unrecognized input format: "' + fname + '"')
        for k, v in seqs.items():
            # yuck, whatever, I'm tired
            if isinstance(v, list):
                seqs[k] = ''.join(v)
        print('Found', len(seqs), 'sequences.')
        return seqs
    except EnvironmentError:
        print('Unable to read from file: "' + fname + '"')


def get_parser():
    parser = argparse.ArgumentParser(description='Renumbering and rename'
            ' chains in a set of PDB files' 'accoring to a MSA.')
    parser.add_argument('--config', '-c', type=str,
            default='findChainsMatchingSeqSettings.txt', help='Name of'
            ' settings file with a sequence archetype')
    parser.add_argument('--list', '-l', type=str, help='List of pdb and fasta'
            ' files to search for matching sequences')
    parser.add_argument('--output', '-o', type=str, required=True,
            help='Name of file to output sequence matches')
    parser.add_argument('pdbs', nargs='*',
            help='As an alternative to list PDB files to search in a file,'
            ' list them on the command line.')
    return parser


AAs = "ACDEFGHIKLMNPQRSTVY"


def findMatches(seqs):
    pass


# Default settings
default_config = {
    'sequenceArchetype': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
    'matchScore': 1,
    'mismatchScore': -0.5,
    'gapScore': -0.25,
    'gapExtensionScore': -0.1,
}


def read_and_check_listfile(listfname):
    if not os.access(listfname, os.R_OK):
        raise ValueError('List file ' + args.list + ' cannot be read!')
    with open(listfname) as listfile:
        fnames = [l.strip() for l in listfile.readlines()]
        for f in fnames:
            if not os.access(f, os.R_OK):
                raise ValueError('File ' + f + ', listed in ' +
                        listfname + ', cannot be read!')
    return fnames


def main(argv=None):
    aparser = get_parser()
    args = aparser.parse_args(argv)

    # check that there are pdb inputs, list or command line

    # Check for args.list or command line files
    # Warn if both present
    # Check that files are readable
    if args.list and args.pdbs:
        warning('Arguments have been given for both command line and list file'
                ' input of filenames to inspect. Both will be used, but this'
                ' may not have been intented.')

    # TODO
    input_fnames = []
    if args.list:
        input_fnames = read_and_check_listfile(args.list)

    for f in args.pdbs:
        if not os.access(f, os.R_OK):
            raise IOError('Positional argument ' + f + ' is not a file that'
                    ' can be read')
    input_fnames.extend(args.pdbs)

    config = OrderedDict()
    if os.access(args.config, os.R_OK):
        config_fname = os.path.abspath(args.config)
        print('Reading pairwise alignment config from ' + config_fname)
        execfile(config_fname, config)
    else:
        config = default_config
        print('Using default configuration for pairwise alignment',
              file=sys.stderr)

    print('Pairwise alignment configuration:', file=sys.stderr)
    maxlen = max(len(k) for k in config.keys())
    for k, v in config.iteritems():
        print(' '*4, k, (maxlen - len(k))*' ' + '=', v, file=sys.stderr)

    # Search every file for relevant chains

    def align_scorer(candidate):
        return max( a[2] for a in pairwise2.align.localms(candidate,
            config['sequenceArchetype'], config['matchScore'],
            config['mismatchScore'], config['gapScore'],
            config['gapExtensionScore'])
        )

    align_scores = {}

    # seqs = {}
    for f in input_fnames:
        seqs_from_file = proteinSequencesFromFilename(f)
        # print('seqs_from_file:', seqs_from_file)
        align_scores[f] = []
        for c, s in seqs_from_file.items():
            score = align_scorer(s)
            protein_name = os.path.splitext(os.path.basename(f))[0].upper()
            align_scores[f].append({'protein_name': protein_name, 'chain': c, 'sequence': s, 'align_score': score})
            # print(' '*3, 'chain', c, 'align score = ', score)

    print('Writing data to ', args.output)
    data = {'args': args.__dict__, 'config': config, 'sequences': align_scores}
    with open(args.output, 'w') as fp:
        json.dump(data, fp, sort_keys=True,
                   indent=4, separators=(',',': '))

    return 0

if __name__ == '__main__':
    exit(main(sys.argv[1:]))
