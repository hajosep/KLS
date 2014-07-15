#!/usr/bin/env python

from __future__ import print_function
import sys
import os
# from itertools import imap
# from itertools import chain
# import subprocess
# import re
# import urllib2
# from string import rjust
# This module is from Biopython but that takes a lot of installing so I just
# plucked out the shared lib I need.  Of course, this isn't portable.
# from Bio import pairwise2
# import pairwise2
import string
import argparse
# import logging
import gzip
from collections import OrderedDict
import itertools
import pairwise2
import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')
import re
import inspect

from res_naming import *

# Given PDB structures, generate an MSA from relevant chains and then rename
# and renumber all the relevant chains accoringly.

# A fasta is required to properly generate the MSA. Generating MSA straight
# from PDB files could leave ambiguities. Suppose every sequence contains two
# alanines in a row and one sequence has missing density at one of the
# alanines.  Which alanine is missing? The first or second? Fastas don't have
# missing density.

# This approach requires that missing denisty be matched to the fasta. This is
# first attempted with REMARK 465. As a fall back Calpha distances are tested.
# The fasta can come from SEQRES annotation inside a PDB or by querying
# RCSB for a fasta file. Using SEQRES requires annotated PDBs and using
# RCSB requires pdbcodes in filenames. If your PDB is properly annotated, it
# also has the PDB code so the second method is less restrictive and is
# what this script uses.

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
                raise ValueError('Malformed FASTA file line: ' + l)
            chain = l[6]
            if chain in seqs:
                raise ValueError('Wtf? Redefinion of chain!')
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
        return longer_names[modres[threeLetter]]
        warning('Interpretting ' + threeLetter + ' as a modified ' + unmod)
    else:
        return None


def renumberChainATOMs(chain_alignment, lead_it, misdens, resNumIndex0, align_ndx, misdens_ndx):
    trail_it = lead_it.__copy__()
    line = next(lead_it)
    chain = line[21]
    chainlines = []
    align_last = len(chain_alignment)
    while (align_last !=0 and chain_alignment[align_last - 1] != '-'):
        align_last -= 1
    #chain_alignment.append('-')
    thisResNum = int(line[22:26].strip())
    lastResNum = thisResNum
    try:
        recordName = line[:6]
        #print('chain_alignment:', chain_alignment)
        #print('misdens:', misdens)
        while (line[21]==chain and (recordName == 'ATOM  ' or recordName == 'HETATM' or
               recordName == 'ANISOU' or recordName == 'TER   ')):
            # fix guard conditions
            thisResNum = int(line[22:26].strip())
            if thisResNum != lastResNum:
                lastResNum = thisResNum
                align_ndx += 1
            while align_ndx < len(chain_alignment):
                # If this is the first line a new residue number, increment the align_ndx
                # alignment index + offset due to reference sequence gives the new residue number
                if chain_alignment[align_ndx] == '-':
                    align_ndx += 1
                elif misdens_ndx < len(misdens) and thisResNum > misdens[misdens_ndx]:
                    misdens_ndx += 1
                    align_ndx += 1
                else:
                    break
            # We're on a line that has a residue we need to renumber, our alignment index is pointing
            # at the relevant alignment, and the aligned residue is not missing density
            newResNum = align_ndx + resNumIndex0
            # check that our aligned seq agrees with the residue that is being renumbered
            resName3 = line[17:20]
            if (align_ndx < len(chain_alignment)) and (threeLetterToOne(resName3) != chain_alignment[align_ndx]):
                print('line#', inspect.currentframe().f_back.f_lineno)
                print('thisResNum:', thisResNum)
                print('newResNum:', newResNum)
                print('Aligned seq doesn\'t match pdb')
                print('alignment expects:', chain_alignment[align_ndx])
                print('line has name3:', resName3)
                print('line has name1:', threeLetterToOne(resName3))
                print('line has:', threeLetterToOne(resName3))
                print('line:', line[:-1])
                print('next:', next(lead_it))
                print('align up to here:\n' + ''.join(chain_alignment[:(align_ndx+1)]))
                print('alignment:\n' + ''.join(chain_alignment))
                print('misdens:', misdens)
                print('misdens_ndx:', misdens_ndx)
                # print('locals():', locals())
                exit(1)
            chainlines.append(line[:22] + '%4s' % str(newResNum) + line[26:])
            trail_it = lead_it.__copy__()
            line = next(lead_it)
            recordName = line[:6]
    except StopIteration:
        pass
    align_ndx += 1
    # return an iterator that will produce the first line not-relevant to this
    # chain and the generated new lines
    return trail_it, chainlines, align_ndx, misdens_ndx


def get_parser():
    parser = argparse.ArgumentParser(description='Renumbering and rename'
            ' chains in a set of PDB files' 'accoring to a MSA.')
    parser.add_argument('--reference_ubq', '-r', type=str, required=True,
            help='Name of ubq model that determines residue number 1'
            ' formatted as <filename>,<chainID>.')
    parser.add_argument('--input', '-i', type=str, required=True,
            help='Name of aligned fasta file.')
    parser.add_argument('--outdir', '-o', type=str, default='.',
            help='Name of directory to output renumbered PDB files.')
    return parser


def main(argv=None):
    aparser = get_parser()
    args = aparser.parse_args(argv)

    if args.reference_ubq.find(',') == -1:
        print('Couldn\'t find a comma in the arg to reference_ubq,'
                ' it should be formatted as <filename>,<chainID>')
        print('Argument was seen as:', args.reference_ubq)
        exit(1)

    if not os.path.isdir(args.outdir):
        print('Argument given for outdir is not a directory')
        print('Argument seen as:', args.outdir)
        print('Interpretted as:', os.path.abspath(args.outdir))
        exit(1)

    alignments = {}
    with open(args.input) as afa:
        finname=None
        chain=None
        for l in afa:
            if l.startswith('>'):
                tokens = l[1:].strip().split('|')
                finname = tokens[0]
                chain = tokens[1]
                if finname not in alignments:
                    alignments[finname] = {chain: []}
                elif chain not in alignments[finname]:
                    alignments[finname][chain] = []
                continue
            alignments[finname][chain].extend(list(l.strip()))

    indexResNum1 = 0
    ref_finname, ref_chain = tuple(args.reference_ubq.split(','))
    if ref_finname not in alignments or ref_chain not in alignments[ref_finname]:
        print('Coudlnt\t find reference ubiquitin in aligned fasta file.')
        print('Searched for ">' + ref_finname + '|' + ref_chain + '"', 'in', args.input)
        exit(1)
    for i, x in enumerate(alignments[ref_finname][ref_chain]):
        if x != '-':
            indexResNum1 = i
            break
    resNumIndex0 = -indexResNum1 + 1

    for finname in alignments.keys():
        misdens = {}
        newResNum = resNumIndex0
        hasRemarks = False
        lastChain = None
        lastResnum = None
        file_alignment = alignments[finname]
        finbname = os.path.basename(finname)
        gzipped_input = (finbname[:-3] == '.gz')
        foutbname = finbname
        if gzipped_input:
            foutbname = foutbname[:-3]
        foutname = os.path.join(args.outdir, foutbname)

        if (os.path.abspath(foutname) == os.path.abspath(finname)):
            print('Input and output paths are the identical. Exiting')
            exit(1)

        print('Reading file:', finname)
        print('Renumbering chains:', ', '.join(file_alignment.keys()))
        print('Writing file:', foutname)

        if gzipped_input:
            fin = gzip.open(finname)
        else:
            fin = open(finname)
        lines = fin.readlines()
        fin.close()

        fout = open(foutname, 'w')
        newlines = []
        lead_it, = itertools.tee(lines, 1)
        try:
            trail_it = lead_it.__copy__()
            line = next(lead_it)
            while not (line.startswith('ATOM  ') or line.startswith('HETATM')):
                if line.startswith('REMARK'):
                    hasHeader = True
                if line.startswith('REMARK 465'):
                    m = re.match('REMARK 465     [A-Z]{3} ([A-Z])   (..[ 0-9])', line)
                    if m:
                        chain = m.group(1)
                        resNum = int(m.group(2).strip())
                        if chain in misdens:
                            misdens[chain].append(resNum)
                        else:
                            misdens[chain] = [resNum]
                newlines.append(line)
                trail_it = lead_it.__copy__()
                line = next(lead_it)

            # rewind lead_it so that next(lead_it) produces first line that begins with ATOM
            # lead_it = trail_it.__copy__()

            align_indices = {}
            misdens_indices = {}
            while True:
                recordName = line[:6]
                # reset numbering on new models
                if recordName == 'MODEL ':
                    for k in align_indices:
                        align_indices[k] = 0
                    for k in misdens_indices:
                        misdens_indices[k] = 0
                chain = line[21]
                # ignore Rosetta entries
                if (chain in file_alignment and line[0] != '#'):
                    if chain not in align_indices:
                        align_indices[chain] = 0
                    if chain not in misdens_indices:
                        misdens_indices[chain] = 0
                    chain_alignment = file_alignment[chain]
                    chain_misdens = misdens[chain] if chain in misdens else []
                    # keep track of the align_ndx and misdens_ndx we stopped at in case we have to continue numbering this chain
                    lead_it, chainlines, align_ndx, misdens_ndx = renumberChainATOMs(chain_alignment, trail_it, chain_misdens, resNumIndex0, align_indices[chain], misdens_indices[chain])
                    align_indices[chain] = align_ndx
                    misdens_indices[chain] = misdens_ndx
                    newlines.extend(chainlines)
                else:
                    newlines.append(line)
                trail_it = lead_it.__copy__()
                line = next(lead_it)
        except StopIteration:
            pass

        fout.writelines(newlines)
        fout.close()

    return 0

if __name__ == '__main__':
    exit(main(sys.argv[1:]))
