#!/usr/bin/env python

import renumFromMuscleFasta
import py.test

def test_main_basic():
    args = ['-r', '1ubq.pdb.gz,A', '-i', 'seqs.afa']
    assert renumFromMuscleFasta.main(args) == 0
