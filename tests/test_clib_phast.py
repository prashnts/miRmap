# -*- coding: utf-8 -*-

import unittest

import mirmap
import mirmap.if_clib_phast

class TestPhast(unittest.TestCase):
    def test_get_exact_prob(self):
        lib = mirmap.if_clib_phast.Phast(library_path='libs/compiled')
        aln = open('tests/input/NM_024573_ts1.fa').read()
        result = lib.phylop(method='SPH', mode='CONACC', mod_fname='tests/input/NM_024573.mod', aln=aln, aln_format='FASTA', prune=True)
        self.assertAlmostEqual(0.590965154, result)
