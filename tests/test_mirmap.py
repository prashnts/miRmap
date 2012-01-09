# -*- coding: utf-8 -*-

import unittest

import mirmap

class TestRNAvienna(unittest.TestCase):
    _mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')

    def test_mirmap(self):
        mim = mirmap.mm(self._mrnas['NM_024573'], self._mirs['hsa-miR-30a-3p'])
        mim.find_potential_targets_with_seed(allowed_lengths=[6,7], allowed_gu_wobbles={6:0,7:0}, allowed_mismatches={6:0,7:0}, take_best=True)
        result = mim.end_sites
        self.assertEqual([931], result)
