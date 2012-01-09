# -*- coding: utf-8 -*-

import unittest

import mirmap
import mirmap.if_clib_rna

class TestRNAvienna(unittest.TestCase):
    def test_fold(self):
        seq = 'AUCGAUGCGAUCGAGGGGCGCCCUUAAAGCUCUGAGGCGGCCCCCCA'
        lib = mirmap.if_clib_rna.RNAvienna(library_path='libs/compiled')
        lib.init_pf_fold(len(seq))
        struc_buffer = lib.get_string_buffer(len(seq))
        result = lib.pf_fold(seq, struc_buffer)
        lib.free_pf_arrays()
        self.assertAlmostEqual(-20.24213027, result)
