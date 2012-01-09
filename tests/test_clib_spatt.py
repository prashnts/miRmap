# -*- coding: utf-8 -*-

import unittest

import mirmap
import mirmap.if_clib_spatt

class TestSpatt(unittest.TestCase):
    def test_get_exact_prob(self):
        lib = mirmap.if_clib_spatt.Spatt(library_path='libs/default')
        result = lib.get_exact_prob('AUUAAAA', 1, 1000, ['A', 'U'], [0.1780821917808219, 0.821917808219178, 0.6593406593406593, 0.34065934065934067], 1, 'o')
        self.assertAlmostEqual(0.3704134091486, result)
