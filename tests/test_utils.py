# -*- coding: utf-8 -*-

import unittest

from collections import namedtuple

from mirmap import utils


class TestUtils(unittest.TestCase):
  def test_grouper(self):
    self.assertEqual(
      list(utils.grouper(3, 'ABCDEFG', 'x')),
      [('A', 'B', 'C'), ('D', 'E', 'F'), ('G', 'x', 'x')]
    )

  def test_flatten(self):
    i1 = range(2)
    i2 = range(10, 12)

    self.assertEqual(
      utils.flatten([i1, i2]),
      [0, 1, 10, 11]
    )

  def test_clean_seq(self):
    dirty = "ATGCPRASHANTAUGC"
    clean = "AGCAAAUGC"
    alpha = "AUGC"

    self.assertEqual(utils.clean_seq(dirty, alpha), clean)

  def test_load_fasta(self):
    with self.assertRaises(FileNotFoundError):
      utils.load_fasta("NONEXISTENT.fa")

    fa = utils.load_fasta('tests/input/NM_024573.fa')

    fa_content = {
      'NM_024573': (
        'CUUGAUUUAGGAGCUCUCAGUUGCAUAGAAAGAUCUGGUGAGCACCUUUUCAUCCCCAGAAAAGGAGCAC'
        'GUGAAUUGAGUCGCCUGGCGGCUCUGUACGCGCUCAGGGAAGCUUAGCUUCUUGGUGCCCAUCUACGUGC'
        'ACUGGAUGAUUUUUCUUUUGAACAUUUUGCCCCACUACACUGUUUUGGGGAUAGCUGGGUUAAGCAAGUU'
        'AAAGAUAUUUACAUUUAUAUUGGAAUUUUAGCAACUUUUUUUCAGGUUAAAUAUAUAAUUUCAAGUGCUU'
        'UUAAUGAACUUAUUUUUAAUUGGCUAGGGAGCAAAAAAUAAGUGAGUUCUGCUUUUAGUUAGUUAACCUU'
        'GUUCUUUUCUUAAAUAGUACACUGCAUGGUAUUUAAUAUUCCAGGAAGCAUGGGAUUUUAUUUUGCUUGA'
        'UUUUGGGCACAUGAAAUAAUAGCUCUAGGAAAAUGCGCAUCUUAAUGACUCUUUGUAAAGAGAGGCAUUU'
        'CUUACAACUGUGAUGUUUGCUUACAUAAAAGUUACCUCAUAAGUUAAUUCUAACUUUUAUUCUUGAAUUU'
        'UAUUUCAUUUCAAUAGCUUGUUUCAUUUGCACGCCUUUGUAUUUUGAUUGACCUGUAGAAUGGAUGUUAG'
        'GAAACUCAAAAUUGAACACAGUGAAACAAAUGGUAUUUGAAGAAAUGUAAUAUCUUUUAUAUUCUAUUUA'
        'UGAUAUCCAUAAUCAAAUGAGAUUAUUUUACCACAUAAAUGUUUUAAAUAUCAGAUUUUUAGUUUGCAGU'
        'UUUAGGAAAAUGCUUUAGAUAGAAAAGGUUCUUAUGCAUUGAAUUUGGAGUACUACCAACAAUGAAUGAA'
        'UUUAUUUUUUAUAUUCUUACACAUUUUAUUGGUCAUUGUCACAGAUAGUAAAUACUAAAAAUUUCAGGUC'
        'AGUUUGUUUUGAAACUGAAAUUGGAAAUAAAUCUGGAAAUGUUUUGUUGCACUAAAAUAAUAAAAUGAAU'
        'UGUACUG'
      )
    }

    self.assertIsInstance(fa, dict)
    self.assertEqual(fa.keys(), fa_content.keys())
    self.assertEqual(fa['NM_024573'], fa_content['NM_024573'])

  def test_reverse_complement(self):
    self.assertEqual(
      utils.reverse_complement("AUGGCCAUUGUAA"),
      "UACCGGUAACAUU"[::-1]
    )
    self.assertEqual(
      utils.reverse_complement("ATGGCCATTGTAA"),
      "TACCGGTAACATT"[::-1]
    )

  def test_rgetattr(self):
    self.assertEqual(self.maxDiff, utils.rgetattr(self, 'maxDiff'))

    t2 = namedtuple('Two', ['c'])("Pass")
    t1 = namedtuple('One', ['a', 'b'])(1, t2)

    self.assertEqual(t1.a, utils.rgetattr(t1, 'a'))
    self.assertEqual(t1.b.c, utils.rgetattr(t1, 'b.c'))
