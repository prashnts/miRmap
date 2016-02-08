# -*- coding: utf-8 -*-

import unittest

from mirmap import prob


class TestProb(unittest.TestCase):
  def test_permutations(self):
    items = ["A", "B", "C"]
    permute2 = [
      "AA", "AB", "AC",
      "BA", "BB", "BC",
      "CA", "CB", "CC",
    ]
    self.assertEqual(list(prob.permutations(items, 1)), items)
    self.assertEqual(list(prob.permutations(items, 2)), permute2)

  def test_get_transitions(self):
    seq = "AUGC"
    alp = ["A", "U", "G", "C"]
    out = [
      [0.0, 1.0, 0.0, 0.0],
      [0.0, 0.0, 1.0, 0.0],
      [0.0, 0.0, 0.0, 1.0],
      [0.0, 0.0, 0.0, 0.0]
    ]
    self.assertEqual(prob.get_transitions(seq, alp, 1), out)

  def test_prob_motif(self):
    seq = "AUGC"
    alp = ["A", "U", "G", "C"]
    t = prob.get_transitions(seq, alp, 1)
    self.assertEqual(prob.prob_motif("AU", alp, 2, t), 1.0)
