# -*- coding: utf-8 -*-

import unittest

from mirmap import prob_binomial, seed, utils


class TestProbBinomial(unittest.TestCase):
  def test_ncr(self):
    num = (4 * 3) // (1 * 2)
    self.assertEqual(prob_binomial.nCr(4, 2), num)
    self.assertEqual(prob_binomial.nCr(4, 0), 1)

  def test_binomial_cdf(self):
    #: Cumulative Probability Distribution Function.
    out = sum([
      prob_binomial.binom_pmf(2, 4, 0.8),
      prob_binomial.binom_pmf(1, 4, 0.8),
      prob_binomial.binom_pmf(0, 4, 0.8),
    ])
    self.assertAlmostEqual(prob_binomial.binom_cdf(2, 4, 0.8), out)

  def test_binomial_pdf(self):
    self.assertAlmostEqual(prob_binomial.binom_pmf(2, 4, 0.8), 0.15359, 4)

class TestMmProbBinomial(unittest.TestCase):
  def setUp(self):
    _mirs = utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = utils.load_fasta('tests/input/NM_024573.fa')
    args = {
      'target_seq': _mrnas['NM_024573'],
      'mirna_seq': _mirs['hsa-miR-30a-3p'],
    }
    self.seed = seed.mmSeed(**args)
    self.seed.find_potential_targets_with_seed()

  def test_init(self):
    ob = prob_binomial.mmProbBinomial(self.seed)
    self.assertEqual(ob.markov_order, 1)
    self.assertEqual(ob.alphabet, list(set(self.seed.target_seq)))

    ob = prob_binomial.mmProbBinomial(
      self.seed, markov_order=10, motif_def='TEST'
    )
    self.assertEqual(ob.markov_order, 10)
    self.assertEqual(ob.motif_def, 'TEST')

  def test_eval_prob_binomial(self):
    ob = prob_binomial.mmProbBinomial(self.seed)

    with self.assertRaises(AttributeError):
      ob.prob_binomials

    t = ob._eval_prob_binomial()
    r = [0.07012894456680518, 0.8369842777406993]

    self.assertIsInstance(ob.prob_binomials, list)
    self.assertEqual(r, t)
    self.assertEqual(ob.prob_binomial, min(r))

  def test_properties(self):
    ob = prob_binomial.mmProbBinomial(self.seed)
    r = [0.07012894456680518, 0.8369842777406993]
    self.assertEqual(ob.prob_binomial, min(r))
