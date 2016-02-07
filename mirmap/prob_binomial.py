# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2013 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Probability based on binomial distribution feature."""

import functools
import operator as op
import math

from mirmap import iseed, prob

try:
  #: Fix for Python 2
  range = xrange
except NameError:
  pass


def nCr(n, r):
  r = min(r, n - r)
  if r <= 0:
    return 1
  numer = functools.reduce(op.mul, range(n, n-r, -1))
  denom = math.factorial(r)
  return numer//denom


def binom_pmf(k, n, p):
  return nCr(n, k) * (p**k) * ((1 - p)**(n - k))


def binom_cdf(k, n, p):
  cump = 0.0
  for ik in range(k + 1):
    cump += binom_pmf(ik, n, p)
  return cump


class mmProbBinomial(object):
  """
  Computes the *P.over binomial* score.

  Args:
    seed (iseed.mmSeed)*: Seed Instance
    markov_order (int): Markov Chain order
    alphabet (list): List of nucleotides to consider in the sequences
      (others get filtered).
    transitions (list): Transition matrix of the Markov Chain model
    motif_def (str): 'seed' or 'seed_extended' or 'site'.
    motif_upstream_extension (int): Upstream extension length.
    motif_downstream_extension (int): Downstream extension length.
  """

  def __init__(self, seed, **kwargs):
    self.seed = seed
    self.__dict__.update({'markov_order': 1})
    self.__dict__.update(kwargs)

  def _eval_prob_binomial(self):
    # Reset
    self.prob_binomials = []
    # Compute
    for its in range(len(self.seed.end_sites)):
      end_site = self.seed.end_sites[its]

      # start_motif and end_motif are sequence coordinates => 1-based
      start_motif, end_motif = iseed.get_motif_coordinates(
        end_site, self.motif_def, self.seed.pairings[its],
        self.motif_upstream_extension, self.motif_downstream_extension,
        self.seed.min_target_length
      )

      motif = self.seed.target_seq[start_motif - 1:end_motif]

      self.prob_binomials.append(sum([
        1.0,
        -1 * binom_cdf(
          self.seed.target_seq.count(motif),
          self.seed.len_target_seq - len(motif) + 1,
          prob.prob_motif(
            motif, self.alphabet, self.markov_order, self.transitions
          )
        )
      ]))
    return self.prob_binomials

  @property
  def prob_binomial(self):
    """*P.over binomial* score with default parameters."""
    try:
      return min(self.prob_binomials)
    except AttributeError:
      return min(self._eval_prob_binomial())
