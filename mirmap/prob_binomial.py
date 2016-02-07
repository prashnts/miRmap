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

from mirmap import seed, prob, utils

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
    seed (seed.mmSeed)*: Seed Instance
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
    alph = list(set(self.seed.target_seq))
    self.__dict__.update({
      'markov_order': 1,
      'alphabet': alph,
      'transitions': prob.get_transitions(self.seed.target_seq, alph, 1),
      'motif_def': None,
      'motif_upstream_extension': 0,
      'motif_downstream_extension': 0,
    })
    self.__dict__.update(kwargs)

  def _eval_prob(self, worker):
    prob_ev = []
    # Compute
    for its in range(len(self.seed.end_sites)):
      end_site = self.seed.end_sites[its]

      # start_motif and end_motif are sequence coordinates => 1-based
      start_motif, end_motif = seed.get_motif_coordinates(
        end_site, self.motif_def, self.seed.pairings[its],
        self.motif_upstream_extension, self.motif_downstream_extension,
        self.seed.min_target_length
      )

      motif = self.seed.target_seq[start_motif - 1:end_motif]
      prob_ev.append(worker(motif))

    return prob_ev

  def _eval_prob_binomial(self):
    def worker(motif):
      return sum([
        1.0,
        -1 * binom_cdf(
          self.seed.target_seq.count(motif),
          self.seed.len_target_seq - len(motif) + 1,
          prob.prob_motif(
            motif, self.alphabet, self.markov_order, self.transitions
          )
        )
      ])
    self.prob_binomials = self._eval_prob(worker)
    return self.prob_binomials

  def _eval_prob_exact(self):
    def worker(motif):
      return get_exact_prob(
        motif=utils.clean_seq(motif, self.alphabet),
        nobs=self.target_seq.count(motif),
        length_seq=self.len_target_seq,
        alphabet=self.alphabet,
        transitions=self.transitions,
        markov_order=self.markov_order,
        direction='o'
      )
    self.prob_binomials = self._eval_prob(worker)

  @property
  def prob_binomial(self):
    """*P.over binomial* score with default parameters."""
    try:
      return min(self.prob_binomials)
    except AttributeError:
      return min(self._eval_prob_binomial())
