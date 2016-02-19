# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

from mirmap.vienna import RNAvienna
from mirmap.utils import gen_dot_bracket_notation


class mmThermo(object):
  """
  Compute the Thermodynamic Properties.

  KwArgs:
    seed (model.seed): Seed
    upstream_rest (int): Upstream unfolding length.
    downstream_rest (int): Downstream unfolding length.
    dg_binding_area (int): Supplementary sequence length to fold
      (applied twice: upstream and downstream).
    temperature (float): Folding temperature.
  """

  def __init__(self, seed, **kwargs):
    self.seed = seed
    defaults = {
      'temperature': 37.0,
      'mirna_start_pairing': 2,
      'upstream_rest': 10,
      'downstream_rest': 15,
      'dg_binding_area': 70,
    }
    self.fold = RNAvienna()
    self.__dict__.update(defaults)
    self.__dict__.update(kwargs)
    self._routine_done = False

  def _eval_dg_duplex(self):
    self.dg_duplex_seeds = []
    self.dg_binding_seeds = []
    self.dg_duplexs = []
    self.dg_duplex_foldings = []
    self.dg_bindings = []
    # Compute
    for its in range(len(self.seed.end_sites)):
      # Target site and seed binding sequences
      a1 = self.seed.end_sites[its] - self.seed.min_target_length
      b1 = self.seed.end_sites[its]
      target_site_seq = self.seed.target_seq[a1:b1]

      a2 = (
        self.seed.end_sites[its] -
        (self.mirna_start_pairing - 1) -
        self.seed.seed_lengths[its]
      )
      b2 = (self.seed.end_sites[its] - (self.mirna_start_pairing - 1))
      target_seed_seq = self.seed.target_seq[a2:b2]

      a3 = self.mirna_start_pairing - 1
      b3 = self.mirna_start_pairing + self.seed.seed_lengths[its] - 1
      mirna_seed_seq = self.seed.mirna_seq[a3:b3]

      # Constraint sequence
      len_no_constraints = (
        self.seed.min_target_length -
        self.seed.seed_lengths[its] -
        (self.mirna_start_pairing - 1)
      )
      constraints_seq = (
        '.' * len_no_constraints +
        gen_dot_bracket_notation(self.seed.pairings[its]) +
        '.' * len_no_constraints
      )
      # Co-folding of seed
      result = self.fold.cofold(
        target_seed_seq,
        mirna_seed_seq,
        partfunc=True,
        temperature=self.temperature
      )

      self.dg_duplex_seeds.append(result['mfe'])
      self.dg_binding_seeds.append(result['efe_binding'])
      # Co-folding of target site

      result = self.fold.cofold(
        target_site_seq,
        self.seed.mirna_seq,
        constraints=constraints_seq,
        partfunc=True,
        temperature=self.temperature
      )
      self.dg_duplexs.append(result['mfe'])
      self.dg_duplex_foldings.append(result['mfe_structure'])
      self.dg_bindings.append(result['efe_binding'])

    return {
      'dg_duplex_seeds': self.dg_duplex_seeds,
      'dg_binding_seeds': self.dg_binding_seeds,
      'dg_duplexs': self.dg_duplexs,
      'dg_duplex_foldings': self.dg_duplex_foldings,
      'dg_bindings': self.dg_bindings,
    }

  def _eval_dg_open(self):
    """
    Computes the *ΔG open* score.
    """
    self.dg_opens = []
    # Compute
    for its in range(len(self.seed.end_sites)):
      len_polya_upstream = 0
      len_polya_downstream = 0
      start_dg_open_targetseq = 0
      end_dg_open_targetseq = 0
      start_theoretic = (
        self.seed.end_sites[its] -
        self.seed.min_target_length -
        self.upstream_rest -
        self.dg_binding_area + 1
      )

      end_theoretic = (
        self.seed.end_sites[its] +
        self.downstream_rest +
        self.dg_binding_area
      )

      if start_theoretic < 1:
        start_dg_open_targetseq = 1
        len_polya_upstream = abs(start_theoretic) + 1
      else:
        start_dg_open_targetseq = start_theoretic
        len_polya_upstream = 0

      if end_theoretic > self.seed.len_target_seq:
        end_dg_open_targetseq = self.seed.len_target_seq
        len_polya_downstream = end_theoretic - self.seed.len_target_seq
      else:
        end_dg_open_targetseq = end_theoretic
        len_polya_downstream = 0

      a4 = start_dg_open_targetseq - 1
      b4 = end_dg_open_targetseq
      seq_for_dg_open = (
        len_polya_upstream * 'A' +
        self.seed.target_seq[a4:b4] +
        len_polya_downstream * 'A'
      )

      # Constraint sequences
      c1 = (self.upstream_rest + self.seed.min_target_length +
            self.downstream_rest)
      constraints_seq = (
        '.' * self.dg_binding_area +
        'x' * c1 +
        '.' * self.dg_binding_area
      )
      # Folding
      # dg0
      result_dg0 = self.fold.fold(
        seq_for_dg_open,
        partfunc=True,
        temperature=self.temperature
      )
      # dg1
      result_dg1 = self.fold.fold(
        seq_for_dg_open,
        constraints=constraints_seq,
        partfunc=True,
        temperature=self.temperature,
      )
      # dg_open
      self.dg_opens.append(result_dg1['efe'] - result_dg0['efe'])
    return self.dg_opens

  def _eval_dg_total(self):
    """
    Computes the *ΔG total* score combining *ΔG duplex* and *ΔG open* scores.
    """
    self.dg_totals = []
    # Compute
    for its in range(len(self.seed.end_sites)):
      self.dg_totals.append(self.dg_duplexs[its] + self.dg_opens[its])

  def routine(self):
    self._eval_dg_duplex()
    self._eval_dg_open()
    self._eval_dg_total()
    self._routine_done = True

  @property
  def dg_duplex(self):
    return min(self.dg_duplexs)

  @property
  def dg_binding(self):
    return min(self.dg_bindings)

  @property
  def dg_duplex_seed(self):
    return min(self.dg_duplex_seeds)

  @property
  def dg_binding_seed(self):
    return min(self.dg_binding_seeds)

  @property
  def dg_open(self):
    return min(self.dg_opens)

  @property
  def dg_total(self):
    return min(self.dg_totals)
