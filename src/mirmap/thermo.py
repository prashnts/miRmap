# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Thermodynamics features."""

import math

from . import seed

def get_pairing_string(pairing):
    """Returns the pairing with the dots and brackets notation"""
    string = ''
    for i in pairing[::-1]:
        if i == 0:
            string += '.'
        else:
            string += '('
    for i in pairing:
        if i == 0:
            string += '.'
        else:
            string += ')'
    return string

class mmThermo(seed.mmSeed):
    def eval_dg_duplex(self, librna=None, mirna_start_pairing=None, temperature=None):
        """Computes the *ΔG duplex* and *ΔG binding* scores.

           :param librna: Link to the Vienna RNA library.
           :type librna: :class:`LibraryLink`
           :param mirna_start_pairing: Starting position of the seed in the miRNA (from the 5').
           :type mirna_start_pairing: int"""
        # Parameters
        if librna is None:
            librna = self.libs.get_library_link('rna')
        if mirna_start_pairing is None:
            mirna_start_pairing = Defaults.mirna_start_pairing
        if temperature is None:
            librna.set_temperature(Defaults.temperature)
        else:
            librna.set_temperature(temperature)
        # Reset
        self.dg_duplexs = []
        self.dg_duplex_foldings = []
        self.dg_bindings = []
        # Compute
        for its in range(len(self.end_sites)):
            # Concatenated nucleotide sequences
            concat_seq = self.target_seq[self.end_sites[its] - self.min_target_length : self.end_sites[its]] + self.mirna_seq
            # Constraint sequence
            len_no_constraints = self.min_target_length - self.seed_lengths[its] - (mirna_start_pairing - 1)
            constraints_seq = '.' * len_no_constraints + get_pairing_string(self.pairings[its]) + '.' * len_no_constraints
            # Init. library
            struc_buffer = librna.get_string_buffer(self.min_target_length * 2 + 1)
            struc_buffer.value = constraints_seq
            librna.set_fold_constrained(1)
            librna.set_cut_point(self.min_target_length + 1)
            # Computing MFE
            self.dg_duplexs.append(librna.cofold(concat_seq, struc_buffer))
            self.dg_duplex_foldings.append(struc_buffer.value)
            # Computing deltaG binding
            librna.init_co_pf_fold(len(concat_seq))
            struc_buffer.value = constraints_seq
            pffold = librna.co_pf_fold(concat_seq, struc_buffer)
            librna.free_co_pf_arrays()
            self.dg_bindings.append(pffold.FcAB - pffold.FA - pffold.FB)

    def eval_dg_open(self, librna=None, upstream_rest=None, downstream_rest=None, dg_binding_area=None, temperature=None):
        """Computes the *ΔG open* score.

           :param librna: Link to the Vienna RNA library.
           :type librna: :class:`LibraryLink`
           :param upstream_rest: Upstream unfolding length.
           :type upstream_rest: int
           :param downstream_rest: Downstream unfolding length.
           :type downstream_rest: int
           :param dg_binding_area: Supplementary sequence length to fold (applied twice: upstream and downstream).
           :type dg_binding_area: int"""
        # Parameters
        if librna is None:
            librna = self.libs.get_library_link('rna')
        if upstream_rest is None:
            upstream_rest = Defaults.upstream_rest
        if downstream_rest is None:
            downstream_rest = Defaults.downstream_rest
        if dg_binding_area is None:
            dg_binding_area = Defaults.dg_binding_area
        if temperature is None:
            librna.set_temperature(Defaults.temperature)
        else:
            librna.set_temperature(temperature)
        # Reset
        self.dg_opens = []
        # Compute
        for its in range(len(self.end_sites)):
            len_dg_open_seq = dg_binding_area * 2 + upstream_rest + downstream_rest + self.min_target_length
            # Sequence (positions are 1-based ; only when cutting the target sequence, they are changed to 0-based)
            len_polya_upstream = 0
            len_polya_downstream = 0
            start_dg_open_targetseq = 0
            end_dg_open_targetseq = 0
            start_theoretic = self.end_sites[its] - self.min_target_length - upstream_rest - dg_binding_area + 1
            end_theoretic = self.end_sites[its] + downstream_rest + dg_binding_area
            if start_theoretic < 1:
                start_dg_open_targetseq = 1
                len_polya_upstream = abs(start_theoretic) + 1
            else:
                start_dg_open_targetseq = start_theoretic
                len_polya_upstream = 0
            if end_theoretic > self.len_target_seq:
                end_dg_open_targetseq = self.len_target_seq
                len_polya_downstream = end_theoretic - self.len_target_seq
            else:
                end_dg_open_targetseq = end_theoretic
                len_polya_downstream = 0
            seq_for_dg_open = len_polya_upstream * 'A' + self.target_seq[start_dg_open_targetseq - 1 : end_dg_open_targetseq] + len_polya_downstream * 'A'
            # Constraints
            constraints_seq = '.' * dg_binding_area + 'x' * (upstream_rest + self.min_target_length + downstream_rest) + '.' * dg_binding_area
            struc_buffer = librna.get_string_buffer(len_dg_open_seq + 1)
            # dg0
            librna.set_fold_constrained(0)
            librna.init_pf_fold(len_dg_open_seq)
            dg0 = librna.pf_fold(seq_for_dg_open, struc_buffer)
            librna.free_pf_arrays()
            # dg1
            librna.set_fold_constrained(1)
            struc_buffer.value = constraints_seq
            librna.init_pf_fold(len_dg_open_seq)
            dg1 = librna.pf_fold(seq_for_dg_open, struc_buffer)
            librna.free_pf_arrays()
            # dg_open
            self.dg_opens.append(dg1 - dg0)

    def eval_dg_total(self):
        """Computes the *ΔG total* score combining *ΔG duplex* and *ΔG open* scores."""
        # Check (and run) dependencies
        if hasattr(self, 'dg_duplexs') is False:
            self.eval_dg_duplex()
        if hasattr(self, 'dg_opens') is False:
            self.eval_dg_open()
        # Reset
        self.dg_totals = []
        # Compute
        for its in range(len(self.end_sites)):
            self.dg_totals.append(self.dg_duplexs[its] + self.dg_opens[its])

    @property
    def dg_duplex(self):
        """*ΔG duplex* score with default parameters."""
        if hasattr(self, 'dg_duplexs') is False:
            self.eval_dg_duplex()
        self._dg_duplex = None
        for dg in self.dg_duplexs:
            if self._dg_duplex > dg or self._dg_duplex is None:
                self._dg_duplex = dg
        return self._dg_duplex

    @property
    def dg_open(self):
        """*ΔG open* score with default parameters."""
        if hasattr(self, 'dg_opens') is False:
            self.eval_dg_open()
        self._dg_open = 0.
        for dg in self.dg_opens:
            self._dg_open += math.exp(dg * -1.0)
        self._dg_open = math.log(self._dg_open) * -1.0
        return self._dg_open

    @property
    def dg_total(self):
        """*ΔG total* score with default parameters."""
        if hasattr(self, 'dg_totals') is False:
            self.eval_dg_total()
        self._dg_total = 0.
        for dg in self.dg_totals:
            self._dg_total += math.exp(dg * -1.0)
        self._dg_total = math.log(self._dg_total) * -1.0
        return self._dg_total

class Defaults(object):
    temperature = 37.
    mirna_start_pairing = 2
    upstream_rest = 10
    downstream_rest = 15
    dg_binding_area = 70
