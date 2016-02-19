# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Evolutionary features."""

import collections
import copy

import dendropy

from mirmap.phast import Phast


def get_coord_vec(seq, alphabet, shift=None):
  if shift is None:
    shift = 1
  coord_vec = []
  i_s = 0
  for s in seq:
    if s in alphabet:
      coord_vec.append(i_s + shift)
    i_s += 1
  return coord_vec


def find_all(s, sub, indices=None, offset=None):
  if indices is None:
    indices = []
  if offset is None:
    offset = 0
  i = s.find(sub, offset)
  while i >= 0:
    indices.append(i)
    i = s.find(sub, i + 1)
  return indices


def remove_gap_column(aln):
  clean_aln = copy.copy(aln)
  keep_cols = []
  for i in range(len(aln.values()[0])):
    only_gap = True
    for seq in aln.values():
      if seq[i] != '-':
        only_gap = False
        continue
    if only_gap is False:
      keep_cols.append(i)
  for seq_name in aln.keys():
    clean_aln[seq_name] = ''.join([aln[seq_name][i] for i in keep_cols])
  return clean_aln


class mmEvolution(object):
  """
  Compute Evolutionary features

  Args:
    aln_fname (str): Alignment filename.
    aln (str): Alignment it-self.
    aln_format (str): Alignment format. Currently supported is FASTA.
    aln_alphabet (list): List of nucleotides to consider in the aligned
      sequences (others get filtered).
    subst_model (str): PhyloFit substitution model (REV...).
    tree (str): Tree in the Newick format.
    fitting_tree (bool): Fitting or not the tree on the alignment.
    use_em (bool): Fitting or not the tree with Expectation-Maximization algorithm.
    motif_def (str): 'seed' or 'seed_extended' or 'site'.
    motif_upstream_extension (int): Upstream extension length.
    motif_downstream_extension (int): Downstream extension length.
  """

  def __init__(self, seed, **kwargs):
    self.seed = seed
    defaults = {
      # PhyloFit
      'subst_model': 'REV',
      'use_em': True,
      # PhyloP
      'method': 'SPH',
      'mode': 'CONACC',
      'aln_alphabet': ['A', 'T', 'C', 'G', 'N'],
    }
    self.phast = Phast()
    self.__dict__.update(defaults)
    self._routine_done = False

  def _eval_routine(self, setup, worker, **kwargs):
    # Parameters
    if 'aln_fname' not in kwargs and 'aln' not in kwargs:
      raise IOError('An alignment is required')
    if 'tree' not in kwargs and 'fitting_tree' not in kwargs:
      raise IOError('A tree is required')

    if 'aln_format' not in kwargs and 'aln_fname' in kwargs:
      aln_format = kwargs.get('aln_fname').split('.')[-1].upper()
      if aln_format == 'FA':
        aln_format = 'FASTA'
    if aln_format is None:
      raise ValueError('Alignment format undetected')

    # Load alignment
    if 'aln_fname' in kwargs:
      seqs = utils.load_fasta(kwargs['aln_fname'])
    else:
      seqs = utils.load_fasta(kwargs['aln'], as_string=True)

    seqs_cleaned = {}
    seqs_coords = {}
    for seq_name, seq in seqs.items():
      seqs_cleaned[seq_name] = utils.clean_seq(seq, self.aln_alphabet)
      seqs_coords[seq_name] = get_coord_vec(seq, self.aln_alphabet)

    kwargs.update({'seqs': seqs})
    args = setup(**kwargs)
    args.update({
      'aln_format': aln_format
    })

    # Reset
    out = []

    # Compute
    for its in range(len(self.seed.end_sites)):
      end_site = self.seed.end_sites[its]
      # Motif
      # start_motif and end_motif are sequence coordinates => 1-based
      start_motif, end_motif = seed.get_motif_coordinates(
        end_site, self.motif_def, self.seed.pairings[its],
        self.motif_upstream_extension, self.motif_downstream_extension,
        self.seed.min_target_length
      )
      motif = self.seed.target_seq[start_motif - 1:end_motif].replace('U', 'T')

      # Species with seed(s)
      species_with_seed = []
      for seq_name, seq in seqs_cleaned.items():
        if seq.find(motif) != -1:
          species_with_seed.append(seq_name)

      args.update({
        'species_with_seed': species_with_seed,
        'start_motif': start_motif,
        'end_motif': end_motif,
        'seqs_coords': seqs_coords
      })
      out.append(worker(**args))

    return out

  def _eval_cons_bls(self, **kwargs):
    def setup(**kwargs):
      return {
        'subst_model': kwargs.get('subst_model', self.subst_model),
        'fitting_tree': kwargs.get('fitting_tree', True),
        'use_em': kwargs.get('use_em', self.use_em),
        'fitting_tree_done': False
      }

    def worker(species_with_seed, fitting_tree_done,
               fitting_tree, subst_model, aln_format, use_em):
      if len(species_with_seed) > 1:
        # Fitting tree if necessary
        if fitting_tree_done is False:
          if fitting_tree:
            if kwargs.get('aln_fname'):
              fitted_tree = self.phast.phylofit(
                subst_model=subst_model,
                aln_fname=kwargs['aln_fname'],
                aln_format=aln_format,
                tree=kwargs['tree'],
                use_em=use_em
              )['tree']
            elif kwargs('aln'):
              fitted_tree = self.phast.phylofit(
                subst_model=subst_model,
                aln=kwargs['aln'],
                aln_format=aln_format,
                tree=kwargs['tree'],
                use_em=use_em
              )['tree']
          else:
            fitted_tree = kwargs['tree']
          fitting_tree_done = True

        # Compute BLS
        dtree = dendropy.Tree.get_from_string(
          fitted_tree, schema='newick', preserve_underscores=True
        )
        dtree.retain_taxa_with_labels(species_with_seed)
        return sum(
          [edge.length for edge in dtree.postorder_edge_iter()][:-1]
        )
      else:
        return 0.0

    self.cons_blss = self._eval_routine(setup, worker, **kwargs)
    return self.cons_blss

  def _eval_selec_phylop(self, **kwargs):
    def setup(**kwargs):
      return {
        'method': kwargs.get('method', self.method),
        'mode': kwargs.get('mode', self.mode),
        'motif_upstream_extension': kwargs.get('motif_upstream_extension', 0),
        'motif_downstream_extension':
          kwargs.get('motif_downstream_extension', 0),
        'ref_species': kwargs['seq'].keys()[0],
        'mod_fname': kwargs.get('mod_fname', None),
      }

    def worker(species_with_seed, ref_species, start_motif, end_motif,
               seqs_coords, seqs, method, mode, mod_fname):
      pval = 1.0
      if len(species_with_seed) > 1:
        # Extract alignment
        start_motif_in_aln = seqs_coords[ref_species][start_motif-1]
        end_motif_in_aln = seqs_coords[ref_species][end_motif-1]
        partial_seqs = collections.OrderedDict()

        for seq_name, seq in seqs.items():
          if seq_name in species_with_seed:
            a = start_motif_in_aln - 1
            b = end_motif_in_aln
            partial_seqs[seq_name] = seq[a:b]

        partial_seqs = remove_gap_column(partial_seqs)
        aln = '\n'.join(['> %s\n%s' % (k, v) for k, v in partial_seqs.items()])

        pval = self.phast.phylop(
          method=method,
          mode=mode,
          mod_fname=mod_fname,
          aln=aln,
          aln_format='FASTA'
        )
      return pval

    self.selec_phylops = self._eval_routine(setup, worker, **kwargs)
    return self.selec_phylops

  def routine(self, **kwargs):
    try:
      self._eval_cons_bls(**kwargs)
      self._eval_selec_phylop(**kwargs)
    except IOError:
      self.cons_blss = [0 for _ in range(len(self.seed.end_sites))]
      self.selec_phylops = [0 for _ in range(len(self.seed.end_sites))]
    self._routine_done = True

  @property
  def cons_bls(self):
    try:
      return min(self.cons_blss)
    except AttributeError:
      return min(self._eval_cons_bls())

  @property
  def selec_phylop(self):
    try:
      return min(self.selec_phylops)
    except AttributeError:
      return min(self._eval_selec_phylop())
