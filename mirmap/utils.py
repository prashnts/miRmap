# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2013 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Common functions."""

import collections
import itertools
import functools

from Bio import SeqIO, Seq
from Bio.Alphabet import generic_rna, generic_dna


try:
  # Python 3 support
  from itertools import zip_longest
  maketrans = ''.maketrans
except (AttributeError, ImportError):
  # fallback for Python 2
  from string import maketrans
  from itertools import izip_longest
  zip_longest = izip_longest


def grouper(n, iterable, fillvalue=None):
  """http://docs.python.org/library/itertools.html#recipes"""
  args = [iter(iterable)] * n
  return zip_longest(fillvalue=fillvalue, *args)


def flatten(l1d):
  return list(itertools.chain.from_iterable(l1d))


def clean_seq(seq, alphabet):
  alphabet_real = list(set(seq))
  for l in alphabet_real:
    if l not in alphabet:
      seq = seq.replace(l, '')
  return seq


def load_fasta(fasta, as_string=False):
  if as_string:
    ff = fasta.split('\n')
    seqs = collections.OrderedDict()
    for line in ff:
      line = line.strip()
      if line.startswith('>'):
        name_seq = line[1:].strip()
        seqs[name_seq] = ''
      else:
        seqs[name_seq] += line.upper()
    return seqs
  else:
    with open(fasta) as minion:
      rec = SeqIO.to_dict(SeqIO.parse(minion, "fasta"))
      out = {k: str(v.seq) for k, v in rec.items()}
      return out


def reverse_complement(seq):
  alphabet = generic_rna if 'U' in seq or 'u' in seq else generic_dna
  s = Seq.Seq(seq, alphabet).reverse_complement()
  return str(s)


def rgetattr(obj, attr):
  return functools.reduce(getattr, [obj]+attr.split('.'))


def gen_dot_pipe_notation(pairing):
  """Returns the pairing with the dots and pipes notation"""
  string = ''
  for i in pairing[::-1]:
    if i == 0:
      string += '.'
    else:
      string += '|'
  return string


def gen_dot_bracket_notation(pairing):
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
