# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""
Interface classes with the `Vienna RNA <http://www.tbi.univie.ac.at/RNA>`
executable programs.
"""

import re
import subprocess
import tempfile

from shutil import which


class RNAvienna(object):
  """Interface class for RNA programs from Vienna."""

  def __init__(self):
    if not which("RNAfold"):
      raise EnvironmentError("RNAfold Vienna is required for Thermodynamics.")

  def fold(self, seq, **kwargs):
    return self._fold([seq], 'RNAfold', **kwargs)

  def cofold(self, seq1, seq2, **kwargs):
    return self._fold([seq1, seq2], 'RNAcofold', **kwargs)

  def _fold(self, seqs, prog, **kwargs):
    cmd = [format(prog), "--noPS"]
    regex = r'.+\n(?P<mfe_structure>\S+) \((?P<mfe>.+)\)'

    if 'constraints' in kwargs:
      cmd.append('--constraint')

    if kwargs.get('partfunc', False):
      cmd.append('--partfunc')
      if prog == 'RNAfold':
        regex += (
          r'\n *(?P<efe_structure>\S+) \[(?P<efe>.+)\]'
          r'\n(?P<cfe_structure>\S+) \{ *(?P<cfe>\S+) '
          r'd=(?P<dist>\S+)\}\n frequency of mfe struc'
          r'ture in ensemble (?P<mfe_frequency>\S+); e'
          r'nsemble diversity (?P<ensemble_diversity>\S+)'
        )
      elif prog == 'RNAcofold':
        regex += (
          r'\n *(?P<efe_structure>\S+) \[(?P<efe>.+)\]'
          r'\n frequency of mfe structure in ensemble '
          r'(?P<mfe_frequency>\S+) , delta G binding= '
          r'*(?P<efe_binding>\S+)'
        )

    if 'temperature' in kwargs:
      cmd.append('--temp=' + str(kwargs.get('temperature')))

    p = subprocess.Popen(
      cmd,
      stdin=subprocess.PIPE,
      stdout=subprocess.PIPE,
      cwd=tempfile.gettempdir()
    )

    stdout, stderr = p.communicate(
      '\n'.join(['&'.join(seqs), kwargs.get('constraints', '')]).encode()
    )

    decoded = re.match(regex, stdout.decode())
    result = {}
    for k, v in decoded.groupdict().items():
      if 'structure' in k:
        result[k] = v
      else:
        result[k] = float(v)

    return result
