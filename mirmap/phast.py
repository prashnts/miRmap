# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""
Interface classes with the `PHAST <http://compgen.bscb.cornell.edu/phast>`
executable programs.
"""

import re
import subprocess
import tempfile

from mirmap.vienna import which


class Phast(object):

  def __init__(self):
    if not which("phast"):
      raise EnvironmentError("PHAST is required for Phylogenetic Models.")

  def phylofit(self, **kwargs):
    cmd = [
      'phast', 'phyloFit', '--precision', 'HIGH',
      '--out-root', '-', '--msa-format', kwargs.get('aln_format')
    ]

    if 'subst_model' in kwargs:
      cmd.append('--subst-mod')
      cmd.append(kwargs.get('subst_model'))

    if kwargs.get('use_em', False):
      cmd.append('--EM')

    tree_file = None
    if 'tree' in kwargs:
      tree_file = tempfile.NamedTemporaryFile(mode='w')
      tree_file.write(kwargs.get('tree'))
      tree_file.flush()
      cmd.append('--tree')
      cmd.append(tree_file.name)

    aln_file = None
    if 'aln_fname' in kwargs and 'aln' in kwargs:
      cmd.append(kwargs.get('aln_fname'))
    elif 'aln_fname' not in kwargs and 'aln' in kwargs:
      aln_file = tempfile.NamedTemporaryFile(mode='w')
      aln_file.write(kwargs.get('aln'))
      aln_file.flush()
      cmd.append(aln_file.name)

    p = subprocess.Popen(
      cmd,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      cwd=tempfile.gettempdir()
    )

    stdout, stderr = p.communicate()

    if tree_file is not None:
      tree_file.close()
    if aln_file is not None:
      aln_file.close()

    reg = (
      r'ALPHABET: (?P<alphabet>[^\n]+)\nORDER: (?P<order>\S+)'
      r'\nSUBST_MOD: (?P<subst_mod>\S+)\nTRAINING_LNL: (?P<tr'
      r'aining_lnl>\S+)\nBACKGROUND: (?P<background>[^\n]+)\n'
      r'RATE_MAT:\n(?P<rate_mat>.+)\nTREE: (?P<tree>.+;)'
    )

    decoded = re.match(reg, stdout.decode(), re.DOTALL)
    result = decoded.groupdict()
    result['mod_raw'] = stdout
    result['training_lnl'] = float(result['training_lnl'])
    return result

  def phylop(self, method, mode, mod_fname, **kwargs):

    cmd = [
      'phast', 'phyloP', '--method', method, '--mode', mode,
      '--msa-format', kwargs.get('aln_format')
    ]

    if 'gff_fname' in kwargs:
      cmd.append('--features')
      cmd.append(kwargs['gff_fname'])
    if 'branch' in kwargs:
      cmd.append('--branch')
      cmd.append(kwargs['branch'])
    if kwargs.get('prune', False):
      cmd.append('--no-prune')

    cmd.append(mod_fname)

    aln_file = None
    if 'aln_fname' in kwargs and 'aln' not in kwargs:
      cmd.append(kwargs['aln_fname'])
    elif 'aln_fname' not in kwargs and 'aln' in kwargs:
      aln_file = tempfile.NamedTemporaryFile(mode='w')
      aln_file.write(kwargs['aln'])
      aln_file.flush()
      cmd.append(aln_file.name)

    p = subprocess.Popen(
      cmd,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      cwd=tempfile.gettempdir()
    )

    stdout, stderr = p.communicate()

    if aln_file is not None:
      aln_file.close()

    reg = r'p-value of conservation: (?P<prob>\S+)'

    decoded = re.search(reg, stdout.decode())
    return float(decoded.groupdict()['prob'])
