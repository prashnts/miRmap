# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""
Interface classe with the `Spatt <http://www.mi.parisdescartes.fr/~nuel/spatt>`
"""

import re
import subprocess
import tempfile

from mirmap.vienna import which


class Spatt(object):
  """
  Interface class for the Spatt program.
  """

  def __init__(self):
    if not which("sspatt"):
      raise EnvironmentError("SPATT is required for Exact Probabilities.")

  def get_exact_prob(self, **kwargs):

    cmd = [
      'sspatt',
      # '--format', '%a',
      '-c', str(kwargs['nobs']),
      '-p', kwargs['motif'],
      '-l', str(kwargs['length_seq']),
      '-a', ''.join(kwargs['alphabet']),
      '-m', str(kwargs['markov_order']),
      '-d'
    ]

    # if kwargs['direction'] == 'o':
    #   cmd.append('--over')
    # elif kwargs['direction'] == 'u':
    #   cmd.append('--under')

    markovf = tempfile.NamedTemporaryFile(mode='w')
    markovf.write('\n'.join(
      [' '.join(map(str, l)) for l in kwargs['transitions']]
    ))
    markovf.flush()
    cmd.append('-M')
    cmd.append(markovf.name)

    seqf = tempfile.NamedTemporaryFile(mode='w')
    seqf.write(kwargs['motif'])
    seqf.flush()
    cmd.append(seqf.name)

    p = subprocess.Popen(
      cmd,
      stdout=subprocess.PIPE,
      cwd=tempfile.gettempdir()
    )
    stdout, stderr = p.communicate()
    markovf.close()
    seqf.close()
    print(stdout)

    reg = r'P\(N>=Nobs\)=(?P<prob>\S+)'
    decoded = re.search(reg, stdout.decode())
    return float.fromhex(decoded.groupdict()['prob'])
