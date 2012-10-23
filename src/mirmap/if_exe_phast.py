# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Interface classes with the `PHAST <http://compgen.bscb.cornell.edu/phast>`_ executable programs."""

import os
import re
import subprocess
import tempfile

class Phast(object):
    """Interface class for RNA programs from Vienna."""
    def __init__(self, exe_path=None):
        if exe_path is None:
            self.exe_path = ''
        else:
            self.exe_path = exe_path

    def phylop(self, method, mode, mod_fname, aln_fname=None, aln=None, aln_format=None, gff_fname=None, branch=None, prune=None):
        # Cmd
        cmd = [os.path.join(self.exe_path, 'phyloP'), '--method', method, '--mode', mode, '--msa-format', aln_format]
        # Options
        if gff_fname is not None:
            cmd.append('--features')
            cmd.append(gff_fname)
        if branch is not None:
            cmd.append('--branch')
            cmd.append(branch)
        if prune is False:
            cmd.append('--no-prune')
        # Model
        cmd.append(mod_fname)
        # MSA
        if aln_fname is not None and aln is None:
            cmd.append(aln_fname)
        elif aln_fname is None and aln is not None:
            aln_file = tempfile.NamedTemporaryFile(mode='w')
            aln_file.write(aln)
            aln_file.flush()
            cmd.append(aln_file.name)
        # Executing programs
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tempfile.gettempdir())
        stdout, stderr = p.communicate()
        aln_file.close()
        # Parsing results
        decoded = re.search(r'p-value of conservation: (?P<prob>\S+)', stdout)
        return float(decoded.groupdict()['prob'])
