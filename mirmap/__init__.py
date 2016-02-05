# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2013 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

""":class:`mm` and :class:`mmPP` base classes of :mod:`miRmap` that inherit their methods from all the modules. Each module define the methods for one category."""

from mirmap import model
from mirmap import prob_binomial
from mirmap import report
from mirmap import targetscan

class mm(model.mmModel, prob_binomial.mmProbBinomial, report.mmReport, targetscan.mmTargetScan):
    """miRNA and mRNA containing class with pure Python methods only."""
    pass
