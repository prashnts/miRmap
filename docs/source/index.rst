:mod:`miRmap` --- Comprehensive prediction of microRNA target repression strength
=================================================================================

:mod:`miRmap` library is Python_ library organized with ...

The :mod:`miRmap` library is a Python_ library predicting the repression strength of microRNA (miRNA) targets. The model combines:

 - **thermodynamic** features: *ΔG duplex*, *ΔG binding*, *ΔG open* and *ΔG total*,
 - **evolutionary** features: *BLS* and *PhyloP*,
 - **probabilistic** features: *P.over binomial* and *P.over exact*, and
 - **sequence-based** features: *AU content*, *UTR position* and *3' pairing*.

.. toctree::
   :maxdepth: 2

Download
========

The stable releases of :mod:`miRmap` are available at http://cegg.unige.ch/mirmap.

.. note::

 Development versions are available at http://dev.vejnar.org/mirmap.

Installation
============

Requirements
------------

The :mod:`miRmap` library has the following requirements:

1. :mod:`miRmap` requires Python_ 2.7 but it can be used with Python_ 2.6 if the :mod:`collections` module is installed (A version compatible with Python_ 2.4-2.6 is available as the `ordereddict <http://pypi.python.org/pypi/ordereddict>`_ module.).

2. For the evolutionary features, the Python_ library :mod:`DendroPy` is needed for tree manipulation. You can install `DendroPy <http://pypi.python.org/pypi/DendroPy/>`_ directly from the `Python Package Index <http://pypi.python.org>`_.

3. C librairies. A compiled version of the 3 libraries (\*.so) is included in the :mod:`miRmap` distribution. If you want/have to compile them, please follow these intructions:

 - For the thermodynamic features, the `Vienna RNA <http://www.tbi.univie.ac.at/RNA>`_ library is required.

  Download the latest `Vienna RNA`_ tarball (Versions 1.8.x were successfully tested), then do:

  ::

   cd ViennaRNA-<version>
   ./configure --without-kinfold --without-forester --without-svm --without-perl
   make
   gcc -shared -Wl,-O2 -o lib/libRNAvienna.so `find lib/ -name "*.o"` -lm

 - For the evolutionary features, the `PHAST <http://compgen.bscb.cornell.edu/phast>`_ library is required (The `CLAPACK <http://www.netlib.org/clapack>`_ has to be compiled first, please follow the instructions in Phast package).

  ::

   svn co http://compgen.bscb.cornell.edu/svnrepo/phast/trunk phast
   cd phast/src

  In the file *make-include.mk*, add the ``-DUSE_PHAST_MEMORY_HANDLER`` parameter to the line starting with ``CFLAGS += -I${INC} -DPHAST_VERSION=${PHAST_VERSION}``. Then replace the path to the CLAPACK_ and compile with:

  ::

   make CLAPACKPATH=../CLAPACK-3.2.1 sharedlib

 - For the *P.over exact* feature, the `Spatt <http://www.mi.parisdescartes.fr/~nuel/spatt>`_ library is required (You will need a working copy of `CMake <http://www.cmake.org>`_ on your system).

  Download the latest Spatt_ tarball (Version 2.0 was successfully tested), then do:

  ::

   cd spatt-<version>
   mkdir build
   cd build
   cmake -DWITH_SHARED_LIB=ON ..
   make

From the directory you compiled the C libraries:

 ::

  mv spatt-<version>/libspatt2/libspatt2.so mirmap/libs/default
  mv ViennaRNA-<version>/lib/libRNAvienna.so mirmap/libs/default
  mv phast/lib/sharedlib/libphast.so mirmap/libs/default

Usage
=====

Example with the pure Python features.

.. code-block:: python

    >>> import mirmap
    >>> seq_target = 'GCUACAGUUUUUAUUUAGCAUGGGGAUUGCAGAGUGACCAGCACACUGGACUCCGAGGUGGUUCAGACAAGACAGAGGGGAGCAGUGGCCAUCAUCC\
    ... UCCCGCCAGGAGCUUCUUCGUUCCUGCGCAUAUAGACUGUACAUUAUGAAGAAUACCCAGGAAGACUUUGUGACUGUCACUUGCUGCUUUUUCUGCGCUUCAGUAACAAGU\
    ... GUUGGCAAACGAGACUUUCUCCUGGCCCCUGCCUGCUGGAGAUCAGCAUGCCUGUCCUUUCAGUCUGAUCCAUCCAUCUCUCUCUUGCCUGAGGGGAAAGAGAGAUGGGCC\
    ... AGGCAGAGAACAGAACUGGAGGCAGUCCAUCUA'
    >>> seq_mirna = 'UAGCAGCACGUAAAUAUUGGCG'
    >>> mim = mirmap.mm(seq_target, seq_mirna)
    >>> mim.find_potential_targets_with_seed(allowed_lengths=[6,7], allowed_gu_wobbles={6:0,7:0},\
    ... allowed_mismatches={6:0,7:0}, take_best=True)
    >>> mim.end_sites                                    # Coordinate(s) (3' end) of the target site on the target sequence
    [186]
    >>> mim.eval_tgs_au(with_correction=False)           # TargetScan features manually evaluated with
    >>> mim.eval_tgs_pairing3p(with_correction=False)    # a non-default parameter.
    >>> mim.eval_tgs_position(with_correction=False)
    >>> mim.prob_binomial                                # mim's attribute: the feature is automatically computed
    0.03311825751646191
    >>> print mim.report()
    155                            186
    |                              |
    CAGGAAGACUUUGUGACUGUCACUUGCUGCUUUUUCUGCGCU
                            |||||||.
              GCGGUUAUAAAUGCACGACGAU
      AU content                     0.64942
      UTR position                   166.00000
      3' pairing                     1.00000
      Probability (Binomial)         0.03312

With the C libraries installed:

.. code-block:: python

    >>> import mirmap.library_link
    >>> mim.libs = mirmap.library_link.LibraryLink('libs/compiled') # Change to the path where you unzipped the *.so files
    >>> mim.dg_duplex
    -13.5
    >>> mim.dg_open
    12.180591583251953
    >>> mim.prob_exact
    0.06798900807193115
    >>> print mim.report()
    155                            186
    |                              |
    CAGGAAGACUUUGUGACUGUCACUUGCUGCUUUUUCUGCGCU
                            |||||||.
              GCGGUUAUAAAUGCACGACGAU
      ΔG duplex (kcal/mol)          -13.50000
      ΔG binding (kcal/mol)         -11.91708
      ΔG open (kcal/mol)             12.18059
      AU content                     0.64942
      UTR position                   166.00000
      3' pairing                     1.00000
      Probability (Exact)            0.06799
      Probability (Binomial)         0.03312

Classes
=======

.. automodule:: mirmap
   :members:
   :inherited-members:
   :show-inheritance:
   :member-order: groupwise

.. autoclass:: mirmap.library_link.LibraryLink
   :members:

Copyright, License and Warranty
===============================

The miRmap library is:

.. container:: licenseblock

    .. centered:: |mirmap_copyright|

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 3 as published by
    the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (in /LICENCE).  If not, see `GNU Licenses <http://www.gnu.org/licenses/>`_.

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
