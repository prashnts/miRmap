# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Interface classes with the `Vienna RNA <http://www.tbi.univie.ac.at/RNA>`_ executable programs."""

import os
import re
import subprocess

class RNAvienna2(object):
    """Interface class for RNA programs from Vienna"""
    def __init__(self, exe_path=None):
        if exe_path is None:
            self.exe_path = ''
        else:
            self.exe_path = exe_path

    def fold(self, seq):
        #'--noClosingGU' '--noLP' '--partfunc'
        p = subprocess.Popen([os.path.join(self.exe_path, 'RNAfold'), '--noPS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        result, stderr = p.communicate(seq)
        print result
        decoded = re.match('.+\n(\S+) \((.+)\)', result)
        return float(decoded.group(2)), decoded.group(1)


sequence = 'ACUCGGGUAGCGCGUACGCGGCGUAGCGCAUC'
rc = RNAvienna2(exe_path='/home/cegg/vejnar/tmp/ViennaRNA-2.0.7/Progs/')
folding = rc.fold(sequence)
print folding[0]
print folding[1]


import os

from ctypes import *

# Classes for structures

class DUPLEX(Structure):
    _fields_ = [('i', c_int), ('j', c_int), ('structure', c_char_p), ('energy', c_float)]

class PU_CONTRIB(Structure):
    _fields_ = [('H', POINTER(POINTER(c_double))), ('I', POINTER(POINTER(c_double))), ('M', POINTER(POINTER(c_double))), ('E', POINTER(POINTER(c_double))), ('length', c_int)]

class PLIST(Structure):
    _fields_ = [('i', c_int), ('j', c_int), ('p', c_float)]

# - part_func_co.h

class COFOLDF(Structure):
    _fields_ = [('F0AB', c_double), ('FAB', c_double), ('FcAB', c_double), ('FA', c_double), ('FB', c_double)]

class RNAvienna2_ctypes(object):
    """Interface class for RNA library from Vienna"""
    def __init__(self, library_path=None, library_name=None):
        if library_name is None:
            lib = 'libRNAvienna.so'
        else:
            lib = library_name
        if library_path is not None:
            lib = os.path.join(library_path, lib)
        self._library = cdll.LoadLibrary(lib)
        # Functions arguments and result types
        self._library.fold.argtypes = [c_char_p, c_char_p]
        self._library.fold.restype = c_float
        self._library.Lfold.argtypes = [c_char_p, c_char_p, c_int]
        self._library.Lfold.restype = c_float
        self._library.init_pf_fold.argtypes = [c_int]
        self._library.init_pf_fold.restype = None
        self._library.pf_fold.argtypes = [c_char_p, c_char_p]
        self._library.pf_fold.restype = c_float
        self._library.pfl_fold.argtypes = [c_char_p, c_int, c_int, c_float, POINTER(c_double)]
        self._library.pfl_fold.restype = POINTER(PLIST)
        self._library.free_pf_arrays.argtypes = []
        self._library.free_pf_arrays.restype = None
        self._library.cofold.argtypes = [c_char_p, c_char_p]
        self._library.cofold.restype = c_float
        self._library.duplexfold.argtypes = [c_char_p, c_char_p]
        self._library.duplexfold.restype = DUPLEX
        self._library.pf_unstru.argtypes = [c_char_p, c_int]
        self._library.pf_unstru.restype = POINTER(PU_CONTRIB)
        self._library.space.argtypes = [c_int]
        self._library.space.restype = c_void_p
        # - part_func_co.h
        self._library.co_pf_fold.argtypes = [c_char_p, c_char_p]
        self._library.co_pf_fold.restype = COFOLDF
        self._library.free_co_pf_arrays.argtypes = []
        self._library.free_co_pf_arrays.restype = None
        self._library.init_co_pf_fold.argtypes = [c_int]
        self._library.init_co_pf_fold.restype = None

    def fold(self, seq, struc_buffer):
        return self._library.fold(seq, struc_buffer)

    def lfold(self, seq, struc_buffer, maxdist):
        return self._library.Lfold(seq, struc_buffer, maxdist)

    def duplexfold(self, seq1, seq2):
        return self._library.duplexfold(seq1, seq2)

    def init_pf_fold(self, length):
        self._library.init_pf_fold(length)

    def pf_fold(self, seq, struc_buffer):
        return self._library.pf_fold(seq, struc_buffer)

    def pfl_fold(self, seq, size_window, size_max_interact, size_free, ener_cutoff):
        pup = cast(self.space((len(seq)+1)*sizeof(c_double)), POINTER(c_double))
        if size_free > 0:
            pup[0] = size_free
        return self._library.pfl_fold(seq, size_window, size_max_interact, ener_cutoff, pup)

    def pfl_fold_with_pup(self, seq, size_window, size_max_interact, size_free, ener_cutoff):
        pup = cast(self.space((len(seq)+1)*sizeof(c_double)), POINTER(c_double))
        if size_free > 0:
            pup[0] = size_free
        pl = self._library.pfl_fold(seq, size_window, size_max_interact, ener_cutoff, pup)
        return pup, pl

    def free_pf_arrays(self):
        self._library.free_pf_arrays()

    def cofold(self, seq, struc_buffer):
        return self._library.cofold(seq, struc_buffer)

    def free_co_pf_arrays(self):
        self._library.free_co_pf_arrays()

    def pf_unstru(self, seq, length):
        return self._library.pf_unstru(seq, length)

    def init_co_pf_fold(self, length):
        self._library.init_co_pf_fold(length)

    def co_pf_fold(self, seq, struc_buffer):
        return self._library.co_pf_fold(seq, struc_buffer)

    def set_temperature(self, temperature):
        c_double.in_dll(self._library, 'temperature').value = temperature

    def set_dangles(self, dangles):
        c_int.in_dll(self._library, 'dangles').value = dangles

    def set_nolonelypairs(self, nolonelypairs):
        c_int.in_dll(self._library, 'noLonelyPairs').value = nolonelypairs

    def set_pf_scale(self, pf_scale):
        c_double.in_dll(self._library, 'pf_scale').value = pf_scale

    def set_fold_constrained(self, fold_constrained):
        c_int.in_dll(self._library, 'fold_constrained').value = fold_constrained

    def set_cut_point(self, cut_point):
        c_int.in_dll(self._library, 'cut_point').value = cut_point

    def get_pf_scale(self):
        return c_double.in_dll(self._library, 'pf_scale').value

    def get_iindx(self):
        return POINTER(c_int).in_dll(self._library, 'iindx')

    def get_pr(self):
        return POINTER(c_double).in_dll(self._library, 'pr')

    def get_string_buffer(self, length):
        return create_string_buffer(length)

    def space(self, size):
        return self._library.space(size)



rc = RNAvienna2_ctypes(library_path='/home/cegg/vejnar/tmp/ViennaRNA-2.0.7/lib/')
struc_buffer = rc.get_string_buffer(len(sequence)+1)
#rc.set_fold_constrained(0)
print rc.fold(sequence, struc_buffer)
print struc_buffer.value