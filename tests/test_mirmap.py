# -*- coding: utf-8 -*-

import unittest

import mirmap
from mirmap import integrated_model
from mirmap import iseed


class TestRNAvienna(unittest.TestCase):
    _mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')

    def test_mirmap(self):
        mim = mirmap.mm(self._mrnas['NM_024573'], self._mirs['hsa-miR-30a-3p'])
        mim.find_potential_targets_with_seed(allowed_lengths=[6,7], allowed_gu_wobbles={6:0,7:0}, allowed_mismatches={6:0,7:0}, take_best=True)
        result = mim.end_sites
        self.assertEqual([931], result)


class TestIntegratedModel(unittest.TestCase):
    def test_init(self):
        with self.assertRaises(TypeError):
            #: No Args
            integrated_model.miRmap()
            integrated_model.miRmap(seq_mir="TEST")
            integrated_model.miRmap(seq_mrn="TEST")

        temp = integrated_model.miRmap(seq_mrn="TEST_MRN", seq_mir="TEST_MIR")
        self.assertIsInstance(temp, integrated_model.miRmap, msg="Instance Ok")
        self.assertEqual(temp.seq_mrn, "TEST_MRN", msg="mRNA Ok")
        self.assertEqual(temp.seq_mir, "TEST_MIR", msg="miRNA Ok")

        default_models = {
            'full_seed6': {},
            'full_seed7': {},
            'python_only_seed6': {},
            'python_only_seed7': {}
        }
        self.assertEqual(temp.models.keys(), default_models.keys())

    def test_seed_property(self):
        temp = integrated_model.miRmap(seq_mrn="augcaugc", seq_mir="augc")
        self.assertIsInstance(temp.seed, iseed.mmSeed)
        self.assertEqual(temp.seed.target_seq, "AUGCAUGC")
        self.assertEqual(temp.seed.mirna_seq, "AUGC")
        self.assertEqual(temp.seed.min_target_length, len("AUGC"))
        temp.seed = {
            'target_seq': 'augcaugcaugc',   #: Wont be updated
            'mirna_seq': 'cgua',            #: Wont be updated
            'min_target_length': 2          #: Will be updated
        }

        self.assertNotEqual(temp.seed.target_seq, "AUGCAUGCAUGC")
        self.assertNotEqual(temp.seed.mirna_seq, "CGUA")

        self.assertEqual(temp.seed.target_seq, "AUGCAUGC")
        self.assertEqual(temp.seed.mirna_seq, "AUGC")

        self.assertNotEqual(temp.seed.min_target_length, len("AUGC"))
        self.assertEqual(temp.seed.min_target_length, 2)

        temp.seed = {
            'min_target_length': 10         #: Will be updated
        }
        self.assertNotEqual(temp.seed.min_target_length, len("AUGC"))
        self.assertEqual(temp.seed.min_target_length, 10)

        temp.seed = {}  #: Reset
        self.assertEqual(temp.seed.min_target_length, len("AUGC"))

class TestISeed(unittest.TestCase):
    _mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')
    maxDiff = None

    def test_init(self):
        with self.assertRaises(KeyError):
            iseed.mmSeed()
            iseed.mmSeed(target_seq="AUGC")
            iseed.mmSeed(mirna_seq="AUGC")

        temp = iseed.mmSeed(target_seq="augcaugc", mirna_seq="augc")
        self.assertIsInstance(temp, iseed.mmSeed)
        self.assertEqual(temp.target_seq, "AUGCAUGC")
        self.assertEqual(temp.mirna_seq, "AUGC")
        self.assertEqual(temp.min_target_length, len("AUGC"))
        self.assertEqual(temp.len_mirna_seq, len("AUGC"))
        self.assertEqual(temp.len_target_seq, len("AUGCAUGC"))

        temp2 = iseed.mmSeed(target_seq="augc", mirna_seq="augc",
                             min_target_length=2)
        self.assertEqual(temp2.min_target_length, 2)

    def test_find_potential_targets_with_seed(self):
        args = {
            'target_seq': self._mrnas['NM_024573'],
            'mirna_seq': self._mirs['hsa-miR-30a-3p'],
        }
        obj = iseed.mmSeed(**args)
        self.assertIsInstance(obj, iseed.mmSeed)

        r1 = obj.find_potential_targets_with_seed()
        ar1 = {
            'end_sites': [931, 931],
            'seed_lengths': [7, 6],
            'pairings': [[0, 2, 3, 4, 5, 6, 7, 8], [0, 2, 3, 4, 5, 6, 7]],
            'nb_mismatches_except_gu_wobbles': [0, 0],
            'nb_gu_wobbles': [0, 0]
        }
        self.assertEqual(r1, ar1)

        tar_1 = {
            'mirna_start_pairing': 2,
            'allowed_lengths': [6, 7, 8],
            'allowed_gu_wobbles': {6: 0, 7: 1, 8: 2},
            'allowed_mismatches': {6: 0, 7: 0, 8: 0},
            'take_best': False,
        }
        tr1 = obj.find_potential_targets_with_seed(**tar_1)
        tr1_out = {
            'end_sites': [931, 931],
            'seed_lengths': [7, 6],
            'pairings': [[0, 2, 3, 4, 5, 6, 7, 8], [0, 2, 3, 4, 5, 6, 7]],
            'nb_mismatches_except_gu_wobbles': [0, 0],
            'nb_gu_wobbles': [0, 0]
        }
        self.assertEqual(tr1, tr1_out)
        self.assertEqual(obj.end_sites, tr1['end_sites'])
        self.assertEqual(obj.seed_lengths, tr1['seed_lengths'])
        self.assertEqual(obj.pairings, tr1['pairings'])
        self.assertEqual(obj.nb_mismatches_except_gu_wobbles,
                         tr1['nb_mismatches_except_gu_wobbles'])
        self.assertEqual(obj.nb_gu_wobbles, tr1['nb_gu_wobbles'])

        tar_2 = {
            'mirna_start_pairing': 0,   #: Changed the pairing start
            'allowed_lengths': [6, 7, 8],
            'allowed_gu_wobbles': {6: 0, 7: 1, 8: 2},
            'allowed_mismatches': {6: 0, 7: 0, 8: 0},
            'take_best': False,
        }
        tr2 = obj.find_potential_targets_with_seed(**tar_2)
        tr2_out = {
            'end_sites': [987, 987, 987],
            'seed_lengths': [8, 7, 6],
            'pairings': [[], [], []],
            'nb_mismatches_except_gu_wobbles': [0, 0, 0],
            'nb_gu_wobbles': [0, 0, 0]
        }
        self.assertEqual(tr2, tr2_out)
        self.assertEqual(obj.end_sites, tr2['end_sites'])
        self.assertEqual(obj.seed_lengths, tr2['seed_lengths'])
        self.assertEqual(obj.pairings, tr2['pairings'])
        self.assertEqual(obj.nb_mismatches_except_gu_wobbles,
                         tr2['nb_mismatches_except_gu_wobbles'])
        self.assertEqual(obj.nb_gu_wobbles, tr2['nb_gu_wobbles'])

        tar_3 = {
            'mirna_start_pairing': 2,
            'allowed_lengths': [6, 7, 8],
            'allowed_gu_wobbles': {6: 0, 7: 1, 8: 2},
            'allowed_mismatches': {6: 0, 7: 0, 8: 0},
            'take_best': True,  #: Here.
        }
        tr3 = obj.find_potential_targets_with_seed(**tar_3)
        tr3_out = {
            'end_sites': [931],
            'seed_lengths': [7],
            'pairings': [[0, 2, 3, 4, 5, 6, 7, 8]],
            'nb_mismatches_except_gu_wobbles': [0],
            'nb_gu_wobbles': [0]
        }
        self.assertEqual(tr3, tr3_out)
        self.assertEqual(obj.end_sites, tr3['end_sites'])
        self.assertEqual(obj.seed_lengths, tr3['seed_lengths'])
        self.assertEqual(obj.pairings, tr3['pairings'])
        self.assertEqual(obj.nb_mismatches_except_gu_wobbles,
                         tr3['nb_mismatches_except_gu_wobbles'])
        self.assertEqual(obj.nb_gu_wobbles, tr3['nb_gu_wobbles'])

        tar_4 = {
            'mirna_start_pairing': 2,
            'allowed_lengths': [6, 7, 8],
            'allowed_gu_wobbles': {6: 0, 7: 1, 8: 2},
            'allowed_mismatches': {6: 0, 7: 0, 8: 1},
            'take_best': False,
        }
        tr4 = obj.find_potential_targets_with_seed(**tar_4)
        tr4_out = {
            'end_sites': [937, 931, 931, 931, 796, 236],
            'seed_lengths': [8, 8, 7, 6, 8, 8],
            'pairings': [
                [0, 2, 3, 4, 5, 6, 7, 8, 0],
                [0, 2, 3, 4, 5, 6, 7, 8, 0],
                [0, 2, 3, 4, 5, 6, 7, 8],
                [0, 2, 3, 4, 5, 6, 7],
                [0, 2, 3, 4, 5, 0, 7, 8, 9],
                [0, 2, 3, 4, 5, 6, 7, 8, 0]
            ],
            'nb_mismatches_except_gu_wobbles': [1, 1, 0, 0, 1, 1],
            'nb_gu_wobbles': [2, 0, 0, 0, 1, 2]
        }
        self.assertEqual(tr4, tr4_out)
        self.assertEqual(obj.end_sites, tr4['end_sites'])
        self.assertEqual(obj.seed_lengths, tr4['seed_lengths'])
        self.assertEqual(obj.pairings, tr4['pairings'])
        self.assertEqual(obj.nb_mismatches_except_gu_wobbles,
                         tr4['nb_mismatches_except_gu_wobbles'])
        self.assertEqual(obj.nb_gu_wobbles, tr4['nb_gu_wobbles'])

class TestMiRmapModel(unittest.TestCase):
    _mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')

    def test_mmmodel_eval(self):
        mim = mirmap.mm(self._mrnas['NM_024573'], self._mirs['hsa-miR-30a-3p'])


