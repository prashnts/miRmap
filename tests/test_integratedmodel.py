# -*- coding: utf-8 -*-

import unittest

import mirmap

from mirmap import integrated_model, iseed


class TestIntegratedModel(unittest.TestCase):
  _mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
  _mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')
  maxDiff = None

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

  def test_model_property(self):
    temp = integrated_model.miRmap(seq_mrn="augcaugc", seq_mir="augc")

    self.assertEqual(temp.model, 'python_only_seed')
    self.assertEqual(temp.model_select(6),
                     temp.models.get('python_only_seed6'))
    self.assertEqual(temp.model_select(7),
                     temp.models.get('python_only_seed7'))
    self.assertEqual(temp.model_select(8),
                     temp.models.get('python_only_seed7'))

    temp.model = 'full_seed'
    self.assertEqual(temp.model, 'full_seed')
    self.assertEqual(temp.model_select(6),
                     temp.models.get('full_seed6'))
    self.assertEqual(temp.model_select(7),
                     temp.models.get('full_seed7'))
    self.assertEqual(temp.model_select(8),
                     temp.models.get('full_seed7'))

    with self.assertRaises(ValueError):
      temp.model_select(2)

  def test_eval_score(self):
    temp = integrated_model.miRmap(seq_mrn="TEST_MRN", seq_mir="TEST_MIR")

    with self.assertRaises(AttributeError):
      temp.scores

    with self.assertRaises(AttributeError):
      temp._eval_score()

    temp.routine()
    self.assertIsInstance(temp.scores, list)

    args = {
      'seq_mrn': self._mrnas['NM_024573'],
      'seq_mir': self._mirs['hsa-miR-30a-3p'],
    }
    obj = integrated_model.miRmap(**args)
    obj.routine()

    self.assertIsInstance(temp.scores, list)
