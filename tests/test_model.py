# -*- coding: utf-8 -*-

import unittest

import mirmap

from mirmap import seed, miRmap, utils


class TestModel(unittest.TestCase):
  _mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
  _mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')
  maxDiff = None

  def test_init(self):
    with self.assertRaises(TypeError):
      #: No Args
      miRmap()
      miRmap(seq_mir="TEST")
      miRmap(seq_mrn="TEST")

    temp = miRmap(seq_mrn="TEST_MRN", seq_mir="TEST_MIR")
    self.assertIsInstance(temp, miRmap, msg="Instance Ok")
    self.assertEqual(temp.seq_mrn, "TEST_MRN", msg="mRNA Ok")
    self.assertEqual(temp.seq_mir, "TEST_MIR", msg="miRNA Ok")

    default_models = {
      'full_seed6': {},
      'full_seed7': {},
      'python_only_seed6': {},
      'python_only_seed7': {}
    }
    self.assertEqual(temp.models.keys(), default_models.keys())

  def test_seed_init(self):
    temp = miRmap(seq_mrn="augcaugc", seq_mir="augc")
    self.assertIsInstance(temp._seed, seed.mmSeed)
    self.assertEqual(temp._seed.target_seq, "AUGCAUGC")
    self.assertEqual(temp._seed.mirna_seq, "AUGC")
    self.assertEqual(temp._seed.min_target_length, len("AUGC"))

    temp_seed_args = {}
    temp_seed_args['seed_args'] = {
      'target_seq': 'augcaugcaugc',   #: Wont be updated
      'mirna_seq': 'cgua',            #: Wont be updated
      'min_target_length': 2          #: Will be updated
    }
    temp = miRmap(
      seq_mrn="augcaugc", seq_mir="augc", **temp_seed_args
    )

    self.assertNotEqual(temp._seed.target_seq, "AUGCAUGCAUGC")
    self.assertNotEqual(temp._seed.mirna_seq, "CGUA")

    self.assertEqual(temp._seed.target_seq, "AUGCAUGC")
    self.assertEqual(temp._seed.mirna_seq, "AUGC")

    self.assertNotEqual(temp._seed.min_target_length, len("AUGC"))
    self.assertEqual(temp._seed.min_target_length, 2)

    temp_seed_args['seed_args'] = {
      'min_target_length': 10          #: Will be updated
    }
    temp = miRmap(
      seq_mrn="augcaugc", seq_mir="augc", **temp_seed_args
    )

    self.assertNotEqual(temp._seed.min_target_length, len("AUGC"))
    self.assertEqual(temp._seed.min_target_length, 10)

  def test_model_property(self):
    temp = miRmap(seq_mrn="augcaugc", seq_mir="augc")

    if getattr(temp, '_thermodynamic', False):
      seed_def = 'full_seed'
    else:
      seed_def = 'python_only_seed'

    self.assertEqual(temp.model, seed_def)
    self.assertEqual(temp.model_select(6),
                     temp.models.get(seed_def + '6'))
    self.assertEqual(temp.model_select(7),
                     temp.models.get(seed_def + '7'))
    self.assertEqual(temp.model_select(8),
                     temp.models.get(seed_def + '7'))

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
    temp = miRmap(seq_mrn="TEST_MRN", seq_mir="TEST_MIR")

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
    obj = miRmap(**args)
    obj.routine()

    self.assertIsInstance(temp.scores, list)

  def test_score(self):
    _mirs = utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = utils.load_fasta('tests/input/NM_024573.fa')
    obj = miRmap(
      seq_mrn=_mrnas['NM_024573'], seq_mir=_mirs['hsa-miR-30a-3p']
    )

