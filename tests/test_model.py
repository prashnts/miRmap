# -*- coding: utf-8 -*-

import unittest
import warnings

import mirmap

from mirmap import seed, miRmap, utils, thermodynamics


class BaseTestModel(unittest.TestCase):
  def assertAlmostEqualList(self, l1, l2, places=5):
    for v1, v2 in zip(*[l1, l2]):
      self.assertAlmostEqual(v1, v2, places=places)


class TestModel(BaseTestModel):
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

  def test_warning_raised(self):
    with warnings.catch_warnings(record=True) as w:
      warnings.simplefilter("default")
      obj = miRmap(seq_mrn="TEST", seq_mir="TEST")
      try:
        self.assertEqual(w[-1].category, RuntimeWarning)
        self.assertFalse(getattr(obj, '_thermodynamic', False))
      except (AssertionError, IndexError):
        self.assertIsInstance(
          getattr(obj, '_thermodynamic'),
          thermodynamics.mmThermo
        )


class TestRealModel(BaseTestModel):
  maxDiff = None

  def test_score(self):
    _mirs = utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = utils.load_fasta('tests/input/NM_024573.fa')
    obj = miRmap(seq_mrn=_mrnas['NM_024573'],
                 seq_mir=_mirs['hsa-miR-30a-3p'])
    obj.routine()

    self.assertAlmostEqualList(obj._target_scan.tgs_aus,
                               [-0.05019, -0.11319], places=4)
    self.assertAlmostEqualList(obj._target_scan.tgs_pairing3ps,
                               [0.00312, 0.05150], places=4)
    self.assertAlmostEqualList(obj._target_scan.tgs_positions,
                               [-0.01526, -0.03190], places=4)
    self.assertAlmostEqualList(obj._prob_binomial.prob_binomials,
                               [0.07013, 0.83698], places=4)

    try:
      self.assertAlmostEqualList(obj._thermodynamic.dg_duplexs,
                                 [-7.50000, -7.70000], places=4)
      self.assertAlmostEqualList(obj._thermodynamic.dg_bindings,
                                 [-8.48000, -8.79000], places=4)
      self.assertAlmostEqualList(obj._thermodynamic.dg_opens,
                                 [14.10000, 14.10000], places=4)
    except AttributeError:
      #: Validate if python-only model is used.
      self.assertEqual(obj.model, 'python_only_seed')
