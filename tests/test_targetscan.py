# -*- coding: utf-8 -*-

import unittest

from collections import namedtuple

import mirmap

from mirmap import itargetscan, iseed


class TestTargetScan(unittest.TestCase):
  def setUp(self):
    _mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')
    args = {
      'target_seq': _mrnas['NM_024573'],
      'mirna_seq': _mirs['hsa-miR-30a-3p'],
    }
    self.seed = iseed.mmSeed(**args)
    self.seed.find_potential_targets_with_seed()

  def test_init(self):
    with self.assertRaises(TypeError):
      itargetscan.mmTargetScan()

    obj = itargetscan.mmTargetScan(seed=self.seed)
    self.assertIsInstance(obj, itargetscan.mmTargetScan)
    self.assertEqual(obj.seed, self.seed)

  def test_init_defaults(self):
    obj = itargetscan.mmTargetScan(seed=self.seed)
    self.assertEqual(obj.with_correction, True)
    self.assertEqual(obj.ca_window_length, 30)
    self.assertIsInstance(obj.ts_types, dict)

    type_key = ['6mer', '7mer-A1', '7mer-m8', '8mer']

    types_keys = [
      'name',
      'up_shift',
      'down_shift',
      'fc_mean',
      'ca_fc_slope',
      'ca_fc_intercept',
      'ca_weights_up',
      'ca_weights_down',
      'po_fc_slope',
      'po_fc_intercept',
      'pa_fc_slope',
      'pa_fc_intercept',
      'pa_mirna_seed_start',
      'pa_mirna_seed_overhang',
    ]
    type_val_tuple = namedtuple('TSTypes', types_keys)

    for key in type_key:
      self.assertIn(key, obj.ts_types.keys())

    for key in obj.ts_types:
      self.assertEqual(obj.ts_types[key]._fields,
                       type_val_tuple._fields)

  def test_init_args(self):
    obj = itargetscan.mmTargetScan(
      seed=self.seed,
      with_correction=False)
    self.assertEqual(obj.with_correction, False)
    self.assertEqual(obj.ca_window_length, 30)

    obj = itargetscan.mmTargetScan(
      seed=self.seed,
      with_correction=False,
      ca_window_length=10)
    self.assertEqual(obj.with_correction, False)
    self.assertEqual(obj.ca_window_length, 10)

    obj = itargetscan.mmTargetScan(
      seed=self.seed,
      ca_window_length=10)
    self.assertEqual(obj.with_correction, True)
    self.assertEqual(obj.ca_window_length, 10)

    obj = itargetscan.mmTargetScan(
      seed=self.seed,
      unwanted=True,
      ca_window_length=20)
    self.assertEqual(obj.with_correction, True)
    self.assertEqual(obj.ca_window_length, 20)

    with self.assertRaises(AttributeError):
      obj.unwanted

  def test_targetscan_ts_type(self):
    obj = itargetscan.mmTargetScan(seed=self.seed)
    self.assertEqual(obj._targetscan_ts_type(6, 'A'), '8mer')
    self.assertEqual(obj._targetscan_ts_type(6, 'T'), '7mer-m8')
    self.assertEqual(obj._targetscan_ts_type(7, 'A'), '7mer-A1')
    self.assertEqual(obj._targetscan_ts_type(8, 'A'), '7mer-A1')
    self.assertEqual(obj._targetscan_ts_type(7, 'T'), '6mer')

    with self.assertRaises(ValueError):
      obj._targetscan_ts_type(0, 'A')

  def test_eval_tgs_au(self):
    obj = itargetscan.mmTargetScan(seed=self.seed)

    with self.assertRaises(AttributeError):
      obj.tgs_aus

    t1 = obj._eval_tgs_au()
    r1 = [-0.0501892775521989, -0.11318565255796981]

    self.assertEqual(t1, r1)

  def test_eval_tgs_position(self):
    obj = itargetscan.mmTargetScan(seed=self.seed)

    with self.assertRaises(AttributeError):
      obj.tgs_positions

    t1 = obj._eval_tgs_position()
    r1 = [0.055499999999999994, 0.0995]

    self.assertEqual(t1, r1)

  def test_eval_tgs_pairing3p(self):
    obj = itargetscan.mmTargetScan(seed=self.seed)

    with self.assertRaises(AttributeError):
      obj.tgs_pairing3ps

    t1 = obj._eval_tgs_pairing3p()
    r1 = [0.003119999999999999, 0.051500000000000004]

    self.assertEqual(t1, r1)

  def test_eval_tgs_score(self):
    obj = itargetscan.mmTargetScan(seed=self.seed)

    with self.assertRaises(AttributeError):
      obj._eval_tgs_score()

    obj._eval_tgs_au()

    with self.assertRaises(AttributeError):
      obj._eval_tgs_score()

    obj._eval_tgs_position()

    with self.assertRaises(AttributeError):
      obj._eval_tgs_score()

    obj._eval_tgs_pairing3p()

    t1 = obj._eval_tgs_score()
    r1 = [-0.006569277552198908, -0.1231856525579698]

    self.assertEqual(t1, r1)
