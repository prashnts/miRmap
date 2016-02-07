# -*- coding: utf-8 -*-

from mirmap import thermodynamics, utils, seed, vienna
from tests.test_model import BaseTestModel

#: ONLY tests for Initialization. Properties are tested in Model test.


class TestMmThermo(BaseTestModel):
  def setUp(self):
    _mirs = utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
    _mrnas = utils.load_fasta('tests/input/NM_024573.fa')
    args = {
      'target_seq': _mrnas['NM_024573'],
      'mirna_seq': _mirs['hsa-miR-30a-3p'],
    }
    self.seed = seed.mmSeed(**args)
    self.seed.find_potential_targets_with_seed()

  def test_init(self):
    try:
      temp = thermodynamics.mmThermo(seed=self.seed)
      self.assertIsInstance(temp.fold, vienna.RNAvienna)
      self.assertEqual(temp.temperature, 37.0)

      temp = thermodynamics.mmThermo(seed=self.seed, temperature=0)
      self.assertEqual(temp.temperature, 0)
    except EnvironmentError:
      pass
