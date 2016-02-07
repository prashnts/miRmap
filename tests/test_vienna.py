# -*- coding: utf-8 -*-

from mirmap import vienna
from tests.test_model import BaseTestModel

#: ONLY tests for Initialization. Properties are tested in Model test.


class TestVienna(BaseTestModel):
  def test_init(self):
    try:
      vienna.RNAvienna()
    except EnvironmentError:
      pass

  def test_fold(self):
    try:
      temp = vienna.RNAvienna()
      folded = temp.fold("CCGCACAGCGGGCAGUGCCC")
      self.assertEqual(folded['mfe'], -5.0)

      folded = temp.fold("CCGCACAGCGGGCAGUGCCC", partfunc=True)
      self.assertEqual(folded['dist'], 4.39)
    except EnvironmentError:
      pass

  def test_cofold(self):
    try:
      temp = vienna.RNAvienna()
      folded = temp.cofold("CCGCACAGCGGGCAGUGCCC", "CCGCACAGCGGGCAGUGCCC")
      self.assertEqual(folded['mfe'], -21.4)

      folded = temp.cofold(
        "CCGCACAGCGGGCAGUGCCC", "CCGCACAGCGGGCAGUGCCC",
        partfunc=True
      )
      self.assertEqual(folded['mfe_frequency'], 0.467001)
    except EnvironmentError:
      pass
