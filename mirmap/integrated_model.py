# -*- coding: utf-8 -*-


class miRmap(object):
  """
  miRmap Evaluation and Prediction model.

  Args:
    seq_mir (str): miRNA sequence
    seq_mrn (str): mRNA sequence
  """

  def __init__(self, **kwargs):
    """
    Initialize the class with miRNA and mRNA Sequences.
    """
    if not ('seq_mir' in kwargs and 'seq_mrn' in kwargs):
      raise TypeError("miRNA and mRNA sequences are Required Parameters.")

    self.__dict__.update(kwargs)
