# -*- coding: utf-8 -*-

from mirmap import iseed, itargetscan
from mirmap.utils import rgetattr


class miRmap(object):
  """
  miRmap Evaluation and Prediction model.

  Args:
    seq_mir (str): miRNA sequence
    seq_mrn (str): mRNA sequence
    seed_args (dict): Seed Class init args
  """

  def __init__(self, **kwargs):
    """
    Initialize the class with miRNA and mRNA Sequences.
    """
    if not ('seq_mir' in kwargs and 'seq_mrn' in kwargs):
      raise TypeError("miRNA and mRNA sequences are Required Parameters.")

    self.__dict__.update(kwargs)
    self.__init_models()
    self.__init_seed(**kwargs.get('seed_args', {}))
    self.__init_targetscan()

  def __init_models(self):
    """
    Current models based on Grimson et al dataset
    """
    self.models = {
      'full_seed6': {
          'tgs_au': -0.275016235769136,
          'tgs_position': 5.43367028065211e-06,
          'tgs_pairing3p': -0.00233278119760994,
          'dg_duplex': 0.00772658898496047,
          'dg_binding': -0.00303683833660696,
          'dg_duplex_seed': 0.0496909801533612,
          'dg_binding_seed': -0.048931930580652,
          'dg_open': 0.000674676164622922,
          'prob_exact': 0.16111635592018,
          'prob_binomial': -0.0388333740708671,
          'cons_bls': -0.00426314077593848,
          'selec_phylop': -0.0112455248228072,
          'intercept': 0.148300586692704,
      },
      'full_seed7': {
          'tgs_au': -0.402470212080983,
          'tgs_position': 6.89249707831041e-05,
          'tgs_pairing3p': -0.0129891251446967,
          'dg_duplex': 0.0141332997802509,
          'dg_binding': -0.0132159175462755,
          'dg_duplex_seed': -0.0814445085121904,
          'dg_binding_seed': 0.115558118311931,
          'dg_open': 0.00331507347139685,
          'prob_exact': 0.792962156550929,
          'prob_binomial': -0.22119499646323,
          'cons_bls': -0.0355840335642203,
          'selec_phylop': -0.0127531995991629,
          'intercept': 0.349448109979275,
      },
      'python_only_seed6': {
          '_target_scan.tgs_au': -0.275594504153219,
          '_target_scan.tgs_position': 9.44582844229299e-06,
          '_target_scan.tgs_pairing3p': -0.0111209267382849,
          # 'prob_binomial': 0.0701619992923641,
          # 'cons_bls': -0.00646548621345819,
          'intercept': 0.121104869645859,
      },
      'python_only_seed7': {
          '_target_scan.tgs_au': -0.443606032336791,
          '_target_scan.tgs_position': 6.34603935320321e-05,
          '_target_scan.tgs_pairing3p': -0.0207672870210752,
          # 'prob_binomial': 0.378665477250754,
          # 'cons_bls': -0.0552713344740971,
          'intercept': 0.150015113841088,
      }
    }
    self.model = 'python_only_seed'

  def __init_seed(self, **args):
    arg_init = {
      'target_seq': self.seq_mrn,
      'mirna_seq': self.seq_mir
    }
    args.update(arg_init)
    self.__seed_m = iseed.mmSeed(**args)

  def __init_targetscan(self, **args):
    self._target_scan = itargetscan.mmTargetScan(seed=self.__seed_m)

  @property
  def seed(self):
    return self.__seed_m

  @seed.setter
  def seed(self, kwargs):
    #: Cannot update the sequences. HAVE to update the top level sequences.
    arg_init = {
      'target_seq': self.seq_mrn,
      'mirna_seq': self.seq_mir
    }
    kwargs.update(arg_init)
    self.__seed_m = iseed.mmSeed(**kwargs)

  @property
  def model(self):
    return self.__selected_model

  @model.setter
  def model(self, val):
    self.__selected_model = val

  def model_select(self, count):
    if count == 6:
      return self.models.get(self.__selected_model + "6")
    elif count >= 7:
      return self.models.get(self.__selected_model + "7")
    raise ValueError("Count gotta be greater or equal to 6.")

  def routine(self, **kwargs):
    seed_targ = kwargs.get('arg_seed_tar', {})
    self.seed.find_potential_targets_with_seed(**seed_targ)
    self._target_scan.routine()
    self._eval_score()

  def _eval_score(self):
    """
    Computes the *miRmap* score(s)
    """
    # Reset
    self.scores = []
    # Compute
    for i in range(len(self.seed.end_sites)):
      model = self.model_select(self.seed.seed_lengths[i])
      score = model['intercept']

      for k, v in model.items():
        if k != 'intercept':
          score += rgetattr(self, k + 's')[i] * model[k]

      self.scores.append(score)
