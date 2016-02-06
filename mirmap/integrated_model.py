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
    self.__init_models()

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
          'tgs_au': -0.275594504153219,
          'tgs_position': 9.44582844229299e-06,
          'tgs_pairing3p': -0.0111209267382849,
          'prob_binomial': 0.0701619992923641,
          # 'cons_bls': -0.00646548621345819,
          'intercept': 0.121104869645859,
      },
      'python_only_seed7': {
          'tgs_au': -0.443606032336791,
          'tgs_position': 6.34603935320321e-05,
          'tgs_pairing3p': -0.0207672870210752,
          'prob_binomial': 0.378665477250754,
          # 'cons_bls': -0.0552713344740971,
          'intercept': 0.150015113841088,
      }
    }

  def _eval_score(self, **kwargs):
    """
    Computes the *miRmap* score(s)
    """

    self.model = model
    if model_name is not None:
        self.model = Defaults.models[model_name]
    # Reset
    self.scores = []
    # Compute
    for its in range(len(self.end_sites)):
        if self.model is None and Defaults.model_select:
            self.model = Defaults.model_select(self, its)
        score = self.model['intercept']

        for k in self.model.keys():
            if k != 'intercept':
                score += getattr(self, k+'s')[its] * self.model[k]
        self.scores.append(score)

    pass
