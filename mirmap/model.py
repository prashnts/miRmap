# -*- coding: utf-8 -*-

import warnings

from mirmap import (seed, targetscan, prob_binomial, thermodynamics,
                    evolution, spatt)
from mirmap.utils import (rgetattr, gen_dot_pipe_notation,
                          rgetattrna, rgetattrze)

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
    self.__init_spatt()
    self.__init_seed(**kwargs.get('seed_args', {}))
    self.__init_targetscan(**kwargs.get('tscan_args', {}))
    self.__init_prob_binomial(**kwargs.get('prob_args', {}))
    self.__init_thermodynamics(**kwargs.get('thermo_args', {}))
    self.__init_evolutionary(**kwargs.get('evol_args', {}))

  def __init_models(self):
    """
    Current models based on Grimson et al dataset
    """
    self.models = {
      'full_seed6': {
          '_target_scan.tgs_au': -0.275016235769136,
          '_target_scan.tgs_position': 5.43367028065211e-06,
          '_target_scan.tgs_pairing3p': -0.00233278119760994,
          '_thermodynamic.dg_duplex': 0.00772658898496047,
          '_thermodynamic.dg_binding': -0.00303683833660696,
          '_thermodynamic.dg_duplex_seed': 0.0496909801533612,
          '_thermodynamic.dg_binding_seed': -0.048931930580652,
          '_thermodynamic.dg_open': 0.000674676164622922,
          '_prob_binomial.prob_exact': 0.16111635592018,
          '_prob_binomial.prob_binomial': -0.0388333740708671,
          '_evolutionary.cons_bls': -0.00426314077593848,
          '_evolutionary.selec_phylop': -0.0112455248228072,
          'intercept': 0.148300586692704,
      },
      'full_seed7': {
          '_target_scan.tgs_au': -0.402470212080983,
          '_target_scan.tgs_position': 6.89249707831041e-05,
          '_target_scan.tgs_pairing3p': -0.0129891251446967,
          '_thermodynamic.dg_duplex': 0.0141332997802509,
          '_thermodynamic.dg_binding': -0.0132159175462755,
          '_thermodynamic.dg_duplex_seed': -0.0814445085121904,
          '_thermodynamic.dg_binding_seed': 0.115558118311931,
          '_thermodynamic.dg_open': 0.00331507347139685,
          '_prob_binomial.prob_exact': 0.792962156550929,
          '_prob_binomial.prob_binomial': -0.22119499646323,
          '_evolutionary.cons_bls': -0.0355840335642203,
          '_evolutionary.selec_phylop': -0.0127531995991629,
          'intercept': 0.349448109979275,
      },
      'python_only_seed6': {
          '_target_scan.tgs_au': -0.275594504153219,
          '_target_scan.tgs_position': 9.44582844229299e-06,
          '_target_scan.tgs_pairing3p': -0.0111209267382849,
          '_prob_binomial.prob_binomial': 0.0701619992923641,
          'intercept': 0.121104869645859,
      },
      'python_only_seed7': {
          '_target_scan.tgs_au': -0.443606032336791,
          '_target_scan.tgs_position': 6.34603935320321e-05,
          '_target_scan.tgs_pairing3p': -0.0207672870210752,
          '_prob_binomial.prob_binomial': 0.378665477250754,
          'intercept': 0.150015113841088,
      }
    }
    self.model_maps = {
      '_thermodynamic.dg_duplex':       'ΔG duplex (kcal/mol)',
      '_thermodynamic.dg_binding':      'ΔG binding (kcal/mol)',
      '_thermodynamic.dg_open':         'ΔG open (kcal/mol)',
      '_thermodynamic.dg_total':        'ΔG total (kcal/mol)',
      '_target_scan.tgs_au':            'AU content',
      '_target_scan.tgs_pairing3p':     '3\' pairing',
      '_target_scan.tgs_position':      'UTR position',
      '_target_scan.tgs_score':         'TargetScan score',
      '_prob_binomial.prob_exact':      'Probability (Exact)',
      '_prob_binomial.prob_binomial':   'Probability (Binomial)',
      '_evolutionary.cons_bl':          'Conservation (BLS)',
      '_evolutionary.selec_phylop':     'Conservation (PhyloP)',
    }

    self.display_order = [
      '_thermodynamic.dg_duplex',
      '_thermodynamic.dg_binding',
      '_thermodynamic.dg_open',
      '_thermodynamic.dg_total',
      '_target_scan.tgs_au',
      '_target_scan.tgs_pairing3p',
      '_target_scan.tgs_position',
      '_target_scan.tgs_score',
      '_prob_binomial.prob_exact',
      '_prob_binomial.prob_binomial',
      '_evolutionary.cons_bl',
      '_evolutionary.selec_phylop',
    ]

    self.model = 'full_seed'

  def __init_seed(self, **args):
    arg_init = {
      'target_seq': self.seq_mrn,
      'mirna_seq': self.seq_mir
    }
    args.update(arg_init)
    self._seed = seed.mmSeed(**args)

  def __init_targetscan(self, **args):
    self._target_scan = targetscan.mmTargetScan(self._seed, **args)

  def __init_prob_binomial(self, **args):
    if hasattr(self, '_spatt'):
      args['spatt'] = self._spatt
    self._prob_binomial = prob_binomial.mmProbBinomial(self._seed, **args)

  def __init_thermodynamics(self, **args):
    try:
      self._thermodynamic = thermodynamics.mmThermo(self._seed, **args)
    except EnvironmentError:
      warnings.warn((
        "RNAVienna not available, falling back to Python Only mode. "
        "Please Note that thermodynamic Values will NOT be available. "
      ), RuntimeWarning)
      self.model = 'python_only_seed'

  def __init_evolutionary(self, **args):
    try:
      self._evolutionary = evolution.mmEvolution(self._seed, **args)
    except EnvironmentError:
      warnings.warn((
        "PHAST not available, falling back to Python Only mode. "
        "Please Note that Evolutionary Values will NOT be available. "
      ), RuntimeWarning)
      self.model = 'python_only_seed'

  def __init_spatt(self, **args):
    try:
      self._spatt = spatt.Spatt()
    except EnvironmentError:
      warnings.warn((
        "SPATT not available, falling back to Python Only mode. "
        "Please Note that Exact Probability Value will NOT be available. "
      ), RuntimeWarning)
      self.model = 'python_only_seed'

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
    self._seed.find_potential_targets_with_seed()
    self._target_scan.routine()
    self._prob_binomial._eval_prob_binomial()
    self._prob_binomial._eval_prob_exact()
    try:
      self._thermodynamic.routine()
    except AttributeError:
      pass
    try:
      self._evolutionary.routine()
    except AttributeError:
      pass
    self._eval_score()

  def _eval_score(self):
    """
    Computes the *miRmap* score(s)
    """
    # Reset
    self.scores = []
    # Compute
    for i in range(len(self._seed.end_sites)):
      model = self.model_select(self._seed.seed_lengths[i])
      score = model['intercept']

      for k, v in model.items():
        if k != 'intercept':
          score += rgetattrze(self, k + 's')[i] * model[k]

      self.scores.append(score)
    return self.scores

  @property
  def score(self):
    try:
      return self.scores
    except AttributeError:
      return self._eval_score()

  @property
  def report(self):
    if not getattr(self, 'score', False):
      self.routine()

    report_lines = []
    for i in range(len(self._seed.end_sites)):
      end_site = self._seed.end_sites[i]
      start = max(0, end_site - self._seed.len_mirna_seq - 10)
      end = end_site + 10
      report_lines.append((
        str(start + 1) +
        ' ' * (end_site - start - len(str(start + 1)) - 1) +
        str(end_site)
      ))
      report_lines.append('|' + ' ' * (end_site - start - 2) + '|')
      report_lines.append(self._seed.target_seq[start:end])
      seed_pairing_string = gen_dot_pipe_notation(self._seed.pairings[i])
      report_lines.append((
        ' ' * (end_site - len(seed_pairing_string) - start) +
        seed_pairing_string
      ))
      report_lines.append((
        ' ' * (end_site - self._seed.len_mirna_seq - start) +
        self._seed.mirna_seq[::-1]
      ))

      model = self.model_select(self._seed.seed_lengths[i])

      for k in self.display_order:
        if k in model:
          report_lines.append(
            '  %-30s% .5f' % (self.model_maps[k], rgetattrna(self, k + 's')[i])
          )

      #FIXME: miRmap Score is NOT representative, yet.
      # report_lines.append(
      #   '  %-30s% .5f' % ("miRmap Score", self.scores[i])
      # )

    return "\n".join(report_lines)
