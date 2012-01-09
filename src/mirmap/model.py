# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""miRmap module"""

import seed

class mmModel(seed.mmSeed):
    def eval_score(self, slope_tgs_au=None, slope_tgs_position=None, slope_tgs_pairing3p=None, slope_dg_duplex=None, slope_dg_binding=None, slope_dg_open=None, slope_prob_exact=None, slope_prob_binomial=None, slope_bls=None, slope_phylop=None, intercept=None):
        # Parameters: if the user doesn't give anything, trying to find the most appropriate model in Defaults
        if slope_tgs_au is None and slope_tgs_position is None and slope_tgs_pairing3p is None and slope_dg_duplex is None and slope_dg_binding is None and slope_dg_open is None and slope_prob_exact is None and slope_prob_binomial is None and slope_bls is None and slope_phylop is None and intercept is None:
            needed_slopes = []
            for feat in ['tgs_aus', 'tgs_positions', 'tgs_pairing3ps', 'dg_duplexs', 'dg_bindings', 'dg_opens', 'prob_exacts', 'prob_binomials', 'blss', 'phylops']:
                if hasattr(self, feat):
                    needed_slopes.append('slope_'+feat[:-1])
            num_feats_largest_model = 0
            for model_name,model in Defaults.models.items():
                if len(set(model.keys()).intersection(set(needed_slopes))) == (len(model) - 1):
                    if num_feats_largest_model < len(model.keys()):
                        num_feats_largest_model = len(model.keys())
            if hasattr(self, 'tgs_aus'):
                slope_tgs_au = model=['slope_tgs_au']
            if hasattr(self, 'tgs_positions'):
                slope_tgs_position = model=['slope_tgs_position']
            if hasattr(self, 'tgs_pairing3ps'):
                slope_tgs_pairing3p = model=['slope_tgs_pairing3p']
            if hasattr(self, 'dg_duplexs'):
                slope_dg_duplex = model=['slope_dg_duplex']
            if hasattr(self, 'dg_bindings'):
                slope_dg_binding = model=['slope_dg_binding']
            if hasattr(self, 'dg_opens'):
                slope_dg_open = model=['slope_dg_open']
            if hasattr(self, 'prob_exacts'):
                slope_prob_exact = model=['slope_prob_exact']
            if hasattr(self, 'prob_binomials'):
                slope_prob_binomial = model=['slope_prob_binomial']
            if hasattr(self, 'blss'):
                slope_blss = model=['slope_blss']
            if hasattr(self, 'phylops'):
                slope_phylops = model=['slope_phylops']
        # Reset
        self.scores = []
        # Compute
        for its in range(len(self.end_sites)):
            score = intercept
            if hasattr(self, 'tgs_aus') and slope_tgs_au is not None:
                score += self.tgs_aus[its] * slope_tgs_au
            if hasattr(self, 'tgs_positions') and slope_tgs_position is not None:
                score += self.tgs_positions[its] * slope_tgs_position
            if hasattr(self, 'tgs_pairing3ps') and slope_tgs_pairing3p is not None:
                score += self.tgs_pairing3ps[its] * slope_tgs_pairing3p
            if hasattr(self, 'dg_duplexs') and slope_dg_duplex is not None:
                score += self.dg_duplexes[its] * slope_dg_duplex
            if hasattr(self, 'dg_bindings') and slope_dg_binding is not None:
                score += self.dg_bindings[its] * slope_dg_binding
            if hasattr(self, 'dg_opens') and slope_dg_open is not None :
                score += self.dg_opens[its] * slope_dg_open
            if hasattr(self, 'prob_exacts') and slope_prob_exact is not None:
                score += self.prob_exacts[its] * slope_prob_exact
            if hasattr(self, 'prob_binomials') and slope_prob_binomial is not None:
                score += self.prob_binomials[its] * slope_prob_binomial
            if hasattr(self, 'blss') and slope_blss is not None:
                score += self.blss[its] * slope_blss
            if hasattr(self, 'phylops') and slope_phylops is not None:
                score += self.phylops[its] * slope_phylops
            self.scores.append(score)

class Defaults(object):
    models = {}
    # Full model (all features)
    models['full'] = {}
    models['full']['slope_tgs_au'] = -3.701e-01
    models['full']['slope_tgs_position'] = 3.931e-05
    models['full']['slope_tgs_pairing3p'] = -2.175e-02
    models['full']['slope_dg_open'] = 3.384e-03
    models['full']['slope_prob_exact'] = 5.683e-01
    models['full']['slope_bls'] = -3.526e-02
    models['full']['slope_phylop'] = -8.838e-03
    models['full']['intercept'] = 1.058e-02
    # Full model minus the conservation features
    models['full_minus_conservation'] = {}
    models['full_minus_conservation']['slope_tgs_au'] = -4.246e-01
    models['full_minus_conservation']['slope_tgs_position'] = 5.083e-05
    models['full_minus_conservation']['slope_tgs_pairing3p'] = -2.350e-02
    models['full_minus_conservation']['slope_dg_open'] = 3.812e-03
    models['full_minus_conservation']['slope_prob_exact'] = 4.856e-01
    models['full_minus_conservation']['intercept'] = -2.872e-02
