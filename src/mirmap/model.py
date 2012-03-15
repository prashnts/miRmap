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
    def eval_score(self, model_name=None, model=None):
        # Parameter: if the user doesn't give anything, trying to find the most appropriate model in Defaults
        if model_name is None:
            if model is None:
                needed_slopes = []
                for feat in ['tgs_aus', 'tgs_positions', 'tgs_pairing3ps', 'dg_duplexs', 'dg_bindings', 'dg_opens', 'prob_exacts', 'prob_binomials', 'cons_blss', 'selec_phylops']:
                    if hasattr(self, feat):
                        needed_slopes.append('slope_'+feat[:-1])
                num_feats_largest_model = 0
                for model_name,model_params in Defaults.models.items():
                    if len(set(model_params.keys()).intersection(set(needed_slopes))) == (len(model_params) - 1):
                        if num_feats_largest_model < len(model_params.keys()):
                            num_feats_largest_model = len(model_params.keys())
                            #print(model_name)
                            self.model = model_params
            else:
                self.model = model
        else:
            self.model = Defaults.models[model_name]
        # Reset
        self.scores = []
        # Compute
        for its in range(len(self.end_sites)):
            score = self.model['intercept']
            if self.model.has_key('slope_tgs_au'):
                score += self.tgs_aus[its] * self.model['slope_tgs_au']
            if self.model.has_key('slope_tgs_position'):
                score += self.tgs_positions[its] * self.model['slope_tgs_position']
            if self.model.has_key('slope_tgs_pairing3p'):
                score += self.tgs_pairing3ps[its] * self.model['slope_tgs_pairing3p']
            if self.model.has_key('slope_dg_duplex'):
                score += self.dg_duplexs[its] * self.model['slope_dg_duplex']
            if self.model.has_key('slope_dg_binding'):
                score += self.dg_bindings[its] * self.model['slope_dg_binding']
            if self.model.has_key('slope_dg_open'):
                score += self.dg_opens[its] * self.model['slope_dg_open']
            if self.model.has_key('slope_prob_exact'):
                score += self.prob_exacts[its] * self.model['slope_prob_exact']
            if self.model.has_key('slope_prob_binomial'):
                score += self.prob_binomials[its] * self.model['slope_prob_binomial']
            if self.model.has_key('slope_cons_bls'):
                score += self.cons_blss[its] * self.model['slope_cons_bls']
            if self.model.has_key('slope_selec_phylop'):
                score += self.selec_phylops[its] * self.model['slope_selec_phylop']
            self.scores.append(score)

class Defaults(object):
    models = {}
    # All features model
    models['full'] = {}
    models['full']['slope_tgs_au'] = -0.4119
    models['full']['slope_tgs_position'] = 5.088e-05
    models['full']['slope_tgs_pairing3p'] = -0.01050
    models['full']['slope_dg_duplex'] = 0.009900
    models['full']['slope_dg_binding'] = -0.002919
    models['full']['slope_dg_open'] = 0.003914
    models['full']['slope_prob_exact'] = 1.022
    models['full']['slope_prob_binomial'] = -0.3582
    models['full']['slope_cons_bls'] = -0.03431
    models['full']['slope_selec_phylop'] = -0.009622
    models['full']['intercept'] = 0.08088
    # 7-main features model
    models['7main'] = {}
    models['7main']['slope_tgs_au'] = -3.701e-01
    models['7main']['slope_tgs_position'] = 3.931e-05
    models['7main']['slope_tgs_pairing3p'] = -2.175e-02
    models['7main']['slope_dg_open'] = 3.384e-03
    models['7main']['slope_prob_exact'] = 5.683e-01
    models['7main']['slope_cons_bls'] = -3.526e-02
    models['7main']['slope_selec_phylop'] = -8.838e-03
    models['7main']['intercept'] = 1.058e-02
    # 7-main features model minus the conservation features
    models['7main_minus_conservation'] = {}
    models['7main_minus_conservation']['slope_tgs_au'] = -4.246e-01
    models['7main_minus_conservation']['slope_tgs_position'] = 5.083e-05
    models['7main_minus_conservation']['slope_tgs_pairing3p'] = -2.350e-02
    models['7main_minus_conservation']['slope_dg_open'] = 3.812e-03
    models['7main_minus_conservation']['slope_prob_exact'] = 4.856e-01
    models['7main_minus_conservation']['intercept'] = -2.872e-02
    # TargetScan model
    models['targetscan'] = {}
    models['targetscan']['slope_tgs_au'] = -0.5368
    models['targetscan']['slope_tgs_position'] = 0.0001234
    models['targetscan']['slope_tgs_pairing3p'] = -0.02573
    models['targetscan']['intercept'] = 0.1586
