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
    def eval_score(self, model=None):
        # Parameter: if the user doesn't give anything, trying to find the most appropriate model in Defaults
        if model is None:
            needed_slopes = []
            for feat in ['tgs_aus', 'tgs_positions', 'tgs_pairing3ps', 'dg_duplexs', 'dg_bindings', 'dg_opens', 'prob_exacts', 'prob_binomials', 'blss', 'phylops']:
                if hasattr(self, feat):
                    needed_slopes.append('slope_'+feat[:-1])
            num_feats_largest_model = 0
            for model_name,model_params in Defaults.models.items():
                if len(set(model_params.keys()).intersection(set(needed_slopes))) == (len(model_params) - 1):
                    if num_feats_largest_model < len(model_params.keys()):
                        num_feats_largest_model = len(model_params.keys())
                        #print(model_name)
                        self.model = model_params
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
                score += self.dg_duplexes[its] * self.model['slope_dg_duplex']
            if self.model.has_key('slope_dg_binding'):
                score += self.dg_bindings[its] * self.model['slope_dg_binding']
            if self.model.has_key('slope_dg_open'):
                score += self.dg_opens[its] * self.model['slope_dg_open']
            if self.model.has_key('slope_prob_exact'):
                score += self.prob_exacts[its] * self.model['slope_prob_exact']
            if self.model.has_key('slope_prob_binomial'):
                score += self.prob_binomials[its] * self.model['slope_prob_binomial']
            if self.model.has_key('slope_bls'):
                score += self.blss[its] * self.model['slope_bls']
            if self.model.has_key('slope_phylop'):
                score += self.phylops[its] * self.model['slope_phylop']
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
