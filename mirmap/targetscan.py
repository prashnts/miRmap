# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2013 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""TargetScan features."""

import operator

from collections import namedtuple

try:
    operator.div = operator.truediv
except AttributeError:
    pass


def align_helper(utr_3p_seq, mir_3p_seq, mir_offset, utr_offset, overhang):
    """
    Helper Function.
    """

    score = 0
    tempscore = 0
    prevmatch = 0
    i = 0
    offset = max(mir_offset, utr_offset)
    cutoff = min(len(mir_3p_seq) - mir_offset, len(utr_3p_seq) - utr_offset)
    while i < cutoff:

        m1 = i + utr_offset
        m2 = i + mir_offset

        if utr_3p_seq[m1] + mir_3p_seq[m2] in ['AU', 'UA', 'GC', 'CG']:
            if prevmatch == 0:
                tempscore = 0
            tempscore += 1 if (4 <= i + mir_offset - overhang <= 7) else 0.5
            prevmatch += 1
        elif prevmatch >= 2:
            if tempscore > score:
                score = tempscore
            tempscore = 0
            prevmatch = 0
        else:
            tempscore = 0
            prevmatch = 0
        i += 1
    if prevmatch >= 2 and tempscore > score:
        score = tempscore
    score -= max(0, ((offset - 2) / 2.0))
    return score


class mmTargetScan(object):
    """
    miRmap TargetScan.

    Args:
       ts_types (object): Parameters by seed-type.
       ca_window_length (int): Sequence length to compute the score with.
       with_correction (bool): Apply the linear regression correction?
    """

    def __init__(self, seed, **kwargs):
        self.seed = seed
        self.__init_defaults()
        self.__init_args(**kwargs)

    def __init_args(self, **kwargs):
        allowed_args = [
            'ca_window_length',
            'with_correction',
        ]
        arg = {k: v for k, v in kwargs.items() if k in allowed_args}
        self.__dict__.update(arg)

    def __init_defaults(self):
        ts_types_builder = namedtuple('TSTypes', [
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
        ])

        ts_t_default = {
            '6mer': {
                'name': '6mer',
                'up_shift': -8,
                'down_shift': 0,
                'fc_mean': -0.015,
                'ca_fc_slope': -0.241,
                'ca_fc_intercept': 0.115,
                'ca_weights_up': [
                  31.0, 30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0,
                  21.0, 20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0,
                  11.0, 10.0,  9.0,  8.0,  7.0,  6.0,  5.0,  4.0,  3.0,  2.0
                ],
                'ca_weights_down': [
                   2.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0,
                  11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                  21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0
                ],
                'po_fc_slope': 0.000049,
                'po_fc_intercept': -0.033,
                'pa_fc_slope': -0.00278,
                'pa_fc_intercept': -0.0091,
                'pa_mirna_seed_start': 7,
                'pa_mirna_seed_overhang': 1,
            },
            '7mer-A1': {
                'name': '7mer-A1',
                'up_shift': -8,
                'down_shift': 1,
                'fc_mean': -0.099,
                'ca_fc_slope': -0.42,
                'ca_fc_intercept': 0.137,
                'ca_weights_up': [
                  31.0, 30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0,
                  21.0, 20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0,
                  11.0, 10.0,  9.0,  8.0,  7.0,  6.0,  5.0,  4.0,  3.0,  2.0
                ],
                'ca_weights_down': [
                   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0,
                  12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0,
                  22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0
                ],
                'po_fc_slope': 0.000072,
                'po_fc_intercept': -0.131,
                'pa_fc_slope': -0.0211,
                'pa_fc_intercept': -0.053,
                'pa_mirna_seed_start': 7,
                'pa_mirna_seed_overhang': 1,
            },
            '7mer-m8': {
                'name': '7mer-m8',
                'up_shift': -8,
                'down_shift': 0,
                'fc_mean': -0.161,
                'ca_fc_slope': -0.5,
                'ca_fc_intercept': 0.108,
                'ca_weights_up': [
                  30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0,
                  20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0,
                  10.0,  9.0,  8.0,  7.0,  6.0,  5.0,  4.0,  3.0,  2.0,  1.0
                ],
                'ca_weights_down': [
                   2.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0,
                  11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                  21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0
                ],
                'po_fc_slope': 0.000091,
                'po_fc_intercept': -0.198,
                'pa_fc_slope': -0.031,
                'pa_fc_intercept': -0.094,
                'pa_mirna_seed_start': 8,
                'pa_mirna_seed_overhang': 0,
            },
            '8mer': {
                'name': '8mer',
                'up_shift': -8,
                'down_shift': 1,
                'fc_mean': -0.31,
                'ca_fc_slope': -0.64,
                'ca_fc_intercept': 0.055,
                'ca_weights_up': [
                  30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0,
                  20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0,
                  10.0,  9.0,  8.0,  7.0,  6.0,  5.0,  4.0,  3.0,  2.0,  1.0
                ],
                'ca_weights_down': [
                   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0,
                  12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0,
                  22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0
                ],
                'po_fc_slope': 0.000172,
                'po_fc_intercept': -0.38,
                'pa_fc_slope': -0.0041,
                'pa_fc_intercept': -0.299,
                'pa_mirna_seed_start': 8,
                'pa_mirna_seed_overhang': 0,
            }
        }
        ts_types = {k: ts_types_builder(**v) for k, v in ts_t_default.items()}

        self.__dict__.update({
            'with_correction': True,
            'ca_window_length': 30,
            'ts_types': ts_types
        })

    def _targetscan_ts_type(self, seed_length, nt1):
        """
        Returns the target site type according to TargetScan terminology.
        """
        if seed_length == 6:
            return '8mer' if nt1 == 'A' else '7mer-m8'
        elif seed_length >= 7:
            return '7mer-A1' if nt1 == 'A' else '6mer'

        raise ValueError("seed_length should be >= 6.")

    def _eval_tgs_au(self):
        """
        Computes the *AU content* score.
        """

        # Reset
        self.tgs_aus = []
        # Helper function
        binarize = lambda x: 1.0 if x in ['U', 'A'] else 0.0

        # Compute
        for its in range(len(self.seed.end_sites)):
            end_site = self.seed.end_sites[its]
            try:
                ts_type = self._targetscan_ts_type(
                    self.seed.seed_lengths[its],
                    self.seed.target_seq[end_site - 1]
                )
                tts = self.ts_types[ts_type]

                sus = max(
                    0,
                    end_site + tts.up_shift - self.ca_window_length
                )
                sue = end_site + tts.up_shift
                seq_up = self.seed.target_seq[sus:sue]

                sds = end_site - 1 + tts.down_shift
                sde = min(
                    self.seed.len_target_seq,
                    end_site + self.ca_window_length - 1 + tts.down_shift
                )
                seq_down = self.seed.target_seq[sds:sde]

                wup = tts.ca_weights_up[len(tts.ca_weights_up) - len(seq_up):]
                wdn = tts.ca_weights_down[:len(seq_down)]

                content = sum([
                    sum(map(operator.div, map(binarize, seq_up), wup)),
                    sum(map(operator.div, map(binarize, seq_down), wdn))
                ])

                content /= sum([
                    sum(map(operator.div, [1.0] * len(wup), wup)),
                    sum(map(operator.div, [1.0] * len(wdn), wdn))
                ])

                if self.with_correction:
                    self.tgs_aus.append(sum([
                        content * tts.ca_fc_slope,
                        tts.ca_fc_intercept - tts.fc_mean
                    ]))
                else:
                    self.tgs_aus.append(content)
            except ValueError:
                self.tgs_aus.append(None)
        return self.tgs_aus

    def _eval_tgs_position(self):
        """
        Computes the *UTR position* score.
        """

        # Reset
        self.tgs_positions = []
        # Compute
        for its in range(len(self.seed.end_sites)):
            end_site = self.seed.end_sites[its]
            try:
                ts_type = self._targetscan_ts_type(
                    self.seed.seed_lengths[its],
                    self.seed.target_seq[end_site - 1]
                )
                tts = self.ts_types[ts_type]
                closest_term = max(
                    min(
                        end_site + tts.up_shift,
                        self.seed.len_target_seq - end_site + tts.down_shift
                    ),
                    1500
                )

                if self.with_correction:
                    self.tgs_positions.append(sum([
                        float(closest_term) * tts.po_fc_slope,
                        tts.po_fc_intercept - tts.fc_mean
                    ]))
                else:
                    self.tgs_positions.append(float(closest_term))
            except ValueError:
                self.tgs_positions.append(None)
        return self.tgs_positions

    def _eval_tgs_pairing3p(self):
        """
        Computes the *3' pairing* score.
        """

        # Reset
        self.tgs_pairing3ps = []
        # Compute
        for its in range(len(self.seed.end_sites)):
            end_site = self.seed.end_sites[its]
            try:
                ts_type = self._targetscan_ts_type(
                    self.seed.seed_lengths[its],
                    self.seed.target_seq[end_site - 1]
                )
                tts = self.ts_types[ts_type]

                uts = max(0, end_site - tts.pa_mirna_seed_start - 15)
                ute = end_site - tts.pa_mirna_seed_start
                utr_3p_seq = self.seed.target_seq[uts:ute][::-1]
                mir_3p_seq = self.seed.mirna_seq[tts.pa_mirna_seed_start:]
                maxscore = max(len(utr_3p_seq), len(mir_3p_seq))
                scores_mir = []
                scores_utr = []
                for offset in range(maxscore):
                    scores_mir.append(align_helper(
                        utr_3p_seq, mir_3p_seq, offset, 0,
                        tts.pa_mirna_seed_overhang)
                    )
                    scores_utr.append(align_helper(
                        utr_3p_seq, mir_3p_seq, 0, offset,
                        tts.pa_mirna_seed_overhang)
                    )
                if self.with_correction:
                    self.tgs_pairing3ps.append(sum([
                        float(max(scores_mir + scores_utr)) * tts.pa_fc_slope,
                        tts.pa_fc_intercept - tts.fc_mean
                    ]))
                else:
                    self.tgs_pairing3ps.append(
                        float(max(scores_mir + scores_utr))
                    )
            except ValueError:
                self.tgs_pairing3ps.append(None)
        return self.tgs_pairing3ps

    def _eval_tgs_score(self):
        """
        Computes the *TargetScan* score combining *AU content*, *UTR position*
        and *3' pairing* scores.
        """

        # Reset
        self.tgs_scores = []
        # Compute
        for its in range(len(self.seed.end_sites)):
            end_site = self.seed.end_sites[its]
            try:
                ts_type = self._targetscan_ts_type(
                    self.seed.seed_lengths[its],
                    self.seed.target_seq[end_site - 1]
                )
                tts = self.ts_types[ts_type]
                if self.with_correction:
                    self.tgs_scores.append(sum([
                        self.tgs_aus[its],
                        self.tgs_positions[its],
                        self.tgs_pairing3ps[its],
                        tts.fc_mean]
                    ))
                else:
                    self.tgs_scores.append(sum([
                        self.tgs_aus[its],
                        self.tgs_positions[its],
                        self.tgs_pairing3ps[its]
                    ]))
            except AttributeError as e:
                if 'mmTargetScan' in str(e):
                    raise AttributeError("Routine Did not Run.")
                self.tgs_scores.append(None)
        return self.tgs_scores

    def routine(self):
        self._eval_tgs_au()
        self._eval_tgs_position()
        self._eval_tgs_pairing3p()
        self._eval_tgs_score()

    @property
    def tgs_au(self):
        """*AU content* score with default parameters."""
        try:
            tgs_au = max(self.tgs_aus)
        except AttributeError:
            self._eval_tgs_au()
            tgs_au = max(self.tgs_aus)
        finally:
            return tgs_au

    @property
    def tgs_position(self):
        """*UTR position* score with default parameters."""
        try:
            tgs_positions = max(self.tgs_positions)
        except AttributeError:
            self._eval_tgs_position()
            tgs_positions = max(self.tgs_positions)
        finally:
            return tgs_positions

    @property
    def tgs_pairing3p(self):
        """*3' pairing* score with default parameters."""
        try:
            tgs_pairing3ps = max(self.tgs_pairing3ps)
        except AttributeError:
            self._eval_tgs_pairing3p()
            tgs_pairing3ps = max(self.tgs_pairing3ps)
        finally:
            return tgs_pairing3ps

    @property
    def tgs_score(self):
        """*TargetScan* score with default parameters."""
        try:
            tgs_scores = max(self.tgs_scores)
        except AttributeError:
            self.routine()
            tgs_scores = max(self.tgs_scores)
        finally:
            return tgs_scores
