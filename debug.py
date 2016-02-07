from mirmap import miRmap

seq_target = (
  'GCUACAGUUUUUAUUUAGCAUGGGGAUUGCAGAGUGACCAGCAC'
  'ACUGGACUCCGAGGUGGUUCAGACAAGACAGAGGGGAGCAGUGG'
  'CCAUCAUCCUCCCGCCAGGAGCUUCUUCGUUCCUGCGCAUAUAG'
  'ACUGUACAUUAUGAAGAAUACCCAGGAAGACUUUGUGACUGUCA'
  'CUUGCUGCUUUUUCUGCGCUUCAGUAACAAGUGUUGGCAAACGA'
  'GACUUUCUCCUGGCCCCUGCCUGCUGGAGAUCAGCAUGCCUGUC'
  'CUUUCAGUCUGAUCCAUCCAUCUCUCUCUUGCCUGAGGGGAAAG'
  'AGAGAUGGGCCAGGCAGAGAACAGAACUGGAGGCAGUCCAUCUA'
)

seq_mirna = 'UAGCAGCACGUAAAUAUUGGCG'

ini = {
  'seq_mir': seq_mirna,
  'seq_mrn': seq_target,
  'seed_args': {
    'allowed_lengths': [6, 7],
    'allowed_gu_wobbles': {6: 0, 7: 0},
    'allowed_mismatches': {6: 0, 7: 0},
    'take_best': True
  },
  'tscan_args': {
    'with_correction': False
  },
  'prob_args': {
  }
}

obj = miRmap(**ini)

print(obj.report)
