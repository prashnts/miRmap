import mirmap

_mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
_mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')

mim = mirmap.mm(_mrnas['NM_024573'], _mirs['hsa-miR-30a-3p'])
mim.find_potential_targets_with_seed(allowed_lengths=[6,7], allowed_gu_wobbles={6:0,7:0}, allowed_mismatches={6:0,7:0}, take_best=True)
result = mim.end_sites
