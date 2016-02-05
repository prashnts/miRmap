import mirmap

_mirs = mirmap.utils.load_fasta('tests/input/hsa-miR-30a-3p.fa')
_mrnas = mirmap.utils.load_fasta('tests/input/NM_024573.fa')

mim = mirmap.mm(_mrnas['NM_024573'], _mirs['hsa-miR-30a-3p'])
mim.find_potential_targets_with_seed()
mim.end_sites
mim.tgs_au
mim.prob_binomial
print(mim.report())
