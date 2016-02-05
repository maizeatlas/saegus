


import simuOpt
simuOpt.setOptions(alleleType='short', optimized=True, quiet=True)
import simuPOP as sim
import pandas as pd
import collections as col
from wgs import breed, operators, selection, helpers, parser, parameterizer, selection
import random
import numpy as np
np.set_printoptions(suppress=True, precision=3)
import matplotlib.pyplot as plt



hapmap = pd.read_csv('clean_hapmap.txt')
genetic_map = hapmap.ix[:, :'cM_pos']
genetic_map = pd.read_csv('nam_prefounders_genetic_map.txt', index_col=None,
                         sep='\t')

chr_cm_positions = col.OrderedDict()
for i in range(1, 11):
    chr_cm_positions[i] = []

for idx in range(len(genetic_map)):
    chr_cm_positions[int(genetic_map.iloc[idx]['chr'])].append(
    float(genetic_map.iloc[idx]['cM_pos']))


cM_positions = []
for k, v in chr_cm_positions.items():
    cM_positions.append(v)


snp_to_integer = {'A':0, 'C':1, 'G':2, 'T':3, '-':4, '+':5}

integral_valued_loci = []
relative_integral_valued_loci = {}
for idx in range(len(genetic_map)):
    if str(genetic_map.iloc[idx]['cM_pos'])[-2:] == '.0':
        integral_valued_loci.append(idx)
        relative_integral_valued_loci[idx] = (genetic_map.iloc[idx]['chr'], genetic_map.iloc[idx]['cM_pos'])

alleles = {i: (snp_to_integer[hapmap.ix[i, 'alleles'][0]], 
               snp_to_integer[hapmap.ix[i, 'alleles'][-1]]) for i in
          range(len(hapmap))}

recombination_rates = []
for chromosome in cM_positions:
    for cM in chromosome:
        if str(cM)[-2:] == '.6':
            recombination_rates.append(0.01)
        else:
            recombination_rates.append(0.0)

allele_names = ['A', 'C', 'T', 'G', 'D', 'I']

flat_cM_positions = []
for cMs in cM_positions:
    flat_cM_positions.extend(cMs)


nam = sim.loadPopulation('nam_prefounders.pop')
nam.setSubPopName('prefounders', 0)
sample_sizes = {i: 100 for i in range(0, 21, 2)}

sim_params = {
                'gens_selection': 10,
                'gens_random_mating': 3,
                'main_pop_size': 2000,
                'proportion_saved': 0.05,
                'overshoot': 0.50,
                'breeding_inds_per_sp': 5,
                'heritability': 0.7,
                'sample_sizes': sample_sizes,
                'replicates': 1,
                'prefounder_file': 'nam_prefounders.pop',
                'allele_effects_file': '',
                'qtl': 20,
                'allele_effects': 1,
                'founders': [(1,5), (1, 8), (1, 4), (1, 11)],
}

population_statistics = {
    'aggregate': {},
    'selected': {},
    'non-selected': {}
}



meta_population_statistics = {
    'aggregate': {},
    'selected': {},
    'non-selected': {}
}

s = selection.Truncation(sim_params['gens_selection'],
                       sim_params['gens_random_mating'],
                       sim_params['main_pop_size'],
                       sim_params['proportion_saved'],
                       sim_params['overshoot'],
                       sim_params['breeding_inds_per_sp'],
                       sim_params['heritability'],
                       sim_params['sample_sizes'],
                       sim_params['replicates'])


sim.tagID(nam, reset=True)

founders = sim_params['founders']
replicated_nam = sim.Simulator(nam, rep=2)
pop = replicated_nam.extract(0)
pop.dvars().statistics = population_statistics
meta = replicated_nam.extract(0)
#meta.removeSubPops(0)


# ### Simulated Breeding Scenario ###

# In[ ]:

s.generate_f_one(pop, recombination_rates, sim_params['founders'])
s.expand_by_selfing(pop, recombination_rates)
s.mate_and_merge(pop, recombination_rates)
s.interim_random_mating(pop, recombination_rates)

sim.stat(pop, numOfSegSites=integral_valued_loci, vars=['numOfSegSites', 'segSites'])


# ## Choose QTL and Assign Effects ##

# In[ ]:

qtl = parameterizer.seg_qtl_chooser(pop, integral_valued_loci, sim_params['qtl'])

triplet_qtl = []
for locus in qtl:
    triplet_qtl.append(locus-1)
    triplet_qtl.append(locus)
    triplet_qtl.append(locus+1)
triplet_qtl = sorted(triplet_qtl)


allele_effects = {}
for tqtl in triplet_qtl:
    for allele in alleles[tqtl]:
        allele_effects[tqtl, allele] = random.expovariate(sim_params['allele_effects'])


        

pop.dvars().qtl = qtl
pop.dvars().triplet_qtl = triplet_qtl
pop.dvars().allele_effects = allele_effects


# In[ ]:

s.recurrent_truncation_selection(pop, meta, triplet_qtl, allele_effects,
                                recombination_rates)
                                
meta.removeSubPops(0)

qtd = helpers.Frq(meta, triplet_qtl, alleles, allele_effects)

af = qtd.allele_frequencies(meta, range(meta.totNumLoci()))
qtalleles = qtd.rank_allele_effects(meta, triplet_qtl, alleles, allele_effects)
ties = [locus for locus in range(meta.totNumLoci()) if af['minor', 'alleles'][locus] == af['major', 'alleles'][locus]]
for t in ties:
    af['major', 'alleles'][t] = list(meta.dvars().alleleFreq[t])[0]
    af['minor', 'alleles'][t] = list(meta.dvars().alleleFreq[t])[1]
sum(np.equal(list(af['minor', 'alleles'].values()), list(af['major', 'alleles'].values())))


# ## Gather Data for Use in GWAS ##

# In[ ]:

pca = helpers.PCA(meta, range(meta.totNumLoci()), af)
minor_ac = pca.calculate_count_matrix(meta, af['minor', 'alleles'], 'sim_minor_allele_count.txt')
eigendata = pca.svd(meta, minor_ac)
ts = pca.test_statistic(meta, eigendata['values'])

integer_to_snp = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '-', 5: '+'}
raw_hmap = pd.read_csv('hapmap3.txt', delimiter='\t', index_col=0)
locus_names = list(raw_hmap['nearest.site'])
pos_column = list(raw_hmap['agp_pos'])
individual_names = {ind.ind_id: 'RS_R'+str(1)+'_G'+str(int(ind.generation)) + '_I'+str(int(ind.ind_id))
                   for ind in meta.individuals()}

cols_for_hapmap = {'locus_names': locus_names, 'pos_column': pos_column}
                   
gwas = helpers.GWAS(meta, individual_names, locus_names, pos_column)
hmap = gwas.hapmap_formatter(integer_to_snp, 'C:\\GWAS\\sim_hapmap.txt')
popstruct = gwas.population_structure_formatter(eigendata, 'C:\\GWAS\\sim_structure.txt')
phenos = gwas.trait_formatter('C:\\GWAS\\input\\sim_trait_vector.txt')
kinship_matrix = gwas.calc_kinship_matrix(minor_ac, af, 'C:\\GWAS\\sim_kinship.txt')


