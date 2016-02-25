
# coding: utf-8

# ### Short Description
# 
#     IPYNB for creating output and developing methods to analyze the output from GWAS.
#     This ipynb will choose the same parameters for many replicates and compare them.
#     The prior version Chose different qtl and allele effects for every replicate.

# # Identical QTL Parameters
#     
# ### Twenty Loci
#     
#     This simulates 20 replicates of recurrent selection.
#     Each replicate has the same QTL and allele effects.
#     

# ### Generating Data for GWAS with TASSEL

# In[ ]:

import simuOpt
simuOpt.setOptions(alleleType='short', optimized=True, numThreads=4, quiet=True)
import simuPOP as sim
import pandas as pd
import collections as col
from saegus import breed, operators, simulate, analyze, parse, parameters
import random
import copy
import yaml
import numpy as np
np.set_printoptions(suppress=True, precision=3)


# In[ ]:

hapmap = pd.read_csv('clean_hapmap.txt')
genetic_map = hapmap.ix[:, :'cM_pos']
genetic_map = pd.read_csv('nam_prefounders_genetic_map.txt', index_col=None,
                         sep='\t')

#raw_hmap = pd.read_csv('hapmap3.txt', delimiter='\t', index_col=0)

with open('dummy_columns_in_hapmap.yaml', 'r') as dcols:
    dummy_colz = yaml.load(dcols)
    
locus_names = dummy_colz['locus_names']
pos_column = dummy_colz['pos_column']


chr_cM_positions = {}
for i in range(1, 11):
    chr_cM_positions[i] = []

for idx in range(len(genetic_map)):
    chrome = str(int())
    chr_cM_positions[int(genetic_map.iloc[idx]['chr'])].append(genetic_map.iloc[idx]['cM_pos'])


cM_positions = []
for i in range(1, 11):
    cM_positions.append(chr_cM_positions[i])


# In[ ]:

allele_names = ['A', 'C', 'T', 'G', 'D', 'I']
snp_to_integer = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-':4, '+':5}
integer_to_snp = {0: 'A', 1:'C', 2: 'G', 3: 'T', 4: '-', 5: '+'}


integral_valued_loci = []
relative_integral_valued_loci = {}
for idx in range(len(genetic_map)):
    if str(genetic_map.iloc[idx]['cM_pos'])[-2:] == '.0':
        integral_valued_loci.append(idx)
        relative_integral_valued_loci[idx] = (genetic_map.iloc[idx]['chr'], genetic_map.iloc[idx]['cM_pos'])

alleles = {i: [snp_to_integer[hapmap.ix[i, 'alleles'][0]], 
               snp_to_integer[hapmap.ix[i, 'alleles'][-1]]] for i in
          range(len(hapmap))}

recombination_rates = []
for chromosome in cM_positions:
    for cM in chromosome:
        if str(cM)[-2:] == '.6':
            recombination_rates.append(0.01)
        else:
            recombination_rates.append(0.0)

flat_cM_positions = []
for cMs in cM_positions:
    flat_cM_positions.extend(cMs)


# In[ ]:

nam = sim.loadPopulation('nam_prefounders.pop')
sim.tagID(nam, reset=True)
nam.setSubPopName('prefounders', 0)
sample_sizes = {i: 100 for i in range(0, 21, 2)}
locus_names = list(range(nam.totNumLoci()))


genetic_structure = {}
#genetic_structure['cM_positions'] = cM_positions
#enetic_structure['chr_cM_positions'] = chr_cM_positions
genetic_structure['allele_names'] = allele_names
genetic_structure['integral_valued_loci'] = integral_valued_loci
genetic_structure['relative_integral_valued_loci'] = relative_integral_valued_loci
genetic_structure['alleles'] = alleles
genetic_structure['recombination_rates'] = recombination_rates


# In[ ]:

sim_params = {
                'generations_of_selection': 10,
                'generations_of_drift': 10,
                'generations_of_random_mating': 3,
                'number_of_replicates': 20,
                'operating_population_size': 200,
                'proportion_of_individuals_saved': 0.05,
                'overshoot_as_proportion': 0.50,
                'individuals_per_breeding_subpop': 5,
                'heritability': 0.7,
                'meta_pop_sample_sizes': sample_sizes,
                'prefounder_file_name': 'nam_prefounders.pop',
                'founders': [(3,18), (2, 13), (7, 14), (1, 19),
                            (14, 17), (1, 20), (17, 21), (9, 22)]
    }


# In[ ]:

qtl_params = {
                'qtl': 20,
                'allele_effects': 1,
}
selection_statistics = {
    'aggregate': {},
    'selected': {},
    'non-selected': {}
}
drift_statistics = {
    'aggregate': {},
    'selected': {},
    'non-selected': {}
}


# In[ ]:

ind_names_for_gwas = {i: {} for i in range(sim_params['number_of_replicates'])}


# In[ ]:

s = simulate.Truncation(sim_params['generations_of_selection'],
                       sim_params['generations_of_random_mating'],
                       sim_params['operating_population_size'],
                       sim_params['proportion_of_individuals_saved'],
                       sim_params['overshoot_as_proportion'],
                       sim_params['individuals_per_breeding_subpop'],
                       sim_params['heritability'],
                       sim_params['meta_pop_sample_sizes'],
                       sim_params['number_of_replicates'])


# In[ ]:

founders = sim_params['founders']
replicated_nam = sim.Simulator(nam, rep=3, stealPops=False)
pop = replicated_nam.extract(0)


# ### Run MAGIC Mating Scheme

# In[ ]:

s.generate_f_one(pop, recombination_rates, sim_params['founders'])
s.recombinatorial_convergence(pop, recombination_rates)
s.expand_by_selfing(pop, recombination_rates)
s.interim_random_mating(pop, recombination_rates)


# ## Adapting QTL and Allele Effects to Multiple Replicate Case

# In[ ]:

multipop = sim.Simulator(pop, sim_params['number_of_replicates'])
multi_meta = sim.Simulator(nam, sim_params['number_of_replicates'], stealPops=False)


# #### Assign Each Replicate Identical Parameters
#     Determines a single random set of QTL/allele effects and assigns the
#     same information to every replicate.

# In[ ]:

triplet_qtl, allele_effects = parameters.assign_identical_qtl_parameters(multipop, alleles, integral_valued_loci, 20, [1])


# In[ ]:

assert type(triplet_qtl[0]) == type([]), "Variables are flip-flopped in return."


# In[ ]:

for repid, pop_rep in enumerate(multipop.populations()):
    pop_rep.dvars().statistics = copy.deepcopy(selection_statistics)


# In[ ]:

s.replicate_selection(multipop, multi_meta, triplet_qtl, allele_effects,
                                recombination_rates)


# In[ ]:

allele_effects


# In[ ]:

for meta_rep in multi_meta.populations():
    assert meta_rep.numSubPop() == 7, "Correct number subpopulations before removal of the dummy population"
    meta_rep.removeSubPops(0)
    assert meta_rep.numSubPop() == 6, "Correct number after removal"


# In[ ]:

for i, meta_rep in enumerate(multi_meta.populations()):
    selection_qtd = analyze.Frq(meta_rep, triplet_qtl[i], alleles, allele_effects[i])
    selection_af = selection_qtd.allele_frequencies(meta_rep, range(meta_rep.totNumLoci()))
    selection_qtalleles = selection_qtd.rank_allele_effects(meta_rep, triplet_qtl[i], alleles, allele_effects[i])
    selection_ties = [locus for locus in range(meta_rep.totNumLoci()) 
                      if selection_af['minor', 'alleles'][locus] == selection_af['major', 'alleles'][locus]]

    for st in selection_ties:
        selection_af['major', 'alleles'][st] = list(meta_rep.dvars().alleleFreq[st])[0]
        selection_af['minor', 'alleles'][st] = list(meta_rep.dvars().alleleFreq[st])[1]
    major_minor_allele_conflicts = sum(np.equal(list(selection_af['minor', 'alleles'].values()), 
                 list(selection_af['major', 'alleles'].values())))
    
    assert major_minor_allele_conflicts == 0, "There is a tie in at least one locus."
    
    pca = analyze.PCA(meta_rep, range(meta_rep.totNumLoci()), selection_qtd)
    meta_rep_id = str(meta_rep.dvars().rep)
    
    prefix = 'rs_R' + str(meta_rep_id) + '_'
    
    minor_ac = pca.calculate_count_matrix(meta_rep, selection_af['minor', 'alleles'], 
                                      prefix + 'minor_allele_count.txt')
    
    eigendata = pca.svd(meta_rep, minor_ac)
    
    
    individual_names = {ind.ind_id: 'RS_R'+ meta_rep_id +'_G' + 
                        str(int(ind.generation)) + 
                        '_I'+str(int(ind.ind_id)) 
                        for ind in meta_rep.individuals()}
    
    ind_names_for_gwas[meta_rep_id] = individual_names
    
    meta_rep.save(prefix + 'metapopulation.pop')
    
    names_filename = prefix + 'individual_names.yaml'
    with open(names_filename, 'w') as name_stream:
        yaml.dump(individual_names, name_stream)
    
    in_dir_prefix = 'C:\\GWAS\\input\\'
    out_dir_prefix = 'C:\\GWAS\\result\\'
    config_prefix = 'C:\\GWAS\\tassel-5-standalone'
    
    
    analyze.generate_tassel_gwas_configs(in_dir_prefix, out_dir_prefix, config_prefix, 
                                         prefix, 'sim_mlm_gwas_pipeline.xml')
    
    gwas = analyze.GWAS(meta_rep, individual_names, locus_names, pos_column)
    hmap = gwas.hapmap_formatter(integer_to_snp, in_dir_prefix + prefix + 'simulated_hapmap.txt')
    phenos = gwas.trait_formatter(in_dir_prefix + prefix + 'phenotype_vector.txt')
    kinship_matrix = gwas.calc_kinship_matrix(minor_ac, selection_af, in_dir_prefix + prefix + 'kinship_matrix.txt')
    pop_struct_matrix = gwas.population_structure_formatter(eigendata, in_dir_prefix + prefix + 'structure_matrix.txt')
    pd.DataFrame(multipop.population(i).dvars().statistics).to_csv(in_dir_prefix + prefix + 'means_and_vars.txt', sep='\t')


# In[ ]:

import os


# In[ ]:

analyze.parameter_set_writer('C:\\', os.getcwd() + '\\RS_run_three_twenty_loci_', sim_params, triplet_qtl, 
                             allele_effects, genetic_structure)


# In[ ]:




# In[ ]:



