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
hapmap = pd.read_csv('clean_hapmap.txt')
genetic_map = pd.read_csv('nam_prefounders_genetic_map.txt', index_col=None,
                         sep='\t')




with open('universal_parameters.yaml', 'r') as uparms:
    u_parameters = yaml.load(uparms)


locus_names = u_parameters['locus_names']
pos_column = u_parameters['pos_column']
allele_names = u_parameters['allele_names']
snp_to_integer = u_parameters['snp_to_integer']
integer_to_snp = u_parameters['integer_to_snp']

with open('general_genetic_map_parameters.yaml', 'r') as ggmap:
    general_genetic_map_params = yaml.load(ggmap)



alleles = general_genetic_map_params['alleles']
chr_cM_positions = general_genetic_map_params['chr_cM_positions']
cM_positions = general_genetic_map_params['cM_positions']
integral_valued_loci = general_genetic_map_params['integral_valued_loci']
relative_integral_valued_loci = general_genetic_map_params['relative_integral_valued_loci']
recombination_rates = general_genetic_map_params['recombination_rates']

proto_prefix = 'C:\\Users\\DoubleDanks\\BISB\\wisser\\code\\rjwlab-scripts' \
               '\\saegus_project\\devel\\data_dump' \
               '\\fourth_generation_simulated_gwas\\'

run_prefix = 'rs_L10_H04\\'


tassel_input_dir_prefix = proto_prefix + run_prefix + 'tassel_input\\'
tassel_output_dir_prefix = proto_prefix + run_prefix + 'tassel_output\\'
tassel_config_prefix = proto_prefix + run_prefix + 'tassel_config_files\\'
various_simulation_info_prefix = proto_prefix + run_prefix + 'simulation_data\\'
populations_prefix = proto_prefix + run_prefix + 'populations\\'
parameter_prefix = proto_prefix + run_prefix + 'simulation_parameters\\'
ind_names_prefix = proto_prefix + run_prefix + 'ind_names\\'


nam = sim.loadPopulation(u_parameters['prefounder_file_name'])
sim.tagID(nam, reset=True)
nam.setSubPopName('maize_nam_prefounders', 0)

selection_statistics = {
    'aggregate': {},
    'selected': {},
    'non-selected': {}
}

ind_names_for_gwas = {i: {} for i in range(u_parameters['number_of_replicates'])}
u_parameters['meta_pop_sample_sizes'] = {i: 100 for i in range(0, u_parameters['generations_of_selection']+1, 2)}

s = simulate.Truncation(u_parameters['generations_of_selection'],
                       u_parameters['generations_of_random_mating'],
                       u_parameters['operating_population_size'],
                       u_parameters['proportion_of_individuals_saved'],
                       u_parameters['overshoot_as_proportion'],
                       u_parameters['individuals_per_breeding_subpop'],
                       u_parameters['heritability'],
                       u_parameters['meta_pop_sample_sizes'],
                       u_parameters['number_of_replicates'])


# In[10]:

founders = u_parameters['founders']
replicated_nam = sim.Simulator(nam, rep=2, stealPops=False)
pop = replicated_nam.extract(0)


# ### Run MAGIC Mating Scheme

# In[11]:

s.generate_f_one(pop, recombination_rates, u_parameters['founders'])
s.recombinatorial_convergence(pop, recombination_rates)
s.expand_by_selfing(pop, recombination_rates)
s.interim_random_mating(pop, recombination_rates)


# ## Adapting QTL and Allele Effects to Multiple Replicate Case

# In[12]:

multipop = sim.Simulator(pop, u_parameters['number_of_replicates'])
multi_meta = sim.Simulator(nam, u_parameters['number_of_replicates'], stealPops=False)


# #### Assign Each Replicate Identical Parameters
#     Determines a single random set of QTL/allele effects and assigns the
#     same information to every replicate.

triplet_qtl, allele_effects = parameters.assign_identical_qtl_parameters(multipop, alleles,
                                                                         integral_valued_loci, u_parameters['number_of_qtl'],
                                                                         u_parameters['allele_effect_parameters'])


u_parameters['meta_pop_sample_sizes'] = {i: 100 for i in range(0,
                                                               u_parameters['generations_of_selection']+1, 2)}



assert type(triplet_qtl[0]) == type([]), "Variables are flip-flopped in return."

for repid, pop_rep in enumerate(multipop.populations()):
    pop_rep.dvars().statistics = copy.deepcopy(selection_statistics)

s.replicate_selection(multipop, multi_meta, triplet_qtl, allele_effects,
                                recombination_rates)


for meta_rep in multi_meta.populations():
    assert meta_rep.numSubPop() == 7, "Correct number subpopulations before removal of the dummy population"
    meta_rep.removeSubPops(0)
    assert meta_rep.numSubPop() == 6, "Correct number after removal"


for i, meta_rep in enumerate(multi_meta.populations()):

    meta_rep_id = str(meta_rep.dvars().rep)
    replicate_prefix = 'rs_L10_H04_R' + str(meta_rep_id) + '_'

    meta_rep.dvars().triplet_qtl = triplet_qtl[i]
    meta_rep.dvars().allele_effects = allele_effects[i]
    frq = analyze.Frq(meta_rep, triplet_qtl[i], alleles, allele_effects[i])
    af = frq.allele_frequencies(meta_rep, range(meta_rep.totNumLoci()))
    #qtalleles = frq.rank_allele_effects(meta_rep, triplet_qtl[i], alleles,
    # allele_effects[i])
    ties = [locus for locus in range(meta_rep.totNumLoci())
            if af['minor', 'alleles'][locus] == af['major', 'alleles'][locus]]

    for st in ties:
        af['major', 'alleles'][st] = list(meta_rep.dvars().alleleFreq[st])[0]
        af['minor', 'alleles'][st] = list(meta_rep.dvars().alleleFreq[st])[1]
    major_minor_allele_conflicts = sum(np.equal(list(af['minor', 'alleles'].values()),
                                                list(af['major', 'alleles'].values())))

    assert major_minor_allele_conflicts == 0, "There is a tie in at least one locus."

 #   af_table = frq.allele_frq_table(meta_rep, meta_rep.numSubPop(), af,
  #
    # recombination_rates, genetic_map)
#    qtaf_table = frq.qt_allele_table(qtalleles, allele_effects[i])

    #af_table.to_csv(various_simulation_info_prefix + replicate_prefix +
    # 'allele_frequency_table.txt', sep=',', index=0)
    #qtaf_table.to_csv(various_simulation_info_prefix + replicate_prefix +
    # 'qt_allele_info.txt', sep=',', index=0)

   # del af_table, qtaf_table



    pca = analyze.PCA(meta_rep, range(meta_rep.totNumLoci()), frq)


    minor_ac = pca.calculate_count_matrix(meta_rep, af['minor', 'alleles'],
                                          various_simulation_info_prefix + replicate_prefix + 'minor_allele_count.txt')

    eigendata = pca.svd(meta_rep, minor_ac)

    individual_names = {ind.ind_id: replicate_prefix +'G' +
                        str(int(ind.generation)) +
                        '_I'+str(int(ind.ind_id))
                        for ind in meta_rep.individuals()}

    ind_names_for_gwas[meta_rep_id] = individual_names

    meta_rep.save(populations_prefix + replicate_prefix + 'metapopulation.pop')

    names_filename = ind_names_prefix + replicate_prefix + 'individual_names.yaml'
    with open(names_filename, 'w') as name_stream:
        yaml.dump(individual_names, name_stream)


    gwas = analyze.GWAS(meta_rep, individual_names, locus_names, pos_column)
    hmap = gwas.hapmap_formatter(integer_to_snp, tassel_input_dir_prefix + replicate_prefix + 'simulated_hapmap.txt')
    phenos = gwas.trait_formatter(tassel_input_dir_prefix + replicate_prefix + 'phenotype_vector.txt')
    kinship_matrix = gwas.calc_kinship_matrix(minor_ac, af, various_simulation_info_prefix + replicate_prefix + 'kinship_matrix.txt')
    pop_struct_matrix = gwas.population_structure_formatter(eigendata,
                                                            tassel_input_dir_prefix + replicate_prefix + 'structure_matrix.txt')

    analyze.generate_tassel_gwas_configs(tassel_input_dir_prefix,
                                         tassel_output_dir_prefix,
                                         tassel_config_prefix,
                                         replicate_prefix,
                                         'sim_mlm_gwas_pipeline.xml')

    #pd.DataFrame(multipop.population(i).dvars().statistics).to_csv(
    # various_simulation_info_prefix + replicate_prefix
     #                                                              +
    # 'means_and_vars.txt', sep='\t')


