import simuOpt
simuOpt.setOptions(alleleType='short', optimized=True, numThreads=4, quiet=True)
import simuPOP as sim
import pandas as pd
from saegus import breed, operators, simulate, analyze, parse, parameters
import shelve
import numpy as np
import random
import xml.etree.ElementTree as ET
import lxml.etree as etree
#import h5py
#import collections as col
np.set_printoptions(suppress=True, precision=3)


input_file_prefix = '/home/vakanas/BISB/rjwlab-scripts/saegus_project/devel/magic/1478/'


mg = analyze.MultiGeneration('epsilon')
run_id = 'epsilon'
generations = 10
heritability = 0.7
number_of_qtl = 50
number_of_replicates = 2
founders = [[2, 26], [3, 25], [4, 24], [5, 23]]
os_per_pair = 500
recombination_rates = [0.01]*1478
prefounders = sim.loadPopulation(input_file_prefix+'bia_prefounders.pop')
config_file_template = input_file_prefix + 'gwas_pipeline.xml'



sim.tagID(prefounders, reset=True)
alleles = np.array(pd.read_hdf(input_file_prefix+'parameters/alleles_at_1478_loci.hdf'))


rdm_populations = sim.Simulator(prefounders, number_of_replicates, stealPops=False)
rdm_magic = breed.MAGIC(rdm_populations, founders, recombination_rates)
sim.tagID(prefounders, reset=27)
rdm_magic.generate_f_one(founders, os_per_pair)

sim.stat(rdm_populations.population(0), alleleFreq=sim.ALL_AVAIL)
af = analyze.allele_data(rdm_populations.population(0), alleles, list(range(len(alleles))))
minor_alleles = np.asarray(af.minor_allele, dtype=np.int8)
rdm_magic.recombinatorial_convergence(rdm_populations, len(founders), os_per_pair)


study = analyze.Study(run_id)

sim.stat(rdm_populations.population(0), numOfSegSites=sim.ALL_AVAIL, vars=['segSites'])
concordant_segregating_loci = rdm_populations.population(0).dvars().segSites


qtl = sorted(random.sample(concordant_segregating_loci, number_of_qtl))



#sets_of_segregating_loci = study.seg_loci_among_samples(rdm_populations.population(0))
#concordant_segregating_loci = list(sets_of_segregating_loci.keys())[0]
#qtl = sorted(random.sample(concordant_segregating_loci, number_of_qtl))


additive_trait = parameters.Trait()
allele_effects = additive_trait.assign_allele_effects(alleles, qtl,
                                                      random.expovariate, 1,
                                                     multiplicity=1)

ae_array = additive_trait.convert_allele_effects_into_array(
    prefounders.totNumLoci(), 6, allele_effects)

sampling_generations = [i for i in range(2, 10, 2)]
sample_sizes = {i: 100 for i in range(generations + 1)}

rdm_meta_populations = {rep: [] for rep in range(number_of_replicates)}
rdm_mating = simulate.RandomMating(generations, 2000, 0.05, 0.5, 5, heritability, sample_sizes)
rdm_mating.replicate_random_mating(rdm_populations,
                                   rdm_meta_populations, qtl,
                                   ae_array, recombination_rates)


repz = rdm_meta_populations[0]
analyze.combine_population_samples(repz)
meta_pop = repz[0]


sim.stat(meta_pop, alleleFreq=sim.ALL_AVAIL)
indir = "/home/vakanas/tassel-5-standalone/input/"
outdir = "/home/vakanas/tassel-5-standalone/output/"
rep_id_name = "0"

gwas = analyze.GWAS(meta_pop, list(range(meta_pop.totNumLoci())), run_id)

ccm = gwas.calculate_count_matrix(minor_alleles, list(range(meta_pop.totNumLoci())))
ps_svd = gwas.pop_struct_svd(ccm)
name = run_id+'_'+rep_id_name
gwas.population_structure_formatter(ps_svd, indir+name+'_structure_matrix.txt')


int_to_snp_map = {0:'A', 1:'C', 2:'G', 3:'T', 4:'-', 5:'+'}
locus_names = list(concordant_segregating_loci)
alleles_column = ['NA']*len(concordant_segregating_loci)
chromosomes = [meta_pop.chromLocusPair(locus)[0]+1 for locus in concordant_segregating_loci]
gwas.hapmap_formatter(concordant_segregating_loci, alleles_column,
                      locus_names, chromosomes,
                      locus_names,
                      indir+name+'_simulated_hapmap.txt')

minor_allele_frequency_table = analyze.minor_allele_frequencies_table(
        meta_pop.dvars().alleleFreq, minor_alleles)

minor_allele_frequency_table.to_csv(outdir+name+'_maf_table.txt', sep='\t')

quantitative_allele_table = analyze.generate_allele_effects_table(meta_pop.dvars().alleleFreq, alleles, ae_array)
quantitative_allele_table.to_csv(outdir+name+'_quant_allele_table.txt', sep='\t')


minor_allele_frequencies = np.array(
    minor_allele_frequency_table.minor_frequency)

gwas.calc_kinship_matrix(ccm,
                         minor_allele_frequencies,
                         indir+name+'_kinship_matrix.txt')

gwas.trait_formatter(indir+name+'_trait_vector.txt')


tree = ET.parse(config_file_template)
root = tree.getroot()
lxml_tree = etree.fromstring(ET.tostring(root))
lxml_root = lxml_tree.getroottree()
lxml_root.find('fork1/h').text = indir+name+'_simulated_hapmap.txt'
lxml_root.find('fork2/t').text = indir+name+'_trait_vector.txt'
lxml_root.find('fork3/q').text = indir+name+'_structure_matrix.txt'
lxml_root.find('fork4/k').text = indir+name+'_kinship_matrix.txt'

lxml_root.find('combine6/export').text = outdir+name+'_out_'
lxml_root.write("/home/vakanas/tassel-5-standalone/"+"R"+rep_id_name+'_'+
                run_id+'_'+"_sim_gwas_pipeline.xml",
                encoding="UTF-8",
                method="xml",
                xml_declaration=True,
                standalone='',
                pretty_print=True)

# End Phase One