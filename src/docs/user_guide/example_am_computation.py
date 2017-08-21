#!/home/vakanas/anaconda43/python3.6
import simuOpt
simuOpt.setOptions(alleleType='short', quiet=True, numThreads=4)
import simuPOP as sim
import numpy as np
import pandas as pd
import random
import h5py
from saegus import analyze, operators, parameters
#np.set_printoptions(suppress=True, precision=5)


example_pop = sim.loadPopulation('example_pop.pop')
example_pop.addInfoFields(['ind_id', 'mother_id', 'father_id', 'g', 'p'])
sim.tagID(example_pop)
sim.stat(example_pop, numOfSegSites=sim.ALL_AVAIL,
         vars=['numOfSegSites','segSites', 'fixedSites'])
sim.stat(example_pop, alleleFreq=sim.ALL_AVAIL)

segregating_loci = example_pop.dvars().segSites
allele_states = analyze.gather_allele_data(example_pop)
allele_frequencies = analyze.gather_allele_frequencies(example_pop, allele_states)
gwas = analyze.GWAS(example_pop, np.array(segregating_loci, dtype=np.int_), allele_states[:, 3], 'example')
count_matrix = gwas.calculate_count_matrix('example_count_matrix.txt')
gwas.hapmap_formatter(hapmap_file_name = 'example_hapmap.txt')
eigenvalues, eigenvectors = gwas.pop_struct_eigendecomp(count_matrix)
gwas.population_structure_formatter(eigenvalues, eigenvectors,
                         pop_struct_file_name='example_structure.txt', number_of_pcs=2)
gwas.calc_kinship_matrix(count_matrix,
                     kinship_matrix_file_name='example_kinship.txt')

qtl = sorted(random.sample(segregating_loci, 20))
trait = parameters.Trait()
allele_effects_table = trait.construct_allele_effects_table(example_pop, qtl, random.expovariate, 1)
allele_effects_array = trait.construct_ae_array(allele_effects_table, qtl)
heritability = 0.7
operators.calculate_g(example_pop, allele_effects_array)
operators.calculate_error_variance(example_pop, heritability)
operators.calculate_p(example_pop)
tassel_trait = gwas.trait_formatter(trait_file_name='example_trait.txt')

example_trait_data = h5py.File('8517example_trait_data.hdf5')
example_trait_data['allele/states'] = allele_states
example_trait_data['allele'].attrs['columns'] = \
list(map(np.string_, ['locus', 'alpha', 'omega', 'minor', 'major']))
example_trait_data['allele/effects'] = allele_effects_table
example_trait_data ['allele/frequencies'] = allele_frequencies
example_trait_data['qtl'] = np.array(qtl)
example_trait_data['segregating_loci'] = segregating_loci


example_trait_data.close()

gwas.tassel_gwas_config(config_template='example_config_file.xml',
                       hapmap_file_name='example_hapmap.txt',
                       kinship_file_name='example_kinship.txt',
                       trait_file_name='example_trait.txt',
                       structure_file_name='example_structure.txt',
                       output_prefix='example_')