#!/usr/bin/python
"""
Typical driver script for a run of artificial selection.
"""
import simuOpt
simuOpt.setOptions(alleleType='short', numThreads=4, optimized=False, quiet=True)
import simuPOP as sim
from wgs import breed, operators, selection, helpers, parser, parameterizer

generations_of_selection = 1
generations_of_random_mating = 1
selection_population_size = 2000
proportion_of_individuals_saved = 0.05
overshoot_as_proportion = 0.50
individuals_per_breeding_subpop = 5
heritability = 0.7
metapop_sample_size = 100

nam = sim.loadPopulation('nam_prefounders.pop')
sim.tagID(prefounders, reset=True)
founders = [1, 5, 7, 8]
pop = prefounders.clone()
ep.arbitrary_ordering_of_founders(prefounders, pop, founders)
ep.generate_f_one(pop)
ep.generate_f_two(pop)
ep.mate_and_merge(pop)
ep.interim_random_mating(pop, generations_of_random_mating)
ep.replicates = sim.Simulator(pop, stealPops=False, rep=number_of_replicates)
ep.meta_replicates = sim.Simulator(pop, rep=number_of_replicates)
"""
# Remove individuals which are initially present in Simulator objects
"""
for meta_rep in ep.meta_replicates.populations():
    meta_rep.removeIndividuals(list(range(meta_rep.popSize())))


ep.pure_selection(ep.replicates, ep.meta_replicates)


syn = helpers.Synbreed()
"""
# Allele frequencies and output
"""
# Note: In PCA I should get the same results from the minor allele frequency
# as the major allele frequency
"""
***
TODO: Write a function to output the genotype matrix in a TASSEL format.
Note: Instead of writing two separate files for genotypes and MAC just
use a conversion vector.
TODO: Output the mean and variances for every generation
TODO: Output the super table
***
"""
for meta_rep in ep.meta_replicates.populations():
    splitlets = helpers.lineage_triplets_to_splitlets(meta_rep)
    splitlet_frequencies = helpers.splitlet_frequencies(meta_rep, splitlets)
    minor_alleles = helpers.min_finder(meta_rep, splitlet_frequencies)
    major_alleles = helpers.max_finder(meta_rep, splitlet_frequencies)
    # Output operations
    minor_genotype_matrix_filename = 'rep_' + str(meta_rep.dvars().rep) + 'minor_genotype_matrix.txt'
    phenotype_filename = 'rep_' + str(meta_rep.dvars().rep) + 'phenotypes.txt'
    pedigree_filename = 'rep_' + str(meta_rep.dvars().rep) + 'pedigree.txt'
     
    syn.write_population_data(meta_rep, splitlets, minor_alleles,
                                     minor_genotype_matrix_filename,
                                     phenotype_filename,
                                     pedigree_filename)



