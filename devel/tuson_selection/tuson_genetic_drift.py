#!/usr/bin/python3.4
import sys
import csv
import json
import numpy as np
import pandas as pd
import scipy.stats as st
import collections as col
import random
import simuOpt

simuOpt.setOptions(alleleType='short', numThreads=10, optimized=True,
                   quiet=True)
from wgs import operators, helpers, parameterizer, breed, parser
import simuPOP as sim
from simuPOP import sampling
from simuPOP import utils


class TusonDrift(object):
    """
    Encapsulates all the methods and parameters to use during the simulation
    genetic drift.
    """

    def __init__(self, population_size=1, number_of_replicates=1,
                 number_of_breeding_males=1, number_of_breeding_females=1,
                 number_of_generations=1, genetic_map_filename='',
                 population_structure_filename='',
                 founder_population_filename=''):
        self.population_size = population_size
        self.number_of_replicates = number_of_replicates
        self.number_of_breeding_males = number_of_breeding_males
        self.number_of_breeding_females = number_of_breeding_females
        self.number_of_generations = number_of_generations
        self.genetic_map_filename = genetic_map_filename
        self.population_structure_filename = population_structure_filename
        self.founder_population_filename = founder_population_filename

    def initialize_tuson_pop(self, chromosome_lengths, info_fields,
                             all_genotypes_handle):
        tuson_founders = sim.Population(size=105, loci=chromosome_lengths,
                                        infoFields=info_fields)
        for ind, genotype in zip(tuson_founders.individuals(),
                                 all_genotypes_handle):
            ind.setGenotype(genotype)
        sim.tagID(tuson_founders, reset=True)
        return tuson_founders

    def tuson_drift_simulation(self, pop, meta_population, meta_sample_size,
                               recombination_rates):
        female_chooser = sim.RandomParentChooser(
            selectionField='female_fitness')
        male_chooser = sim.RandomParentChooser(selectionField='male_fitness')
        print("Beginning simulation of genetic drift with parameters:")
        breeding_params = {
            'n_breeding_females': self.number_of_breeding_females,
            'n_breeding_males': self.number_of_breeding_males,
            'sample_sizes': meta_sample_size}
        print(
            "Breeding females: {n_breeding_females}, breeding males: {n_breeding_males}, sample sizes: {sample_sizes}".format(
                **breeding_params))
        return pop.evolve(
            initOps=[
                operators.MetaPopulation(meta_population, meta_sample_size),
                sim.PyEval(
                    r'"Initial: Sampled %d individuals from generation %d\n" % (ss, gen_sampled_from)'),
                ],
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
                sim.InfoExec('generation=gen'),
                operators.MetaPopulation(meta_population, meta_sample_size,
                                         at=[2, 4, 6, 8]),
                # Evaluation specifying the generations should be the same as the evaluation at every generation.
                sim.PyEval(
                    r'"Sampled %d individuals from generation %d.\n" % (ss, gen_sampled_from)',
                    at=[2, 4, 6, 8]),
                operators.RandomlyAssignFemaleFitness(
                    self.number_of_breeding_females),
                operators.RandomlyAssignMaleFitness(
                    self.number_of_breeding_males),
            ],
            matingScheme=sim.HomoMating(
                sim.CombinedParentsChooser(female_chooser, male_chooser,
                                           allowSelfing=True),
                sim.OffspringGenerator(ops=[
                    sim.ParentsTagger(),
                    sim.IdTagger(),
                    sim.PedigreeTagger(),
                    sim.Recombinator(rates=recombination_rates)],
                ),
            ),
            finalOps=[
                sim.InfoExec('generation=gen'),
                operators.MetaPopulation(meta_population, meta_sample_size),
                sim.PyEval(
                    r'"Final: Sampled %d individuals from generation %d\n" % (ss, gen_sampled_from)')
            ],
            gen=self.number_of_generations)

    def replicate_tuson_drift_simulation(self, pop, meta_population,
                                         meta_sample_size,
                                         recombination_rates):
        for replicate in pop.populations():
            replicate.dvars().gen = 0
        female_chooser = sim.RandomParentChooser(
            selectionField='female_fitness')
        male_chooser = sim.RandomParentChooser(selectionField='male_fitness')
        print("Beginning simulation of genetic drift with parameters:")
        breeding_params = {
            'n_breeding_females': self.number_of_breeding_females,
            'n_breeding_males': self.number_of_breeding_males,
            'sample_sizes': meta_sample_size}
        print("Breeding females: {n_breeding_females}, breeding males: "
              "{n_breeding_males}, "
              "sample sizes: {sample_sizes}".format(**breeding_params))
        return pop.evolve(
            initOps=[
                operators.ReplicateMetaPopulation(meta_population,
                                                  meta_sample_size),
                sim.PyEval(
                    r'"Initial: Sampled %d individuals from generation %d Replicate: %d.\n" % (ss, gen_sampled_from, rep)'),
            ],
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
                sim.InfoExec('generation=gen'),
                operators.ReplicateMetaPopulation(meta_population,
                                                  meta_sample_size,
                                                  at=[2, 4, 6, 8]),
                # Evaluation specifying the generations should be the same as the evaluation at every generation.
                sim.PyEval(
                    r'"Sampled %d individuals from generation %d from replicate: %d.\n" % (ss, gen_sampled_from, rep)',
                    at=[2, 4, 6, 8]),
                operators.RandomlyAssignFemaleFitness(
                    self.number_of_breeding_females),
                operators.RandomlyAssignMaleFitness(
                    self.number_of_breeding_males),
            ],
            matingScheme=sim.HomoMating(
                sim.CombinedParentsChooser(female_chooser, male_chooser,
                                           allowSelfing=True),
                sim.OffspringGenerator(ops=[
                    sim.ParentsTagger(),
                    sim.IdTagger(),
                    sim.PedigreeTagger(),
                    sim.Recombinator(rates=recombination_rates)],
                ),
            ),
            finalOps=[
                sim.InfoExec('generation=gen'),
                operators.ReplicateMetaPopulation(meta_population,
                                                  meta_sample_size),
                sim.PyEval(
                    r'"Final: Sampled %d individuals from generation %d\n" % (ss, gen_sampled_from)')
            ],
            gen=self.number_of_generations)

    def population_structure_guided_expansion(self, pop, recombination_rates):
        """
        Uses a population structure matrix to determine the probability of
        selecting a second parent given the first parent's probability mass
        function.
        """
        ps_pc = breed.ForcedPopulationStructureParentChooser(
            self.population_size)
        print(
            "Executing population expansion using estimated population structure.")
        return pop.evolve(
            initOps=sim.InitSex(),
            preOps=[
                sim.InfoExec('generation=gen'),
            ],
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(ps_pc.forced_structure_parent_chooser),
                sim.OffspringGenerator(ops=[
                    sim.IdTagger(),
                    sim.PedigreeTagger(),
                    sim.Recombinator(rates=recombination_rates),
                ]),
                subPopSize=self.population_size),
            gen=1)

    def expansion_through_random_mating(self, pop, expanded_pop_size,
                                        recombination_rates):

        # The purpose of this function is to use the simuPOP pre-defined mating scheme
        # RandomMating to grow the population to an arbitrary size.
        # Self-pollination occurs frequently in maize so we need use HermaphroditicMating
        # instead of RandomMating.
        return pop.evolve(
            initOps=sim.InitSex(),
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
                sim.InfoExec('generation=gen'),
            ],
            matingScheme=sim.HermaphroditicMating(
                ops=[sim.Recombinator(rates=recombination_rates),
                     sim.IdTagger(),
                     sim.PedigreeTagger()],
                subPopSize=expanded_pop_size),
            gen=1)

    def find_fixed_sites(self, founder_population, threshold, error):
        cloned_pop = founder_population.clone()
        clone_ps = parameterizer.PopulationStructure(cloned_pop,
                                                     self.population_structure_filename,
                                                     threshold, error)
        # st.setup_mating_structure()
        valid_inds = list(cloned_pop.indInfo('ind_id'))
        sim.stat(cloned_pop, numOfSegSites=True,
                 vars=['numOfFixedSites', 'fixedSites'])
        num_fixed = cloned_pop.dvars().numOfFixedSites
        fixed_sites = list(cloned_pop.dvars().fixedSites)
        return fixed_sites, valid_inds

    def insert_allele_frequencies_into_aggregate_matrix(self, number_of_reps,
                                                        meta_population,
                                                        minor_allele):
        """
        number_of_reps: integer specifying number of replicates to include
        inside of one allele frequency matrix.

        meta_population: Multiple replicate meta population.
        minor_allele_list: List or numpy.array of minor alleles.
        :param number_of_reps: Replicates of the population
        :type number_of_reps: int
        :param meta_population: Multi-replicate population container
        :type meta_population: sim.Simulator
        :param minor_allele: Locus: Minor Allele key-value pairs
        :type minor_allele: dict
        :return: Minor allele frequencies of multiple replicates
        :rtype: np.array
        """
        number_of_rows = 6 * number_of_reps
        aggregate_frequency_matrix = np.zeros((number_of_rows, 44447))
        row_indices = list(range(number_of_rows))
        print(
            'Calculating allele frequencies for {number_reps} replicates and writing them to an aggregate matrix.'.format(
                number_reps=number_of_reps))
        for replicate in meta_population.populations():
            print("Replicate: {rep_id}".format(rep_id=replicate.dvars().rep))
            sim.stat(replicate, alleleFreq=sim.ALL_AVAIL,
                     vars=['alleleFreq_sp'])
            subpops_and_gens = [(1, 0), (2, 2), (3, 4), (4, 6), (5, 8),
                                (6, 10)]
            for sp, gen in subpops_and_gens:
                row_index = row_indices.pop(0)
                print("Row index: {row_index}".format(row_index=row_index))
                rep = replicate.dvars().rep
                aggregate_frequency_matrix[row_index, 0] = gen
                aggregate_frequency_matrix[row_index, 1] = rep
                aggregate_frequency_matrix[row_index, 2:] = [
                    replicate.dvars(sp).alleleFreq[locus][ma] for ma, locus in
                    zip(minor_allele.values(), range(44445))]
        return aggregate_frequency_matrix

    def insert_minor_genotype_frequencies_into_aggregate_matrix(self,
                                                                minor_homozygotes,
                                                                minor_heterozygotes,
                                                                number_of_reps=1,
                                                                meta_population=None):
        """
        number_of_reps: integer specifying number of replicates to include inside of one allele frequency matrix.
        meta_population: Multiple replicate meta population.
        minor_allele_list: List or numpy.array of minor alleles.
        """

        number_of_rows = 6 * number_of_reps
        homozygote_frequency_matrix = np.zeros((number_of_rows, 44445 + 2))
        heterozygote_frequency_matrix = np.zeros((number_of_rows, 44445 + 2))
        row_indices = list(range(number_of_rows))
        print(
            "Calculating minor-allele genotype frequencies for {number_reps} replicates and writing them to an aggregate matrix.".format(
                number_reps=number_of_reps))
        for replicate in meta_population.populations():
            print("Replicate: {rep_id}".format(rep_id=replicate.dvars().rep))
            sim.stat(replicate, genoFreq=sim.ALL_AVAIL, vars=['genoFreq_sp'])
            subpops_and_gens = [(1, 0), (2, 2), (3, 4), (4, 6), (5, 8),
                                (6, 10)]
            for sp, gen in subpops_and_gens:
                row_index = row_indices.pop(0)
                print("Row index: {row_index}".format(row_index=row_index))
                rep = replicate.dvars().rep
                # Homozygote
                homozygote_frequency_matrix[row_index, 0] = gen
                homozygote_frequency_matrix[row_index, 1] = rep
                homozygote_frequency_matrix[row_index, 2:] = [
                    replicate.dvars(sp).genoFreq[locus][
                        minor_homozygotes[locus]]
                    for locus in range(44445)]
                # Heterozygote
                heterozygote_frequency_matrix[row_index, 0] = gen
                heterozygote_frequency_matrix[row_index, 1] = rep
                heterozygote_frequency_matrix[row_index, 2:] = [(
                                                                replicate.dvars(
                                                                    sp).genoFreq[
                                                                    locus][
                                                                    minor_heterozygotes[
                                                                        locus][
                                                                        0]] +
                                                                replicate.dvars(
                                                                    sp).genoFreq[
                                                                    locus][
                                                                    minor_heterozygotes[
                                                                        locus][
                                                                        1]])
                                                                for locus in
                                                                range(44445)]
        return homozygote_frequency_matrix, heterozygote_frequency_matrix

    def randomly_convert_fixed_sites(self, pop, fixed_sites):
        """
        Randomly converts fixed sites in pop to nonfixed_sites by changing the allele state at that site.
        """
        alleles = [0, 1, 2, 3]
        random.shuffle(fixed_sites)
        for site in fixed_sites:
            random_id = random.choice(pop.indInfo('ind_id'))
            random_individual = pop.indByID(random_id)
            current_allele_state = random_individual.allele(site)
            possible_replacements = [allele for allele in alleles if
                                     allele != current_allele_state]
            replacement_allele = random.choice(possible_replacements)
            random_individual.setAllele(replacement_allele, site)

    def read_tuson_genetic_map(self):
        self.genetic_map = pd.read_csv(self.genetic_map_filename, sep='\t',
                                       usecols=[1, 3], skip_header=True)

    def initialize_meta_population(self, pop, number_of_reps=1):
        if number_of_reps == 1:
            sim_meta = sim.Simulator(pop, rep=1, stealPops=False)
            meta_pop = sim_meta.extract(0)
            meta_pop.removeSubPops(0)
            return meta_pop
        else:
            proto_meta = pop.clone()
            proto_meta.setSubPopName('meta', 0)
            proto_meta.removeSubPops(0)
            sim_meta = sim.Simulator(proto_meta, rep=number_of_reps,
                                     stealPops=False)
        return sim_meta


server_prefix = '/home/jjdoc/tuson/parameter_files/'
map_filename = server_prefix + 'genetic_map.txt'
popst_filename = server_prefix + 'population_structure_matrix.xlsx'
founder_filename = server_prefix + 'tuson.pop'
popsize = 10000
females = 400
males = 400
generations = 10
replicates = 100

td = TusonDrift(genetic_map_filename=map_filename, population_size=popsize,
                number_of_breeding_females=females,
                number_of_breeding_males=males,
                number_of_generations=generations,
                founder_population_filename=founder_filename,
                population_structure_filename=popst_filename)

recom_rates = parser.parse_recombination_rates(map_filename)

tuson = sim.loadPopulation(td.founder_population_filename)
tuson.setSubPopName('tuson', 0)

sim.tagID(tuson, reset=True)
# tuson_meta = td.initialize_meta_population(tuson, number_of_reps=replicates)


# sites, inds = td.find_fixed_sites(tuson, 0.15, 0.03)
sim.stat(tuson, numOfSegSites=sim.ALL_AVAIL,
         vars=['numOfFixedSites', 'fixedSites'])

sites = tuson.dvars().fixedSites




# Assigns popoulation structure
ps = parameterizer.PopulationStructure(tuson, td.population_structure_filename,
                                       0.15, 0.03)
structure = ps.generate_population_structure()
tuson.dvars().mating_pmfs = {i: st.rv_discrete(
    values=([0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            structure[i]))
                             for i in range(1, tuson.popSize() + 1)}
tuson.dvars().assigned_structure = ps.assign_population_structure(structure)

td.randomly_convert_fixed_sites(tuson, sites)
sim.stat(tuson, alleleFreq=sim.ALL_AVAIL)
sim.stat(tuson, genoFreq=sim.ALL_AVAIL)
minor_alleles = helpers.find_minor_alleles(tuson)
major_alleles = helpers.find_major_alleles(tuson)
minor_homozygotes_by_locus = {i: (minor_alleles[i], minor_alleles[i]) for i in
                              range(tuson.totNumLoci())}
minor_heterozygotes_by_locus = {
i: [(minor_alleles[i], major_alleles[i]), (major_alleles[i], minor_alleles[i])]
for i in range(tuson.totNumLoci())}

primary_subpopulations = [tuson.dvars().assigned_structure[int(indid)][0] for
                          indid in list(tuson.indInfo('ind_id'))]
tuson.setIndInfo(values=primary_subpopulations, field='primary')
primary_subpopulation_splitter = sim.InfoSplitter(field='primary',
                                                  values=[0.0, 1.0, 2.0, 3.0,
                                                          4.0, 5.0])
tuson.setVirtualSplitter(primary_subpopulation_splitter)

experimental_meta_population_sample_size = {
    0: 105,
    2: 56,
    4: 55,
    6: 54,
    8: 56,
    10: 55
}

uniform_meta_population_sample_size = {
    0: 400,
    2: 400,
    4: 400,
    6: 400,
    8: 400,
    10: 400,
}

tuson.dvars().sample_sizes = experimental_meta_population_sample_size

td.population_structure_guided_expansion(tuson, recom_rates)

sim_params = {
    'genetic_map': map_filename,
    'population_structure': popst_filename,
    'founder_population_filename': founder_filename,
    'simulated_pop_size': popsize,
    'number_breeding_males': males,
    'number_breeding_females': females,
    'number_generations': generations,
    'number_of_reps': replicates,
    'sample_sizes': experimental_meta_population_sample_size,
}

with open('simulation_parameters.json', 'w') as simparams:
    json.dump(sim_params, simparams, indent=5)

for batch in range(92, 101):
    print("Executing batch: {batch_id}.".format(batch_id=batch))
    subpops_and_gens = [(1, 0), (2, 2), (3, 4), (4, 6), (5, 8), (6, 10)]
    tuson_replicates = sim.Simulator(tuson, rep=replicates, stealPops=False)
    tuson_meta_replicates = td.initialize_meta_population(tuson,
                                                          number_of_reps=replicates)
    td.replicate_tuson_drift_simulation(tuson_replicates,
                                        tuson_meta_replicates,
                                        experimental_meta_population_sample_size,
                                        recom_rates)
    print("Multi-replicate drift simulation complete.")
    allelefrequencies = td.insert_allele_frequencies_into_aggregate_matrix(
        replicates, tuson_meta_replicates,
        minor_alleles)
    np.savetxt('batch_' + str(batch) + '_af.txt', allelefrequencies,
               fmt='%.3f', delimiter='\t')
    homofrequencies, heterofrequencies = td.insert_minor_genotype_frequencies_into_aggregate_matrix(
        minor_homozygotes_by_locus,
        minor_heterozygotes_by_locus,
        number_of_reps=replicates,
        meta_population=tuson_meta_replicates)
    np.savetxt('batch_' + str(batch) + '_homofrequencies.txt', homofrequencies,
               fmt='%.3f', delimiter='\t')
    np.savetxt('batch_' + str(batch) + '_heterofrequencies.txt',
               heterofrequencies, fmt='%.3f', delimiter='\t')
# if batch % 3 == 0:
#        del tuson_replicates, tuson_meta_replicates, frqmatrix
