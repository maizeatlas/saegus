#! usr/bin/python
# -*- coding: utf-8 -*-
import simuPOP as sim
from simuPOP import sampling
import pandas as pd
import numpy as np
import csv
import random

####*


class CalculateErrorVariance(sim.PyOperator):
    def __init__(self, heritability, *args, **kwargs):
        self.heritability = heritability
        sim.PyOperator.__init__(self, func=self.calculate_error_variance,
                                *args, **kwargs)

    def calculate_error_variance(self, pop):
        """
        Calculates the parameter ``epsilon`` to be used as the variance
        of the error distribution. The error distribution generates noise
        found in real experiments.
        """
        variance_of_g = np.var(pop.indInfo('g'))
        epsilon = variance_of_g*(1/self.heritability - 1)
        pop.dvars().epsilon = epsilon
        return True


class GenoAdditive(sim.PyOperator):
    def __init__(self, absolute_qtl, allele_effects, *args, **kwargs):
        self.absolute_qtl = absolute_qtl
        self.allele_effects = allele_effects
        sim.PyOperator.__init__(self,
                                func=self.additive_model,
                                *args, **kwargs)

    def additive_model(self, pop):
        """
        Calculates genotypic contribution ``g`` by summing the effect of each
        allele at each QTL triplet.
        """
        rep_id = pop.dvars().rep
        for ind in pop.individuals():
            genotypic_contribution = \
                sum([self.allele_effects[rep_id][locus][ind.genotype(ploidy=0)[
                    locus]] +
                     self.allele_effects[rep_id][locus][ind.genotype(ploidy=1)[
                         locus]]
                     for
                     locus
                     in self.absolute_qtl[rep_id]])
            ind.g = genotypic_contribution
        return True


class PhenotypeCalculator(sim.PyOperator):
    """
    Under a purely additive model for the time being: P = G + error.
    ``error`` is a random draw from a normal distribution with mean 0 and
    variance determined by the variance in the pre-selection population.
    The variance ``epsilon`` is calculated by another operator:
    CalculateErrorVariance.
    """
    def __init__(self, proportion_selected, *args, **kwargs):
        self.proportion_selected = proportion_selected
        sim.PyOperator.__init__(self, func=self.phenotypic_effect_calculator,
                                *args, **kwargs)

    def phenotypic_effect_calculator(self, pop):
        """
        Simulate measurement error by adding random error to genotypic
        contribution.
        """
        for ind in pop.individuals():
            ind.p = ind.g + random.normalvariate(0, pop.dvars().epsilon)
        return True


class CullPopulation(sim.PyOperator):
    def __init__(self, proportion_selected, *args, **kwargs):
        self.proportion_selected = proportion_selected
        sim.PyOperator.__init__(self, func=self.assign_fitness, *args, **kwargs)

    def assign_fitness(self, pop, subPop):
        pop.sortIndividuals('p', reverse=True)
        p = pop.indInfo('p')
        cutoff_index = -1*(int(self.proportion_selected*pop.popSize()) + 1)
        cutoff = sorted(p)[cutoff_index]
        pop.dvars().cutoff = cutoff
        pop.setIndInfo([x > cutoff for x in p], 'fitness')
        return True


class Sorter(sim.PyOperator):
    """
    Simple wrapper for Population.sortIndividuals() method.
    """
    def __init__(self, info_field, *args, **kwargs):
        self.info_field = info_field
        sim.PyOperator.__init__(self, func=self.sort_on_info_field, *args, **kwargs)

    def sort_on_info_field(self, pop):
        pop.sortIndividuals(self.info_field, reverse=True)
        return True

class MetaPopulation(sim.PyOperator):
    """
    Adds individuals to a meta population.
    """
    def __init__(self, meta_population, sample_size, *args, **kwargs):
        self.meta_population = meta_population
        self.sample_size = sample_size
        sim.PyOperator.__init__(self, func=self.add_to_meta_pop, *args,
                                **kwargs)

    def add_to_meta_pop(self, pop):
        sampled = sampling.drawRandomSample(pop,
                                            sizes=self.sample_size[
                                                pop.dvars().gen],
                                            subPops=sim.ALL_AVAIL)
        pop.dvars().ss = self.sample_size[pop.dvars().gen]
        pop.dvars().gen_sampled_from = pop.dvars().gen
        self.meta_population.addIndFrom(sampled)
        return True


class ReplicateMetaPopulation(sim.PyOperator):
    """
    Operator to sample individuals from the current generations to be added to a pre-initialized meta-population.
    In the replicate case a population is mapped onto its corresponding replicate population by dictionary variable
    'rep' which is assigned by simuPOP automatically during instantiation of a Simulator object.
    param: Both replicates are assigned as attributes to a wgs.EnhancedPopulation object.
    """
    def __init__(self, meta_replicates, sample_size, *args, **kwargs):
        self.meta_replicates = meta_replicates
        self.sample_size = sample_size
        sim.PyOperator.__init__(self, func=self.add_to_meta_pop, *args, **kwargs)

    def add_to_meta_pop(self, pop):
        rep_id = pop.dvars().rep
        sampled = sampling.drawRandomSample(pop, sizes=self.sample_size[pop.dvars().gen])
        pop.dvars().ss = self.sample_size[pop.dvars().gen]
        pop.dvars().gen_sampled_from = pop.dvars().gen
        self.meta_replicates.population(rep_id).addIndFrom(sampled)
        return True


class SaveMetaPopulations(sim.PyOperator):
    """
    Operator to sample individuals from the current generations to be added to a pre-initialized meta-population.
    In the replicate case a population is mapped onto its corresponding replicate population by dictionary variable
    'rep' which is assigned by simuPOP automatically during instantiation of a Simulator object.
    param: Both replicates are assigned as attributes to a wgs.EnhancedPopulation object.
    """
    def __init__(self, meta_replicates, *args, **kwargs):
        self.meta_replicates = meta_replicates
        sim.PyOperator.__init__(self, func=self.save_meta_pop_in_native_format, *args, **kwargs)

    def save_meta_pop_in_native_format(self, pop):
        rep_id = pop.dvars().rep
        meta_filename = 'replicate_' + str(rep_id) + '_meta_pop.pop'
        self.meta_replicates.population(rep_id).save(meta_filename)
        return True


class InfoAndGenotypeWriter(sim.PyOperator):
    """
    Operator to output values of individual infoFields and genotype matrix to file. Very similar
    to simuPOP.utils.Exporter; however, allows for greater developmental flexibility.
    """
    def __init__(self, output_file_name: str, *args, **kwargs):
        """
        output_file_name should not have a file extension.
        """
        self.output_file_name = output_file_name
        sim.PyOperator.__init__(self, func=self.info_and_genotype_writer, *args, **kwargs)

    def info_and_genotype_writer(self, pop):
        full_file_name = self.output_file_name + "_" + str(pop.dvars().gen) + ".txt"
        header = ['ind_id', 'mother_id', 'father_id', 'g', 'p']
        genotype_header = list(range(2*pop.totNumLoci()))
        header.extend(genotype_header)
        with open(full_file_name, 'w') as pop_info:
            info_writer = csv.writer(pop_info, delimiter=',')
            info_writer.writerow(header)
            for ind in pop.individuals():
                info_writer.writerow([ind.ind_id, ind.mother_id,
                                      ind.father_id, ind.g, ind.p,
                                      ind.genotype()])
        return True


class RandomlyAssignFemaleFitness(sim.PyOperator):
    """
    Chooses ``size_breeding_subpopulation`` individuals to be eligible for
    mating. 0 is the default value for ``female_fitness``. Individuals who
    have ``female_fitness`` = 0 cannot be picked as mates.
    Individuals who have female fitness = 1 can be chosen as a 'female'.
    """
    def __init__(self, size_breeding_subpopulation, *args, **kwargs):
        self.size_breeding_subpopulation = size_breeding_subpopulation
        sim.PyOperator.__init__(self, func=self.choose_breeding_individuals_randomly, *args, **kwargs)

    def choose_breeding_individuals_randomly(self, pop):
        random_individual_ids = random.sample(list(pop.indInfo('ind_id')), self.size_breeding_subpopulation)
        for id in random_individual_ids:
            pop.indByID(id).female_fitness = 1
        return True


class RandomlyAssignMaleFitness(sim.PyOperator):
    """
    Operator which parallels the purpose of RandomlyAssignFemaleFitness. Because we are working in a plant speices there
    is the possibility of selfing. Thus breeding females are also counted as breeding males. In this case we supplement
    the number of breeding females with an additional number of individuals who will only be used as males.
    """
    def __init__(self, additional_breeding_males, *args, **kwargs):
        self.additonal_breeding_males = additional_breeding_males
        sim.PyOperator.__init__(self, func=self.randomly_choose_breeding_males, *args, **kwargs)

    def randomly_choose_breeding_males(self, pop):
        breeding_female_ids = [ind.ind_id for ind in pop.individuals() if ind.female_fitness == 1.0]
        # Now every breeding female is also considered as a breeding male.
        for female_id in breeding_female_ids:
            pop.indByID(female_id).male_fitness = 1.0
        non_breeding_ids = [ind.ind_id for ind in pop.individuals() if ind.ind_id not in breeding_female_ids]
        additional_breeding_male_ids = random.sample(non_breeding_ids, self.additonal_breeding_males)
        for male_id in additional_breeding_male_ids:
            pop.indByID(male_id).male_fitness = 1.0
        return True


class DiscardRandomOffspring(sim.PyOperator):
    """
    Operator to choose ``number_to_remove`` individuals at random to remove
    from the offspring population. Simulates the effect of randomly picking
    seed to plant from the larger population.
    """
    def __init__(self, number_to_remove, *args, **kwargs):
        self.number_to_remove = number_to_remove
        sim.PyOperator.__init__(self, func=self.remove_indices, *args, **kwargs)

    def remove_indices(self, pop):
        inds = list(range(pop.popSize()))
        removed_inds = random.sample(inds, self.number_to_remove)
        pop.removeIndividuals(removed_inds)
        return True

def calculate_error_variance(pop, heritability):
    """
    Calculates the parameter ``epsilon`` to be used as the variance
    of the error distribution. The error distribution generates noise
    found in real experiments.
    """
    variance_of_g = np.var(pop.indInfo('g'))
    epsilon = variance_of_g*(1/heritability - 1)
    pop.dvars().epsilon = epsilon

def phenotypic_effect_calculator(pop):
    """
    Simulate measurement error by adding random error to genotypic
    contribution.
    """
    for ind in pop.individuals():
        ind.p = ind.g + random.normalvariate(0, pop.dvars().epsilon)

def assign_additive_g(pop, qtl, allele_effects):
    """
    Calculates genotypic contribution ``g`` by summing the effect of each
    allele at each quantitative trait locus. Assumes a bi-allelic case.
    """
    for ind in pop.individuals():
        genotypic_contribution = \
            sum([
                    allele_effects[locus][ind.genotype(ploidy=0)[locus]] +\
                    allele_effects[locus][ind.genotype(ploidy=1)[locus]]
                 for locus
                 in qtl])
        ind.g = genotypic_contribution