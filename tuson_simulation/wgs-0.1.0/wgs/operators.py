#! usr/bin/python
# -*- coding: utf-8 -*-
import simuPOP as sim
from simuPOP import sampling
import pandas as pd
import numpy as np
import itertools as ite
import csv
import random
import collections as col
import statistics as stat
from . import helpers

####*


class MinorAlleleFrequencyWriter(sim.PyOperator):

    def __init__(self, rows_per_rep, minor_allele_reference, *args, **kwargs):
        self.rows_per_rep = rows_per_rep
        self.minor_allele_reference = minor_allele_reference
        sim.PyOperator.__init__(self, func=self.af_per_rep_writer, *args, **kwargs)

    def af_per_rep_writer(self, pop):
        repid = pop.dvars().rep
        aggr_frq_matrix = np.zeros((self.rows_per_rep, 44444), dtype=np.float16)
        subpops_and_gens = [(1, 0, 0), (2, 2, 1), (3, 4, 2), (4, 6, 3), (5, 8, 4), (6, 10, 5)]
        for sp, curgen, rowidx in subpops_and_gens:
            aggr_frq_matrix[rowidx, 0] = curgen
            aggr_frq_matrix[rowidx, 1] = repid
            aggr_frq_matrix[rowidx, 2:] = [pop.dvars(sp).alleleFreq[locus][ma] for ma, locus in
                                           zip(self.minor_allele_reference.values(), range(44442))]
        np.savetxt('metapop_replicate_{replicate_id}_minor_allele_frq.txt'.format(replicate_id=repid), aggr_frq_matrix,
                   fmt='%.3f', delimiter='\t')
        return True



class GenotypicEffectCalculator(sim.PyOperator):
    def __init__(self, absolute_qtl, ae_matrix, *args, **kwargs):
        self.absolute_qtl = absolute_qtl
        self.ae_matrix = ae_matrix
        sim.PyOperator.__init__(self, func=self.genotypic_effect_calculator, *args, **kwargs)

    def genotypic_effect_calculator(self, ind):
        """
        Calculates the genotypic effect by summing over all qtl triplets which are contained
        in the lineage values of an individual.
        """
        triplets = [ind.lineage()[qtl_idx] for qtl_idx in self.absolute_qtl]
        triplet_effects = [self.ae_matrix[idx-1, int(str(trip)[1])] + self.ae_matrix[idx, int(str(trip)[2])] +
                           self.ae_matrix[idx+1, int(str(trip)[3])] for idx, trip in zip(self.absolute_qtl, triplets)]
        ge = sum(triplet_effects)
        ind.ge = ge
        return True


class PhenotypicEffectCalculator(sim.PyOperator):
    def __init__(self, proportion_selected, *args, **kwargs):
        self.proportion_selected = proportion_selected
        sim.PyOperator.__init__(self, func=self.phenotypic_effect_calculator, *args, **kwargs)

    @staticmethod
    def phenotypic_effect_calculator(pop):
        # Sort population in ascending order
        pop.sortIndividuals('ge', reverse=True)
        ge = np.array(pop.indInfo('ge'))
        sigma_g = np.var(ge)
        if pop.dvars().h_squared_calc == 'fixed':
            heritability = pop.dvars().initial_heritability
            v_e = -1*sigma_g + sigma_g/heritability
            pop.dvars().v_e = v_e
            pop.setIndInfo([g_effect + random.normalvariate(0, v_e) for g_effect in ge], 'pe')
            return True
        elif pop.dvars().h_squared_calc == 'dynamic':
            heritability = pop.dvars().heritability
            v_e = -1*sigma_g + sigma_g/heritability
            pop.dvars().v_e = v_e
            pop.setIndInfo([g_effect + random.normalvariate(0, v_e) for g_effect in ge], 'pe')
            var_g = np.var(np.array(pop.indInfo('ge')))
            var_p = np.var(np.array(pop.indInfo('pe')))
            h_squared = var_g / var_p
            pop.dvars().heritability = h_squared
        return True


class CullPopulation(sim.PyOperator):
    def __init__(self, proportion_selected, *args, **kwargs):
        self.proportion_selected = proportion_selected
        sim.PyOperator.__init__(self, func=self.assign_fitness, *args, **kwargs)

    def assign_fitness(self, pop, subPop):
        pop.sortIndividuals('pe', reverse=True)
        pe = pop.indInfo('pe')
        cutoff_index = -1*(int(self.proportion_selected*pop.popSize()) + 1)
        cutoff = sorted(pe)[cutoff_index]
        pop.dvars().cutoff = cutoff
        pop.setIndInfo([x > cutoff for x in pe], 'fitness')
        return True


class GenerateMinorAlleleList(sim.PyOperator):
    """
    Operator to determine the minor allele at each locus.
    """
    def __init__(self, *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.minor_allele_list_generator, *args, **kwargs)

    @staticmethod
    def minor_allele_list_generator(pop):
        pre_minor_alleles = [pop.dvars().alleleFreq[i] for i in range(pop.totNumLoci())]
        minor_allele_map = map(min, pre_minor_alleles)
        minor_alleles = [m for m in minor_allele_map]
        pop.dvars().minorAlleleList = minor_alleles
        return True


class MinorAlleleCounter(sim.PyOperator):
    """
    Operator to create a minor allele frequency matrix for a single generation.
    """
    def __init__(self, *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.minor_allele_counter, *args, **kwargs)

    @staticmethod
    def minor_allele_counter(pop):
        minor_allele_list = pop.dvars().minorAlleleList
        # option to add the index of 'ind_id' -- may be more suitable for databases
        minor_allele_matrix = pd.DataFrame(np.zeros((pop.popSize(), pop.totNumLoci())))
        for i, ind in enumerate(pop.individuals()):
            cmp_one = np.equal(minor_allele_list, ind.genotype(0))
            cmp_two = np.equal(minor_allele_list, ind.genotype(1))
            minor_allele_count = np.add(cmp_one, cmp_two, dtype=int)
            minor_allele_matrix.ix[i, :] = minor_allele_count
            pop.dvars().minor_allele_matrix = minor_allele_matrix
        return True


class MinorAlleles(sim.PyOperator):
    def __init__(self, *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.minor_allele_finder, *args, **kwargs)

    def minor_allele_finder(self, pop):
        """
        Builds a dictionary of minor alleles.
        This operator only measures frequency of 'central' allele of each triplet.
        """
        allele_frequencies = pop.dvars().alleleFreq
        minor_allele_frequencies = col.OrderedDict()
        minor_alleles = col.OrderedDict()
        for locus in range(pop.totNumLoci()):
            minor_allele_frequencies[locus] = min(allele_frequencies[locus].values())
            for k in allele_frequencies[locus].keys():
                if allele_frequencies[locus][k] == minor_allele_frequencies[locus]:
                    minor_alleles[locus] = k
        pop.minor_alleles = minor_alleles
        return True


class MajorAlleles(sim.PyOperator):
    def __init__(self, *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.major_allele_finder, *args, **kwargs)

    def major_allele_finder(self, pop):
        allele_frequencies = pop.allele_frequencies
        major_allele_frequencies = col.OrderedDict()
        major_alleles = col.OrderedDict()
        for locus in range(pop.totNumLoci()):
            major_allele_frequencies[locus] = min(allele_frequencies[locus].values())
            for k in allele_frequencies[locus].keys():
                if allele_frequencies[locus][k] == major_allele_frequencies[locus]:
                    major_alleles[locus] = k
        pop.major_alleles = major_alleles
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


class CalcTripletFrequencies(sim.PyOperator):
    def __init__(self, *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.calc_triplet_frequencies, *args, **kwargs)

    def calc_triplet_frequencies(self, pop):
        lineage = np.array(pop.lineage()).reshape(pop.popSize(), 2*pop.totNumLoci())
        splitform = np.split(lineage, 2, 1)
        alpha = splitform[0]
        omega = splitform[1]
        triplet_frequencies = col.defaultdict(int)
        for locus in range(pop.totNumLoci()):
            alpha_slice = alpha[:, locus]
            omega_slice = omega[:, locus]
            combined = np.append(alpha_slice, omega_slice)
            triplet_frequencies[locus] = col.defaultdict(int, col.Counter(combined))
            for triplet in list(triplet_frequencies[locus].keys()):
                triplet_frequencies[locus][triplet] = triplet_frequencies[locus][triplet]/(2*pop.popSize())
            pop.dvars().tripletFreq = triplet_frequencies
        return True


class MetaPopulation(sim.PyOperator):
    """
    Adds individuals to a meta population.
    """
    def __init__(self, meta_population, sample_size, *args, **kwargs):
        self.meta_population = meta_population
        self.sample_size = sample_size
        sim.PyOperator.__init__(self, func=self.add_to_meta_pop, *args, **kwargs)

    def add_to_meta_pop(self, pop):
        sampled = sampling.drawRandomSample(pop, sizes=self.sample_size[pop.dvars().gen], subPops=sim.ALL_AVAIL)
        pop.dvars().ss = self.sample_size[pop.dvars().gen]
        pop.dvars().gen_sampled_from = pop.dvars().gen
        self.meta_population.addIndFrom(sampled)
        return True

class TestReplicateMetaPopulation(sim.PyOperator):
    """
    Operator which writes native simuPOP .pop files of each group of individuals sampled. This allows
    us to data in the final metapopulation against the data as it is sampled at run-time. This operator
    is to mainly address the concern that individuals are not being sampled from diferent generations through
    an undetected bug in the code.
    """
    def __init__(self, meta_replicates, sample_size, *args, **kwargs):
        self.meta_replicates = meta_replicates
        self.sample_size = sample_size
        sim.PyOperator.__init__(self, func=self.add_to_meta_pop, *args, **kwargs)

    def add_to_meta_pop(self, pop):
        rep_id = pop.dvars().rep
        sampled = sampling.drawRandomSample(pop, sizes=self.sample_size[pop.dvars().gen])
        sampled_file_name = 'sampled_rep_' + str(rep_id) + '_gen_' + str(pop.dvars().gen) + '_metapop.pop'
        sampled.save(sampled_file_name)
        pop.dvars().ss = self.sample_size[pop.dvars().gen]
        pop.dvars().gen_sampled_from = pop.dvars().gen
        self.meta_replicates.population(rep_id).addIndFrom(sampled)
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
        header = ['ind_id', 'mother_id', 'father_id', 'ge', 'pe']
        genotype_header = list(range(2*pop.totNumLoci()))
        header.extend(genotype_header)
        with open(full_file_name, 'w') as pop_info:
            info_writer = csv.writer(pop_info, delimiter=',')
            info_writer.writerow(header)
            for ind in pop.individuals():
                info_writer.writerow([ind.ind_id, ind.mother_id, ind.father_id, ind.ge, ind.pe, ind.genotype()])
        return True


class RandomlyAssignFemaleFitness(sim.PyOperator):
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
    Operator to choose <number_to_remove> individuals at random to remove
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


