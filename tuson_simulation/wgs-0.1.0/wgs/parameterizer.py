# -*- coding: utf-8 -*-
__docformat__='restructuredtext en'

import simuPOP as sim
import pandas as pd
import itertools as ite
import collections as col
import random
import numpy as np
from scipy import stats



class PopulationStructure(object):

    def __init__(self, pop, population_structure_matrix_filename, threshold, error):
        self.pop = pop
        self.population_structure_matrix_filename = population_structure_matrix_filename
        self.threshold = threshold
        self.error = error

    def generate_population_structure(self):
        """
        Parses a population structure matrix from a file and converts it into a native Python dictionary.
        Extra steps needed because the population structure file is not in the same order as the genotype file.
        """
        popst = pd.read_excel(self.population_structure_matrix_filename)
        indid_to_sampleid = {popst.index[i]: popst['sample_id'][i] for i in range(len(popst))}
        popst_proportions = {popst['sample_id'][i]: list(popst.ix[i, 1:7]) for i in range(len(popst))}
        structure = {indid_to_sampleid[i]: popst_proportions[indid_to_sampleid[i]] for i in range(len(popst))}
        return structure

    def population_structure_filter(self, population_structure: dict):
        """
        In populations with population structure spread over several subpopulations takes individuals whose genomes are
        derived from at most 2 subpopulations plus/minus error. This function also removes the invalid individuals
        referenced by absolute index. Resulting population is initial pop size - len(invalids.keys())
        """
        invalids = col.defaultdict()
        for k, v in population_structure.items():
            primary = max(v)
            secondary = sorted(v)[-2]
            if len(set(v)) > 3:
                if (sum(v) - primary - secondary) > (self.threshold + self.error):
                    invalids[k] = v
                else:
                    pass
        self.pop.removeIndividuals(indexes=list(invalids.keys()))
        return list(invalids.keys())

    def assign_population_structure(self, population_structure):
        """
        Assigns every individual a primary and secondary subpopulation according to the population structure matrix.
        :param population_structure_dict:
        :param number_of_individuals:
        :return:
        """
        assigned_structure = col.OrderedDict()
        for ind in list(self.pop.indInfo('ind_id')):
            assigned_structure[int(ind)] = col.OrderedDict()
            for i, prop in enumerate(population_structure[int(ind)]):
                if max(population_structure[int(ind)]) == prop:
                    assigned_structure[int(ind)][0] = i
                elif population_structure[int(ind)][i] == sorted(population_structure[int(ind)])[-2]:
                    assigned_structure[int(ind)][1] = i
                else:
                    pass
            if min(population_structure[int(ind)]) == 0.000000:
                assigned_structure[int(ind)][1] = assigned_structure[int(ind)][0]
        return assigned_structure

    def assign_structured_mating_probabilities(self, population_structure, assigned_primary_secondary_structure):
        """
        Sums the proportions which are non-primary and non-secondary and adds half of the sum to the primary and secondary.
        The primary and secondary probabilities is the probability of mating with an individual from that primary
        subpopulation (or selfing).
        :param population_structure_dict:
        :param assigned_structure_dict:
        :return:
        """
        mating_probabilities = col.OrderedDict()
        for ind in self.pop.indInfo('ind_id'):
            ind = int(ind)
            mating_probabilities[ind] = col.OrderedDict()
            primary_proportion = population_structure[ind][assigned_primary_secondary_structure[ind][0]]
            secondary_proportion = population_structure[ind][assigned_primary_secondary_structure[ind][1]]
            remainder = (1 - primary_proportion - secondary_proportion)/2
            mating_probabilities[ind][0] = primary_proportion + remainder
            mating_probabilities[ind][1] = secondary_proportion + remainder
        return mating_probabilities

    def generate_mating_probability_mass_functions(self, assigned_primary_secondary_structure: dict,
                                                   assigned_mating_probabilities: dict):
        """
        Assigns mating probabilities as dictionary values keyed by individual ids.
        :param pop:
        :param population_structure:
        :param assigned_structure:
        :return:
        """
        mating_probability_mass_functions = col.OrderedDict()
        for ind in list(self.pop.indInfo('ind_id')):
            ind_id = int(ind)
            if assigned_primary_secondary_structure[ind_id][0] == assigned_primary_secondary_structure[ind_id][1]:
                single_subpopulation = (assigned_primary_secondary_structure[ind_id][0])
                mating_probability = 1.0
                mating_probability_mass_functions[ind_id] = stats.rv_discrete(values=(single_subpopulation,
                                                                                      mating_probability))
            else:
                primary_and_secondary_subpopulations = (assigned_primary_secondary_structure[ind_id][0],
                                                        assigned_primary_secondary_structure[ind_id][1])
                mating_probabilities = (float(assigned_mating_probabilities[ind_id][0]),
                                        float(assigned_mating_probabilities[ind_id][1]))
                mating_probability_mass_functions[ind_id] = stats.rv_discrete(
                    values=(primary_and_secondary_subpopulations, mating_probabilities))
        return mating_probability_mass_functions

    def setup_mating_structure(self):
        """
        Function which applies all the functions necessary to create the mating_pmfs. The mating_pmfs are
        assigned to self.pop's local namespace where they can be accessed by a parent_chooser function.
        """
        pop_structure = self.generate_population_structure()
        self.pop.dvars().invalids = self.population_structure_filter(pop_structure)
        self.pop.dvars().assigned_structure = self.assign_population_structure(pop_structure)
        self.pop.dvars().mating_probabilities = self.assign_structured_mating_probabilities(pop_structure,
                                                                                            self.pop.dvars().assigned_structure)
        self.pop.dvars().mating_pmfs = self.generate_mating_probability_mass_functions(
            self.pop.dvars().assigned_structure, self.pop.dvars().mating_probabilities,
        )


class Parameterizer(object):
    @staticmethod
    def nucleotide_translator(nucleotide):
        """
        :param nucleotide:
        In some cases a genotype matrix file will list SNP markers as nucleotides. nucleotide_translator is a basic
        callback function to convert the string nucleotides into integers suitable for use in simuPOP.
        """
        if nucleotide == 'A':
            return 0
        elif nucleotide == 'C':
            return 1
        elif nucleotide == 'G':
            return 2
        elif nucleotide == 'T':
            return 3
        elif nucleotide == '-':
            return 4
        elif nucleotide == '+':
            return 5

    def genotype_setter(self, pop, map_file_name):
        genotype_table = pd.read_table(map_file_name, index_col=0)
        central_marker_indices = [i for i in range(len(genotype_table))
                                  if float(genotype_table['cM_pos'][i]) == float(int(genotype_table['cM_pos'][i]))]
        founder_names = genotype_table.columns[3:]
        pop.dvars().founderNames = founder_names

        all_genotypes_handle = [
            [self.nucleotide_translator(genotype_table.ix[idx, name][0]) for idx in central_marker_indices] +
            [self.nucleotide_translator(genotype_table.ix[idx, name][1]) for idx in central_marker_indices]
            for name in founder_names]

        all_lineages_handle = [[int(str(1) + str(self.nucleotide_translator(genotype_table.ix[idx - 1, name][0])) +
                                    str(self.nucleotide_translator(genotype_table.ix[idx, name][0])) +
                                    str(self.nucleotide_translator(genotype_table.ix[idx + 1, name][0]))) for idx in
                                central_marker_indices] +
                               [int(str(1) + str(self.nucleotide_translator(genotype_table.ix[idx - 1, name][1])) +
                                    str(self.nucleotide_translator(genotype_table.ix[idx, name][1])) +
                                    str(self.nucleotide_translator(genotype_table.ix[idx + 1, name][1]))) for idx in
                                central_marker_indices] for name in founder_names]

        for ind, genotype, lineotype in zip(pop.individuals(), all_genotypes_handle, all_lineages_handle):
            ind.setGenotype(genotype)
            ind.setLineage(lineotype)


class AE(object):
    """
    Allele effects of the simulator are assigned as random draws of a statistical distribution. The user specifies the
    type of distribution as a string when creating an instance of the class. The types of distributions are:
    exponential, normal, poisson
    """

    def exponential(self, pop, parameter_of_exponential):
        """

        :param pop:
        :param parameter_of_exponential:
        :return:
        Draws allele effects an exponential distribution.
        Creates a copy of the AE map convienient for plotting.
        """
        allele_effects = col.OrderedDict()
        for_plot_allele_effects = col.OrderedDict()
        for idx in pop.dvars().properQTL:
            idxtwo = idx + pop.totNumLoci()
            for nucleotide in range(6):
                allele_effects[idx - 1, nucleotide] = random.expovariate(parameter_of_exponential)
                allele_effects[idx, nucleotide] = random.expovariate(parameter_of_exponential)
                allele_effects[idx + 1, nucleotide] = random.expovariate(parameter_of_exponential)
                allele_effects[idxtwo - 1, nucleotide] = allele_effects[idx - 1, nucleotide]
                allele_effects[idxtwo, nucleotide] = allele_effects[idx, nucleotide]
                allele_effects[idxtwo + 1, nucleotide] = allele_effects[idx + 1, nucleotide]
                for_plot_allele_effects[float(idx) - 0.2, nucleotide] = allele_effects[idx - 1, nucleotide]
                for_plot_allele_effects[float(idx), nucleotide] = allele_effects[idx, nucleotide]
                for_plot_allele_effects[float(idx) + 0.2, nucleotide] = allele_effects[idx + 1, nucleotide]
        return allele_effects, for_plot_allele_effects


class QTL(object):
    """
    Class which has several different QTL choosers.
    """

    @staticmethod
    def chooser(pop, number_qtl):
        """
        Chooses QTL from all possible loci
        :param pop:
        :param number_qtl:
        :return:
        """
        alpha_qtl = sorted(random.sample(list(range(1, pop.totNumLoci() - 1)), number_qtl))
        omega_qtl = sorted([pop.totNumLoci() + qtl_index for qtl_index in alpha_qtl])
        all_qtl = alpha_qtl + omega_qtl
        proper_qtl = alpha_qtl
        return all_qtl, proper_qtl

    @staticmethod
    def seg_chooser(pop, number_qtl):
        """
        Chooses QTL from segregating loci.
        :param pop:
        :param number_qtl:
        :return:
        """
        sim.stat(pop, numOfSegSites=sim.ALL_AVAIL, vars=['segSites'])
        alpha_qtl = sorted(random.sample(pop.dvars().segSites, number_qtl))
        omega_qtl = sorted([pop.totNumLoci() + qtl_index for qtl_index in alpha_qtl])
        all_qtl = alpha_qtl + omega_qtl
        proper_qtl = alpha_qtl
        return all_qtl, proper_qtl

    @staticmethod
    def fixed_chooser(pop, number_qtl):
        """
        Chooses QTL from only fixed loci.
        :param pop:
        :param number_qtl:
        :return:
        """
        sim.stat(pop, numOfSegSites=sim.ALL_AVAIL, vars=['fixedSites'])
        alpha_qtl = sorted(random.sample(pop.dvars().fixedSites, number_qtl))
        omega_qtl = sorted([pop.totNumLoci() + qtl_index for qtl_index in alpha_qtl])
        all_qtl = alpha_qtl + omega_qtl
        proper_qtl = alpha_qtl
        return all_qtl, proper_qtl


class LD(object):
    """
    Class which has methods to handle anything dealing with linkage
    disequilibrium.
    """

    @staticmethod
    def loci_pairs(pop):
        """
        Creates a list of pairs of loci for all chromosomes. Specific to
        maize genome because the chromosome number is 10.
        """
        chromosome_loci_lists = [list(range(pop.chromBegin(i), pop.chromEnd(i))) for i in range(10)]
        loci_pairs = []
        for i in range(10):
            loci_pairs.extend(list(ite.combinations(chromosome_loci_lists[i], 2)))
        return loci_pairs

    @staticmethod
    def loci_pairs_replacement(pop):
        chrom_locs = [list(range(pop.chromBegin(i), pop.chromEnd(i))) for i in range(10)]
        loci_pairs_replacement = []
        for i in range(10):
            loci_pairs_replacement.extend(list(ite.combinations_with_replacement(chrom_locs[i], 2)))
        return loci_pairs_replacement

    @staticmethod
    def qtl_pairs(pop):
        list_of_qtl = pop.dvars().properQTL
        qtl_pairs = list(ite.combinations_with_replacement(list_of_qtl, 2))
        return qtl_pairs


class GenotypeData(object):
    """
    Code not in working state.
    10/14/15
    """
    def __init__(self, genotype_matrix_filename):
        self.genotype_matrix_filename = genotype_matrix_filename

    def parse_genotype_matrix(self, columns_to_drop='popdata'):
        genotype_matrix = pd.read_csv(self.genotype_matrix_filename, sep='\t', index_col=0, low_memory=False)
        droppable_individuals = list(genotype_matrix.index[105:])
        genotype_matrix = genotype_matrix.drop(droppable_individuals, axis=0)
        genotype_matrix = genotype_matrix.drop(columns_to_drop, axis=1)
        return genotype_matrix

    def genotype_counts_to_frequencies(self, genotype_counts: dict, missing_loci: list):
        """
        Converts a the dictionaries of genotype: count for each locus into their
        frequency equivalents by dropping and missing data and dividing by the adjusted
        total.

        :param genotype_counts:
        :param missing_loci:
        :type genotype_counts str:


        :param missing_loci:
        :return:
        """
        geno_frq = {}
        for mlocus in missing_loci:
            geno_frq[mlocus] = {}
            if np.nan in genotype_counts[mlocus]:
                del genotype_counts[mlocus][np.nan]
            inds_counted = sum(genotype_counts[mlocus].values())
            for genotype, cnt in genotype_counts[mlocus].items():
                geno_frq[mlocus][genotype] = cnt/inds_counted
        return geno_frq

    @staticmethod
    def centralized_genotype_pmfs(genotype_frequencies):
        """
        For the time being all of the information required to compute a custom
        probability mass function for each locus is stored a dictionary keyed by locus.
        The values are tuples:
        0: genotype: frequency
        1: integer: genotype
        2: density
        3: genotype: integer
        """
        centralized_pmfs = col.OrderedDict()
        for locus, frq_map in genotype_frequencies.items():
            pre_density = {genotype: frequency for genotype, frequency in frq_map.items()}
            genotype_to_int_map = {genotype: i for i, genotype in list(enumerate(frq_map.keys()))}
            density = {genotype_to_int_map[genotype]: frequency for genotype, frequency in frq_map.items()}
            int_to_genotype_map = {i: genotype for i, genotype in list(enumerate(frq_map.keys()))}
            centralized_pmfs[locus] = (pre_density, genotype_to_int_map, density, int_to_genotype_map)
        return centralized_pmfs

    @staticmethod
    def individual_missing_data(genotype_matrix):
        """
        Each individual has a particular set of loci for which they are missing data. For each individual we need
        to know what loci are missing. Given the missing locus we can replace the 'NA' with a random draw
        from the genotype pmf of that locus.
        :param genotype_matrix:
        """
        nan_dict = {}
        nan_array = np.array(pd.isnull(genotype_matrix))
        for individual, row in enumerate(nan_array):
            nan_dict[individual] = [locus for locus, val in enumerate(row) if val == True]
        return nan_dict

    @staticmethod
    def replace_missing_genotypes(genotype_matrix, population_genotype_pmfs):
        """
        A function to replace each individuals missing genotype data with random draws from a dictionary of
        genotype pmfs. Parameter missing_loci_per_individual is a dictionary of individual: list_of_missing_loci pairs.
        population_genotype_pmfs is a nested dictionary which provides all the necessary mapping data to create the
        replacement data.
        Note: Assumes that genotype_matrix has rows=individuals and columns=genotypes.
        """
        for ind in range(genotype_matrix.shape[0]):
            individuals_missing_loci = [genotype_matrix[ind, i] for i in range(genotype_matrix.shape[1])
                                        if genotype_matrix[ind, i] == np.nan]
            for locus in individuals_missing_loci:
                integer_genotype = population_genotype_pmfs[locus]['pmf'].rvs()
                geno_state = population_genotype_pmfs[locus]['integer_to_state'][integer_genotype]
                genotype_matrix[ind, locus] = geno_state
        return genotype_matrix
