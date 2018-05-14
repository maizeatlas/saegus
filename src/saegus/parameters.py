# -*- coding: utf-8 -*-


import simuPOP as sim
import pandas as pd
import collections as col
import random
import numpy as np
import copy
import yaml
from scipy import stats



class PopulationStructure(object):

    def __init__(self, pop):
        self.pop = pop

    def parse_and_remap_population_structure(self, population_struct_matrix_file_name):
        """
        :parameter str population_structure_matrix_file_name: Pop strct filename

        Parses a population structure matrix from a file and converts it into
        Python dictionary. The :file:`population_structure_matrix.xlsx` file is
        not in the same order as :file:`genotype_matrix.txt`. This function
        remaps the inheritance proportions to it matches the genotype
        matrix of the population.

        """
        popst = pd.read_excel(population_struct_matrix_file_name)
        indid_to_sampleid = {popst.index[i]: popst['sample_id'][i] for i in range(len(popst))}
        popst_proportions = {popst['sample_id'][i]: list(popst.ix[i, 1:7]) for i in range(len(popst))}
        structure = {indid_to_sampleid[i]: popst_proportions[indid_to_sampleid[i]] for i in range(len(popst))}
        return structure

    def generate_mating_pmfs(self, population_structure_dict):
        """
        Converts a dictionary of lists of probabilities into scipy.rv_discrete
        customized probability mass functions.

        :param population_structure_dict: Dict of lists of probabilities
        :return: Dictionary of scipy probability mass functions
        """
        mating_pmfs = {}
        for ind, probabilities in population_structure_dict.items():
            for i, prob in enumerate(probabilities):
                values = []
                probabilites = []
                for i, prob in enumerate(population_structure_dict[ind]):
                    values.append(i)
                    probabilites.append(prob)
                pmf_values = (values, probabilites)
                mating_pmfs[ind] = stats.rv_discrete(values=pmf_values)
        return mating_pmfs

    def assign_primary_subpopulations(self, struct_mating_probabilities):
        """

        :param struct_mating_probabilities: Dict of lists of probabilities

        Assigns the primary subpopulation to each individual according to
        ``ind_id``. Primary subpopulation is the population from which the
        individual derives most of its genome.

        """
        primary_subpop = {}
        for ind_id, inheritance_proportions in struct_mating_probabilities.items():
            primary_subpop[ind_id] = float(np.argmax(inheritance_proportions))
        for ind in self.pop.individuals():
            ind.primary = primary_subpop[ind.ind_id]

    # todo Documentation for correct_rounding_error

    def correct_rounding_error(self, structure_array):
        """
        Rows of ``structure_array`` may not sum to 1 due to rounding error or
        error propagation. Examines the rows which do not sum to 1 and adds
        the small difference to the largest proportion.

        :return:
        """

        for i in range(structure_array.shape[0]):
            if sum(structure_array[i]) < 1:
                greatest_probability_index = np.argmax(structure_array[i])
                difference = 1 - sum(structure_array[i])
                structure_array[i, greatest_probability_index] = \
                    structure_array[i, greatest_probability_index] + difference
            if sum(structure_array[i]) > 1:
                greatest_probability_index = np.argmax(structure_array[i])
                difference = sum(structure_array[i]) - 1
                structure_array[i, greatest_probability_index] = \
                    structure_array[i, greatest_probability_index] - difference

        return structure_array



class MissingGenotypeData(object):
    """
    A class to handle all genotype data.
    """

    def __init__(self, raw_genotype_array, number_of_individuals,
                 number_of_markers, missing_genotype_token):

        if raw_genotype_array.shape[0] == number_of_markers:
            raw_genotype_array = raw_genotype_array.T
            assert raw_genotype_array.shape[0] == number_of_individuals, \
                "Genotype data does not have proper shape"

        self.genotype_array = raw_genotype_array
        self.number_of_individuals = number_of_individuals
        self.number_of_markers = number_of_markers
        self.missing_genotype_token = missing_genotype_token

    def genotype_counter(self, genos):
        """

        :param missing_genotype_token:
        :return:
        """
        genotype_counts = {}
        for genos, locus in zip(self.genotype_array.T, range(
                self.number_of_markers)):
            genotype_counts[locus] = dict(col.Counter(genos))

        loci_missing_data = []
        for locus, gcounts in genotype_counts.items():
            if self.missing_genotype_token in gcounts.keys():
                loci_missing_data.append(locus)
                del gcounts[self.missing_genotype_token]

        return genotype_counts, loci_missing_data

    def determine_individuals_missing_data(self, loci_missing_data,
                                           missing_genotype_token):
        """
        Finds the individuals who are missing genotype data given the list of
        for which data are missing.
        :param loci_missing_data: List of loci which have missing genotype data
        :param missing_genotype_token: A string representing missing data
        :return: Dictionary of lists of individual indexes
        """
        location_of_missing_values = {locus: [] for locus in
                                      loci_missing_data}
        for locus in loci_missing_data:
            for idx, ind in enumerate(self.genotype_array):
                if ind[locus] == missing_genotype_token:
                    location_of_missing_values[locus].append(idx)

        return location_of_missing_values

    def convert_counts_to_frq(self, genotype_counts):
        """
        Converts a the dictionaries of genotype: count for each locus into their
        frequency equivalents by dropping and missing data and dividing by the adjusted
        total.

        :param genotype_counts:
        :return: Dictionary of genotype frequencies
        """
        genotype_frqs = {locus: {} for locus in range(self.number_of_markers)}
        for locus, g_counts in genotype_counts.items():
            markers_counted = sum(g_counts.values())
            for genotype, count in g_counts.items():
                genotype_frqs[locus][genotype] = \
                genotype_counts[locus][genotype]/markers_counted
        return genotype_frqs

    def generate_pmf_mappings(self, genotype_frequencies):
        """
        Collects all information for creating genotype probability mass
        functions to replace missing genotype data.
        :param genotype_frequencies: Dictionary[locus][genotype]: frequency
        :return:
        """
        pmf_mappings = \
            {locus: {} for locus in genotype_frequencies.keys()}

        for locus, frq_map in genotype_frequencies.items():
            genotype_to_int_map = \
                {genotype: i for i, genotype in list(enumerate(frq_map.keys()))}
            int_to_genotype_map = \
                {i: genotype for i, genotype in list(enumerate(frq_map.keys()))}
            int_to_frq_map = \
                {genotype_to_int_map[genotype]:
                     frequency for genotype, frequency in frq_map.items()}
            pmf_mappings[locus]['geno_to_int'] = \
                genotype_to_int_map
            pmf_mappings[locus]['int_to_geno'] = \
                int_to_genotype_map
            pmf_mappings[locus]['int_to_frq'] = \
                int_to_frq_map
        return pmf_mappings

    def generate_genotype_pmfs(self, empirical_pmf_mappings):
        """
        A set of nested dictionaries which enable the use of
        scipy.stats.rv_discrete to generate a genotype probability mass
        functions.
        :param empirical_pmf_mappings:
        :return:
        """
        empirical_genotype_pmfs = {}
        for locus in range(self.number_of_markers):
            xk = list(empirical_pmf_mappings[locus]['int_to_frq'].keys())
            pk = list(empirical_pmf_mappings[locus]['int_to_frq'].values())
            empirical_genotype_pmfs[locus] = \
                stats.rv_discrete(values=(xk, pk))

        return empirical_genotype_pmfs


    def replace_missing_genotypes(self):
        """
        A function to utilize all of :class:`Genotype`'s methods to fill in
        missing genotype values. :func:`replace_missing_genotypes` uses
        empirically observed genotype frequencies.


        :param population_genotype_pmfs:
        :param locations_missing_data:
        :param pmass_func_mappings:
        """
        genotype_counts, loci_missing_data = self.genotype_counter()
        locations_of_missing_data = self.determine_individuals_missing_data(
            loci_missing_data, self.missing_genotype_token)
        empirically_observed_frequencies = \
            self.convert_counts_to_frq(genotype_counts)
        pmf_mappings = \
            self.generate_pmf_mappings(empirically_observed_frequencies)
        genotype_pmfs = self.generate_genotype_pmfs(pmf_mappings)


        for locus, individuals in locations_of_missing_data.items():
            for ind in individuals:
                int_genotype = genotype_pmfs[locus].rvs()
                genotype = \
                    pmf_mappings[locus]['int_to_geno'][int_genotype]
                self.genotype_array[ind, locus] = genotype

        return self.genotype_array


class Trait(object):
    """
    This class carries functions responsible for assigning and handling
    trait models: Decision making process simulator uses to assign phenotypes.

    For the time being we assume certain things about the integers representing
    allele states and the ploidy of the organism.
    Ploidy = 2 and allele states coded to integers 1, 2, 3, 4

    6/20/2017 JJD

    """

    def __init__(self):
        pass

    # todo Create documentation for construct_allele_effects_table

    def construct_allele_effects_table(self, pop: sim.Population, qtl: list,
                                       distribution_function,
                                       *distribution_function_parameters):

        allele_effects_table = np.zeros((pop.totNumLoci(), 5))
        alpha_alleles = []
        omega_alleles = []

        sim.stat(pop, alleleFreq=sim.ALL_AVAIL)

        for locus in range(pop.totNumLoci()):
            alpha_alleles.append(list(pop.dvars().alleleFreq[locus])[0])
            omega_alleles.append(list(pop.dvars().alleleFreq[locus])[-1])
            if alpha_alleles[locus] == omega_alleles[locus]:
                omega_alleles[locus] = 0

        allele_effects_table[:, 0] = list(range(pop.totNumLoci()))
        allele_effects_table[:, 1] = alpha_alleles
        allele_effects_table[:, 3] = omega_alleles

        for locus in qtl:
            allele_effects_table[locus, 2] = \
                distribution_function(*distribution_function_parameters)
            allele_effects_table[locus, 4] = \
                distribution_function(*distribution_function_parameters)

        return allele_effects_table

    # todo Create documentation with example for construct_ae_array

    def construct_ae_array(self, allele_effects_table, qtl):
        """
        Conversion of allele effects table into an array made expressly for
        the purpose of computing ``g`` and ``p`` values. Each column of the
        array corresponds to the allele at that locus. It makes computation
        of genotypic value and phenotypic value very fast.
        """
        allele_effects_array = np.zeros(allele_effects_table.shape)
        for row in allele_effects_table[qtl]:
            allele_effects_array[int(row[0]), int(row[1])] = row[2]
            allele_effects_array[int(row[0]), int(row[3])] = row[4]

        return allele_effects_array

    # todo Create example use for construct_minor_major_effects

    def construct_minor_major_effects(self, allele_state_array,
                                      allele_effects_table,
                                      allele_effects_array,
                                      qtl):
        """
        Defines allele effects at each locus in terms of the minor/major
        alleles. ``allele_state_array`` assumes that allele state data is
        organized in particular way.
        Columns of ``allele_state_array``: locus, alpha, omega, minor, major

        Returns array which has columns:
        locus, minor_allele_state, minor_allele_effects, major_allele_state,
        major_allele_effect

        Allows for easy integration with other tabular data by joining on the
        locus column.

        :param allele_state_array: Array which gives minor and major alleles by locus
        :param qtl: Loci which have non-zero allele effects
        :return: Array which relates minor and major allele states to their effect
        """

        minor_major_allele_effects = np.zeros((allele_effects_array.shape[0],
                                               5))

        minor_alleles = allele_state_array[:, 3]
        major_alleles = allele_state_array[:, 4]

        minor_major_allele_effects[:, 0] = list(range(
            allele_effects_array.shape[0]))
        minor_major_allele_effects[:, 1] = minor_alleles
        minor_major_allele_effects[:, 3] = major_alleles

        minor_alpha_loci = np.where(allele_effects_table[:, 1] ==
                                    minor_alleles)[0]
        minor_omega_loci = np.where(allele_effects_table[:, 3] ==
                                    minor_alleles)[0]

        major_alpha_loci = np.where(allele_effects_table[:, 1] ==
                                       major_alleles)[0]
        major_omega_loci = np.where(allele_effects_table[:, 3] ==
                                       major_alleles)[0]

        for locus in qtl:
            if locus in minor_alpha_loci:
                minor_major_allele_effects[locus, 2] = \
                    allele_effects_array[locus, np.int(minor_alleles[locus])]
            if locus in minor_omega_loci:
                minor_major_allele_effects[locus, 2] = \
                    allele_effects_array[locus, np.int(minor_alleles[locus])]
            if locus in major_alpha_loci:
                minor_major_allele_effects[locus, 4] = \
                    allele_effects_array[locus, np.int(major_alleles[locus])]
            if locus in major_omega_loci:
                minor_major_allele_effects[locus, 4] = \
                    allele_effects_array[locus, np.int(major_alleles[locus])]

        return minor_major_allele_effects







    def seg_qtl_chooser(self, pop: sim.Population, loci_subset: list, number_qtl: int):
        """
        Chooses a random sample of ``number_qtl`` loci to be designated as QTL.
        Only chooses from among loci which are segregating in ``pop``.
        Determines which loci are segregating in ``loci_subset``.
        ``loci_subset`` can be all loci or a subset of them.
        :param number_qtl:
        :type number_qtl:
        :return:
        :rtype:
        """
        sim.stat(pop, numOfSegSites=loci_subset, vars=['numOfSegSites',
                                                       'numOfSegSites_sp',
                                                       'segSites', 'segSites_sp'])

        permissible_qtl = [locus for locus in pop.dvars().segSites if locus in
                           loci_subset]

        qtl = sorted(random.sample(permissible_qtl, number_qtl))
        return qtl

    def load_alleles(self, allele_file_name):
        """
        Loads alleles for all loci formatted in a way appropriate for assign_allele_effects.


        :param str allele_filename: HDF Filename containing alleles at each locus.
        :return:
        """

        alleles = np.array(pd.read_hdf(allele_file_name))
        return alleles

    def assign_allele_effects(self, alleles, qtl, distribution_function,
                                  *distribution_function_parameters):

        allele_effects = {}
        for locus in qtl:
            allele_effects[locus] = {}
            for allele in alleles[locus]:
                allele_effects[locus][allele] = \
                    distribution_function(*distribution_function_parameters)
        return allele_effects

    def convert_allele_effects_into_array(self, total_number_loci,
                                          total_number_alleles, allele_effects):
        """
        Convenience function to turn an allele effect dictionary to an array
        where the row is the locus and column corresponds to allele state.
        Locate effect by allele_effect_array[locus, allele_at_locus].
        """
        allele_effects_array = np.zeros(
            (total_number_loci, total_number_alleles))
        for locus, qt_alleles in allele_effects.items():
            for allele, effect in qt_alleles.items():
                allele_effects_array[locus, allele] = effect

        return allele_effects_array


    def assign_geometric_series(self, allele_effects, base, power):
        """
        Assigns the terms of a geometric series determined by base and power to
        one set of alleles at qtl. The alternate alleles are assigned effect 0.

        :param allele_effects:
        :param base:
        :param power:
        :return:
        """

        geometric_series = [base ** i_power for i_power in range(1, power+1)]
        geometric_allele_effects = copy.deepcopy(allele_effects)
        for locus, i in zip(sorted(allele_effects.keys()),
                            range(len(geometric_series))):
            ael = list(allele_effects[locus])
            geometric_allele_effects[locus][ael[0]] = geometric_series[i]
            geometric_allele_effects[locus][ael[1]] = 0
        return geometric_allele_effects




def count_qtl_concordance(array_of_seg_loci, qtl):
    qtl_agreement_counts = col.defaultdict(int, default=0)
    for i, row in enumerate(array_of_seg_loci):
        for locus in qtl:
            if locus in row:
                qtl_agreement_counts[i] += 1
    return qtl_agreement_counts

def test_qtl_concordance(agreement_counts, qtl):
    qtl_concordance = True
    for k, v in agreement_counts.items():
        if v != len(qtl) and k != 'default':
            qtl_concordance = False
            print("Disagrement of QTL at sample {}".format(k))
    return qtl_concordance

def async (array_of_seg_loci):
    segregating_loci_concordance_counts = col.defaultdict(int, default=0)
    for row in array_of_seg_loci:
        segregating_loci_concordance_counts[tuple(row)] += 1
    return segregating_loci_concordance_counts

def test_segregating_loci_concordance(seg_loci_agreement_counts):
    seg_loci_agreement = True
    if len(seg_loci_agreement_counts) > 2:
        seg_loci_agreement = False
    return seg_loci_agreement


def randomly_convert_fixed_sites(pop, fixed_sites, alleles=(0, 1, 2, 3)):

    """
    Randomly converts fixed sites in pop to
    nonfixed_sites by changing the allele state at that site.
    """
    random.shuffle(fixed_sites)
    for site in fixed_sites:
        random_id = random.choice(pop.indInfo('ind_id'))
        random_individual = pop.indByID(random_id)
        current_allele_state = random_individual.allele(site)
        possible_replacements = [allele for allele in alleles if
                                 allele != current_allele_state]
        replacement_allele = random.choice(possible_replacements)
        random_individual.setAllele(replacement_allele, site)

