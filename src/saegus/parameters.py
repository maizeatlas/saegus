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
                mating_probability_mass_functions[ind_id] = \
                    stats.rv_discrete(values=(single_subpopulation,
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

    """

    def __init__(self):
        pass

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
                                  *distribution_function_parameters,
                                  multiplicity):

        allele_effects = {}
        for locus in qtl:
            allele_effects[locus] = {}
            for allele in alleles[locus]:
                allele_effects[locus][allele] = \
                    sum([distribution_function(*distribution_function_parameters)
                         for i in range(multiplicity)])
        return allele_effects

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

def count_segregating_site_concordance(array_of_seg_loci):
    segregating_loci_concordance_counts = col.defaultdict(int, default=0)
    for row in array_of_seg_loci:
        segregating_loci_concordance_counts[tuple(row)] += 1
    return segregating_loci_concordance_counts

def test_segregatng_loci_concordance(seg_loci_agreement_counts):
    seg_loci_agreement = True
    if len(seg_loci_agreement_counts) > 2:
        seg_loci_agreement = False
    return seg_loci_agreement


def randomly_convert_fixed_sites(pop, fixed_sites):

    """
    Randomly converts fixed sites in pop to
    nonfixed_sites by changing the allele state at that site.
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

