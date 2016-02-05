# -*- coding: utf-8 -*-


import simuPOP as sim
import pandas as pd
import collections as col
import random
import numpy as np
import json
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


def seg_qtl_chooser(pop: sim.Population, loci_subset: list, number_qtl: int):
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

class ReadWrite(object):
    """
    Embodies functions to read and write parameter files for simulations to
    a Python shelve object. At present it is the easiest and fastest way to
    generate storeable data sets.
    """

    def write_trunc_selection_parameters(self, trunc_sel_parameters,
                                         truncation_selection_filename,
                                         qtl_parameters, qtl_filename,
                                         gen_struct_parameters,
                                         genetic_structure_filename):
        """Writes a json file for the parameters required for the simulation."""
        with open(truncation_selection_filename, 'w') as sel:
            json.dump(trunc_sel_parameters, sel, indent=2)
        with open(qtl_filename, 'w') as qf:
            json.dump(qtl_parameters, qf, indent=2)
        with open(genetic_structure_filename, 'w') as gs:
            json.dump(gen_struct_parameters, gs)
        print("Writing parameters to : {}, {}, {}.".format(truncation_selection_filename,
                                                                 qtl_filename,
                                                             genetic_structure_filename))

    def load_trunc_selection_parameters(self, truncation_selection_filename,
                                        qtl_filename,
                                        genetic_structure_filename):
        """Loads a json file for the parameters required for the simulation."""
        with open(truncation_selection_filename, 'r') as sel:
            selection_parameters = json.load(sel)
        with open(qtl_filename, 'r') as qf:
            qtl_parameters = json.load(qf)
        with open(genetic_structure_filename, 'r') as gs:
            gen_stru_parameters = json.load(gs)

        return selection_parameters, qtl_parameters, gen_stru_parameters