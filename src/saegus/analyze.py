# -*- coding: utf-8 -*-
import simuPOP as sim
import math
import numpy as np
from numpy.core import defchararray
import pandas as pd
import collections as col
import os
import copy
import random
import shelve
import h5py
from scipy import linalg
from functools import singledispatch
import xml.etree.ElementTree as ET
import lxml.etree as etree
from . import operators, parameters


def gather_allele_data(pop):
    """
    Constructs a numpy.array with columns:
    locus   alpha   omega   minor   major

    Loci at 0.5 frequency have the minor allele set as the alpha allele
    and the major allele set as the omega allele

    :param sim.Population pop: diploid simuPOP population
    :return: Array labeling alpha, omega, minor and major alleles
    """

    allele_table = np.zeros((pop.totNumLoci(), 5))
    alpha_alleles = []
    omega_alleles = []

    for locus in range(pop.totNumLoci()):
        alpha_alleles.append(list(pop.dvars().alleleFreq[locus])[0])
        omega_alleles.append(list(pop.dvars().alleleFreq[locus])[-1])
        if alpha_alleles[locus] == omega_alleles[locus]:
            omega_alleles[locus] = 0

    minor_table = np.ones((pop.totNumLoci(), 5))
    major_table = np.zeros((pop.totNumLoci(), 5))

    for locus, alpha, omega in zip(range(pop.totNumLoci()), alpha_alleles,
                                   omega_alleles):
        minor_table[locus, alpha] = pop.dvars().alleleFreq[locus][alpha]
        minor_table[locus, omega] = pop.dvars().alleleFreq[locus][omega]
        major_table[locus, alpha] = pop.dvars().alleleFreq[locus][alpha]
        major_table[locus, omega] = pop.dvars().alleleFreq[locus][omega]

    minor_alleles = np.array(
        [np.argmin(minor_table[locus]) for locus in range(pop.totNumLoci())])
    major_alleles = np.array(
        [np.argmax(major_table[locus]) for locus in range(pop.totNumLoci())])
    tied_loci = np.where(minor_alleles == major_alleles)
    minor_alleles[tied_loci[0]] = np.array(alpha_alleles)[tied_loci[0]]
    major_alleles[tied_loci[0]] = np.array(omega_alleles)[tied_loci[0]]

    assert sum(minor_alleles == major_alleles) == 0, "At least one allele is "\
                                                        "classified as a minor"\
                                                        "and major allele."

    allele_table[:, 0] = np.array(range(pop.totNumLoci()))
    allele_table[:, 1] = alpha_alleles
    allele_table[:, 2] = omega_alleles
    allele_table[:, 3] = minor_alleles
    allele_table[:, 4] = major_alleles

    return allele_table

# todo Modify user_guide to reflect the separation of steps

def gather_allele_frequencies(pop, allele_state_table):
    """
    Constructs an array of allele frequencies with columns:
    locus   alpha_frequency omega_frequency minor_frequency major_frequency


    :param sim.Population pop: diploid population bi-allelic loci
    :return: numpy.array of allele frequencies
    """

#    allele_data_table = gather_allele_data(pop)

    allele_frequency_table = np.zeros((pop.totNumLoci(), 5))

    allele_frequency_table[:, 0] = allele_state_table[:, 0]
    for i in range(1, 5):
        allele_frequency_table[:, i] = [
            pop.dvars().alleleFreq[locus][allele] for
                locus, allele in zip(allele_state_table[:, 0],
                                     allele_state_table[:, i])]

    return allele_frequency_table


# todo Add documentation for gather_genotype_frequencies into analyze.rst
# todo Add usage example for gather_genotype_frequencies into analyze.rst
# todo Update documentation of gather_genotype_frequencies to reflect separation of steps

def gather_genotype_frequencies(pop):
    """
    Constructs a 3D array with dimensions number_of_loci x 5 x 5
    Genotypes are treated as 2D coordinates along the locus axis

    .. note:: To get the genotypes by locus we can use a single line of code.

    .. code-block:: python

        >>> genotypes_by_locus = np.array(np.ndarray.nonzero(genotype_frequency_array)).T

    :param sim.Population pop: Diploid population with bi-allelic loci
    :return: 3D genotype frequency array
    """

    genotype_frequency_array = np.ndarray((pop.totNumLoci(), 5, 5))

    for locus in range(pop.totNumLoci()):
        for genotype in pop.dvars().genoFreq[locus].keys():
            genotype_frequency_array[locus][genotype] = \
            pop.dvars().genoFreq[locus][genotype]

    return genotype_frequency_array


def tassel_results_tables(gwas_file_name, q_values_file_name,
                          minor_allele_frequency_table,
                          quantitative_allele_table):
    raw_gwas_results = pd.read_csv(gwas_file_name, sep='\t')
    raw_gwas_results.drop(0, axis=0, inplace=True)
    raw_gwas_results.drop('Trait', axis=1, inplace=True)
    raw_gwas_results.index = np.array(list(map(int, raw_gwas_results.Marker)))
    q_values = pd.read_csv(q_values_file_name, sep='\t')
    q_values.index = np.array(list(map(int, raw_gwas_results.Marker)))
    raw_gwas_results = raw_gwas_results.join(q_values)

    assert minor_allele_frequency_table.index.dtype == raw_gwas_results.index.dtype, "Indexes of these tables are different"

    raw_gwas_results = raw_gwas_results.join(
        minor_allele_frequency_table.ix[raw_gwas_results.index, :])

    assert quantitative_allele_table.index.dtype == raw_gwas_results.index.dtype, "Indexes of these tables are different"

    raw_gwas_results = raw_gwas_results.join(
        quantitative_allele_table.ix[raw_gwas_results.index, :])
    return raw_gwas_results


def generate_allele_effects_table(population_allele_frequencies, allele_array,
                                  allele_effect_array):
    """
    Creates a pandas DataFrame with the columns:
    + alpha allele
    + alpha allele effect
    + alpha allele frequency
    + beta allele
    + beta allele effect
    + beta allele frequency
    + difference

    :warning:`Assumes di-allelic case`
    """
    column_labels = ['alpha', 'alpha_effect', 'alpha_frequency',
                     'beta', 'beta_effect', 'beta_frequency', 'difference']
    number_of_loci = len(population_allele_frequencies)
    alpha_alleles, beta_alleles = allele_array[:, 0], allele_array[:, 1]
    alpha_effects, beta_effects = allele_effect_array[
                                      range(number_of_loci), alpha_alleles], \
                                  allele_effect_array[
                                      range(number_of_loci), beta_alleles]
    alpha_frequencies = np.asarray([population_allele_frequencies[locus][allele]
                                    for locus, allele in
                                    enumerate(alpha_alleles)])
    beta_frequencies = np.asarray([population_allele_frequencies[locus][allele]
                                   for locus, allele in
                                   enumerate(beta_alleles)])
    differences = np.abs(alpha_effects - beta_effects)
    allele_effects_table = pd.DataFrame(
        np.asarray([alpha_alleles, alpha_effects, alpha_frequencies,
                    beta_alleles, beta_effects, beta_frequencies,
                    differences]).T,
        columns=column_labels)

    return allele_effects_table


def minor_allele_frequencies_table(population_allele_frequencies, minor_alleles):
    """
    Creates a pd.DataFrame of minor allele frequencies with columns for the
    minor allele and its frequency.

    :param population_allele_frequencies:
    :param minor_alleles:
    :return:
    """
    column_labels = ['minor_allele', 'minor_frequency']
    minor_allele_frequencies = np.array([population_allele_frequencies[locus][allele]
                                for locus, allele in enumerate(minor_alleles)])
    return pd.DataFrame(np.array([minor_alleles, minor_allele_frequencies]).T,
                        columns=column_labels)


@singledispatch
def combine_population_samples(population_sample_list):
    """
    Given a list of separate simuPOP.Populations this function will add all
    the individuals to the first population of the list. Order of the
    individuals is preserved.
    :param population_sample_list: List of simuPOP.Populations
    :return:
    """
    for sample in population_sample_list[1:]:
        population_sample_list[0].addIndFrom(sample)
    return population_sample_list[0]

@combine_population_samples.register(list)
def list_of_populations(population_sample_list):
    """
    Given a list of separate simuPOP.Populations this function will add all
    the individuals to the first population of the list. Order of the
    individuals is preserved.
    :param population_sample_list: List of simuPOP.Populations
    :return:
    """
    for sample in population_sample_list[1:]:
        population_sample_list[0].addIndFrom(sample)
    return population_sample_list[0]


@combine_population_samples.register(dict)
def dictionary_of_populations(population_sample_library):
    for sample_list in population_sample_library.values():
        for sample in sample_list[1:]:
            sample_list[0].addIndFrom(sample)




class MetaGeneration(object):
    """
    A collection of functions to handle the special case of meta-populations.
    A meta-population is a population of individuals from different generations.
    The functions attached to this class will store individual generations
    then combine the individual samples into a meta-population. The output
    will be stored in an HDF5 file using the h5py package.
    """

    def __init__(self, hdf5):
        """
        hdf5 is a string of an HDF5 file name. Results of any of the methods
        are stored in this file. Methods for each type of data will be stored
        in sub-groups of the same HDF5 file.

        :param str hdf5: HDF5 file name
        """
        self.hdf5 = hdf5


    def extract_genotype_matrix(self, population):
        """
        Converts the genotype of the population into a numpy array for ease
        of use with other modules. Assumes that each locus has less than ~100
        alleles.

        :param population: simuPOP.Population
        """
        genotype_matrix = np.zeros((population.popSize(),
                                    population.genoSize()), dtype=np.int8)
        for row, ind in enumerate(population):
            genotype_matrix[row, ...] = ind.genotype()

        return genotype_matrix


    def store_allele_frequency_data(self, meta_population_library, minor_alleles,
                                    hdf_file_name):
        """
        Collects minor allele frequency data of a multiple generation
        population library. Stores the allele frequency data in an
        HDF5 file.

            af/replicate_id/generation_id

        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    gen_id = int(max(sample.indInfo('generation')))
                    hdf_file.create_dataset('af/' + str(rep_id) + '/' + str(gen_id),
                                            data=np.asarray(list(
                                                sample.dvars().alleleFreq[locus][
                                                    allele]
                                                for locus, allele in
                                                enumerate(minor_alleles))))


    def store_genotype_phenotype_data(self, meta_population_library, hdf_file_name):
        """
        Collects the genotype and phenotype data of a multiple replicate
        multiple sample population dictionary. Stores the results in
        an HDF5 file.

        Stored by
        geno_pheno/replicate_id/generation_id
        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    gen_id = int(max(sample.indInfo('generation')))
                    sample.setIndInfo(rep_id, 'replicate')
                    g_and_p = np.asarray((sample.indInfo('ind_id'),
                                          sample.indInfo('replicate'),
                                          sample.indInfo('generation'),
                                          sample.indInfo('g'),
                                          sample.indInfo('p'))).T
                    hdf_file.create_dataset(
                        'geno_pheno/' + str(rep_id) + '/' + str(gen_id),
                        data=g_and_p)


    def store_heterozygote_frequency_data(self, meta_population_library,
                                          hdf_file_name):
        """
        Collects minor allele frequency data of a multiple generation
        population library.Stores the allele frequency data in an
        HDF5 file.

        hetf replicate_id/generation_id
        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    gen_id = int(max(sample.indInfo('generation')))
                    hdf_file.create_dataset(
                        'hetf/' + str(rep_id) + '/' + str(gen_id),
                        data=np.asarray(list(sample.dvars().heteroFreq.values())))


    def store_genotype_frequency_data(self, meta_population_library,
                                      minor_alleles, hdf_file_name):
        """
        Collects the frequency of the minor allele homozygote data
        of a multiple replicate multiple sample population dictionary. The minor
        allele genotypes are created using the ``minor_alleles`` parameter.
        Stores the results in an HDF5 file.

        Keyed by

            homf/replicate_id/generation_id
        """
        minor_homozygote_genotypes = tuple(zip(minor_alleles, minor_alleles))
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    gen_id = int(max(sample.indInfo('generation')))
                    hdf_file.create_dataset(
                        'homf/' + str(rep_id) + '/' + str(gen_id), data=np.asarray(
                            tuple(sample.dvars().genoFreq[locus][genotype]
                                  for locus, genotype in
                                  enumerate(minor_homozygote_genotypes))))


    def multiple_sample_analyzer(self, meta_population_library, qtl, allele_effects,
                                 minor_alleles, loci):
        int_to_snp_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '-', 5: '+'}
        indir = "/home/vakanas/tassel-5-standalone/"
        minor_allele_frequency_file = h5py.File(hdf_file_name)
        run_id = self.run_id

        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                gen_id_name = str(int(max(sample.indInfo('generation'))))
                rep_id_name = str(rep_id)

                operators.assign_additive_g(sample, qtl,
                                            allele_effects)
                operators.calculate_error_variance(sample, 0.7)
                operators.calculate_p(sample)

                name = run_id + '_' + rep_id_name + '_' + gen_id_name

                af_access_name = '/'.join(['af', rep_id_name, gen_id_name])
                minor_allele_frequencies = \
                minor_allele_frequency_file[af_access_name][list(loci)]

                gwas = GWAS(sample, loci, run_id)

                ccm = gwas.calculate_count_matrix(minor_alleles, loci)
                ps_svd = gwas.pop_struct_svd(ccm)
                gwas.population_structure_formatter(ps_svd,
                                                    indir + name + '_structure_matrix.txt')
                gwas.hapmap_formatter(int_to_snp_map,
                                      indir + name + '_simulated_hapmap.txt')
                gwas.trait_formatter(indir + name + '_phenotype_vector.txt')
                gwas.calc_kinship_matrix(ccm, minor_allele_frequencies,
                                         indir + name + '_kinship_matrix.txt')

                gwas.replicate_tassel_gwas_configs(rep_id_name,
                                                   gen_id_name,
                                                   indir + name + '_simulated_hapmap.txt',
                                                   indir + name + '_kinship_matrix.txt',
                                                   indir + name + '_phenotype_vector.txt',
                                                   indir + name + '_structure_matrix.txt',
                                                   "C:\\tassel\\output\\" + name + '_out_',
                                                   "C:\\Users\DoubleDanks\\BISB\\wisser\\code\\rjwlab-scripts\\"
                                                   "saegus_project\\devel\\magic\\1478\\gwas_pipeline.xml")

        minor_allele_frequency_file.close()





class MultiGeneration(object):

    def __init__(self, run_id):
        self.run_id = run_id

    def collect_allele_frequency_data(self, meta_population_library, minor_alleles):
        """
        Collects minor allele frequency data of a multiple generation
        population library. Columns of the resulting array are:

        + replicate
        + generation
        + locus1
        + locus2
        + so on and so forth
        """
        minor_allele_frequencies = []
        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                gen_id = int(max(sample.indInfo('generation')))
                minor_allele_frequencies.append(np.asarray(
                    ([rep_id, gen_id] + list(sample.dvars().alleleFreq[locus][allele]
                                          for locus, allele in
                                          enumerate(minor_alleles)))))
        return np.asarray(minor_allele_frequencies)

    def collect_genotype_phenotype_data(self, meta_population_library):
        """
        Collects the genotype and phenotype data of a multiple replicate
        multiple sample population dictionary. The resulting data is
        a single array. Each row has ind_id, replicate, generation, g and p.
        """
        genotype_phenotype_data = []
        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                sample.setIndInfo(rep_id, 'replicate')
                genotype_phenotype_data.extend(np.asarray((sample.indInfo('ind_id'),
                                    sample.indInfo('replicate'),
                                    sample.indInfo('generation'),
                                    sample.indInfo('g'),
                                    sample.indInfo('p'))).T)
        return np.asarray(genotype_phenotype_data)

    def collect_heterozygote_frequency_data(self, meta_population_library):
        """
        Collects minor allele frequency data of a multiple generation
        population library.
        """
        heterozygote_frequencies = []
        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                gen_id = int(max(sample.indInfo('generation')))
                heterozygote_frequencies.append(np.asarray(
                    ([rep_id, gen_id] + list(sample.dvars().heteroFreq.values()))))
        return np.asarray(heterozygote_frequencies)

    def collect_genotype_frequency_data(self, meta_population_library, minor_alleles):
        """
        Collects the frequency of the minor allele homozygote data
        of a multiple replicate multiple sample population dictionary. The minor
        allele homozygote genotypes are constructed from the ``minor_alleles``
        parameter. Results are returned as a single array with columns

        + replicate
        + generation
        + locus1
        + locus2
        + so on and so forth

        """

        minor_allele_homozygotes = tuple(zip(minor_alleles, minor_alleles))
        minor_allele_homozygote_frequencies = []
        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                gen_id = int(max(sample.indInfo('generation')))
                minor_allele_homozygote_frequencies.append(np.asarray(([rep_id,
                                                                        gen_id] + list(
                    sample.dvars().genoFreq[locus][genotype]
                    for locus, genotype in enumerate(minor_allele_homozygotes)))))
        return np.asarray(minor_allele_homozygote_frequencies)

    def store_allele_frequency_data(self, meta_population_library, minor_alleles,
                                    hdf_file_name):
        """
        Collects minor allele frequency data of a multiple generation
        population library. Stores the allele frequency data in an
        HDF5 file.

            af/replicate_id/generation_id

        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    gen_id = int(max(sample.indInfo('generation')))
                    hdf_file.create_dataset('af/' + str(rep_id) + '/' + str(gen_id),
                        data=np.asarray(list(sample.dvars().alleleFreq[locus][allele]
                             for locus, allele in enumerate(minor_alleles))))


    def store_genotype_phenotype_data(self, meta_population_library, hdf_file_name):
        """
        Collects the genotype and phenotype data of a multiple replicate
        multiple sample population dictionary. Stores the results in
        an HDF5 file.

        Stored by
        geno_pheno/replicate_id/generation_id
        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    gen_id = int(max(sample.indInfo('generation')))
                    sample.setIndInfo(rep_id, 'replicate')
                    g_and_p = np.asarray((sample.indInfo('ind_id'),
                                          sample.indInfo('replicate'),
                                          sample.indInfo('generation'),
                                          sample.indInfo('g'),
                                          sample.indInfo('p'))).T
                    hdf_file.create_dataset('geno_pheno/'+str(rep_id)+'/'+str(gen_id),
                                            data=g_and_p)

    def store_heterozygote_frequency_data(self, meta_population_library, hdf_file_name):
        """
        Collects minor allele frequency data of a multiple generation
        population library.Stores the allele frequency data in an
        HDF5 file.

        hetf replicate_id/generation_id
        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    gen_id = int(max(sample.indInfo('generation')))
                    hdf_file.create_dataset(
                        'hetf/'+str(rep_id)+'/'+str(gen_id),
                        data=np.asarray(list(sample.dvars().heteroFreq.values())))

    def store_genotype_frequency_data(self, meta_population_library,
                                      minor_alleles, hdf_file_name):
        """
        Collects the frequency of the minor allele homozygote data
        of a multiple replicate multiple sample population dictionary. The minor
        allele genotypes are created using the ``minor_alleles`` parameter.
        Stores the results in an HDF5 file.

        Keyed by

            homf/replicate_id/generation_id
        """
        minor_homozygote_genotypes = tuple(zip(minor_alleles, minor_alleles))
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    gen_id = int(max(sample.indInfo('generation')))
                    hdf_file.create_dataset('homf/'+str(rep_id)+'/'+str(gen_id), data=np.asarray(
                    tuple(sample.dvars().genoFreq[locus][genotype]
                         for locus, genotype in enumerate(minor_homozygote_genotypes))))

    def multiple_sample_analyzer(self, meta_population_library, qtl, allele_effects,
                                 minor_alleles, loci):


        int_to_snp_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '-', 5: '+'}
        indir = "/home/vakanas/tassel-5-standalone/"
        minor_allele_frequency_file = h5py.File(hdf_file_name)
        run_id = self.run_id

        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:

                gen_id_name = str(int(max(sample.indInfo('generation'))))
                rep_id_name = str(rep_id)

                operators.assign_additive_g(sample, qtl,
                                    allele_effects)
                operators.calculate_error_variance(sample, 0.7)
                operators.calculate_p(sample)

                name = run_id+'_'+rep_id_name+'_'+gen_id_name

                af_access_name = '/'.join(['af', rep_id_name, gen_id_name])
                minor_allele_frequencies = minor_allele_frequency_file[af_access_name][list(loci)]

                gwas = GWAS(sample, loci, run_id)

                ccm = gwas.calculate_count_matrix(minor_alleles, loci)
                ps_svd = gwas.pop_struct_svd(ccm)
                gwas.population_structure_formatter(ps_svd,
                                        indir + name + '_structure_matrix.txt')
                gwas.hapmap_formatter(int_to_snp_map,
                                      indir + name + '_simulated_hapmap.txt')
                gwas.trait_formatter(indir + name + '_phenotype_vector.txt')
                gwas.calc_kinship_matrix(ccm, minor_allele_frequencies,
                                 indir + name + '_kinship_matrix.txt')

                gwas.replicate_tassel_gwas_configs(rep_id_name,
                                                   gen_id_name,
                                       indir + name + '_simulated_hapmap.txt',
                                       indir + name + '_kinship_matrix.txt',
                                       indir + name + '_phenotype_vector.txt',
                                       indir + name + '_structure_matrix.txt',
                                       "C:\\tassel\\output\\" + name + '_out_',
                                                   "C:\\Users\DoubleDanks\\BISB\\wisser\\code\\rjwlab-scripts\\"
                                                   "saegus_project\\devel\\magic\\1478\\gwas_pipeline.xml")

        minor_allele_frequency_file.close()


    def calculate_additive_genetic_variance(self, qtl, founder_allele_data,
                                            minor_allele_frequency_data,
                                            allele_effect_array):
        """
        Calculates additive genetic variance according to the formula:

            V_a = p * (1 - p) * (a^2)
            where a is the average effect (major + minor)/2

        :param qtl: numpy array of integers of qtl
        :param founder_allele_data: pandas DataFrame with columns \
                                        minor_allele, major_allele
        :param minor_allele_frequency_data: numpy array with columns rep, gen, \
                                                    locus1, locus2, ...
        :param allele_effect_array: Array of allele effects
        """
        quant_minor_alleles = np.asarray(founder_allele_data.minor_allele)[qtl]
        quant_major_alleles = np.asarray(founder_allele_data.major_allele)[qtl]
        minor_allele_effects = allele_effect_array[qtl, quant_minor_alleles]
        major_allele_effects = allele_effect_array[qtl, quant_major_alleles]
        quant_minor_allele_frequencies = minor_allele_frequency_data[:, qtl+2]
        reps = max(minor_allele_frequency_data[:, 0]) + 1
        gens = (max(minor_allele_frequency_data[:, 1]))/2 + 1


        average_effects = (1/2)*(minor_allele_effects + major_allele_effects)
        additive_variance = np.zeros((int(reps*gens), qtl.shape[0]))
        additive_variance[:, ...] = quant_minor_allele_frequencies*\
                                    (1 - quant_minor_allele_frequencies)\
                                    *(average_effects**2)
        additive_variance = np.hstack((minor_allele_frequency_data[:, :2],
                                       additive_variance))
        additive_variance = pd.DataFrame(additive_variance,
                                 columns=['rep', 'gen']+list(map(str, qtl)))
        additive_variance.rep = np.asarray(additive_variance.rep, dtype=np.int64)
        additive_variance.gen = np.asarray(additive_variance.gen, dtype=np.int64)

        return additive_variance


class SingleGeneration(object):
    def __init__(self):
        pass

    def collect_allele_frequency_data(self, meta_population_library, minor_alleles):
        """
        Collects minor allele frequency data of a multiple generation
        population library. Columns of the resulting array are:

        + replicate
        + generation
        + locus1
        + locus2
        + so on and so forth
        """
        minor_allele_frequencies = []
        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                #gen_id = int(max(sample.indInfo('generation')))
                minor_allele_frequencies.append(np.asarray(
                    ([rep_id, sample.popSize()] + list(
                        sample.dvars().alleleFreq[locus][allele]
                        for locus, allele in
                        enumerate(minor_alleles)))))
        return np.asarray(minor_allele_frequencies)

    def collect_genotype_phenotype_data(self, meta_population_library):
        """
        Collects the genotype and phenotype data of a multiple replicate
        multiple sample population dictionary. The resulting data is
        a single array. Each row has ind_id, replicate, generation, g and p.
        """
        genotype_phenotype_data = []
        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                size_identifier = np.array([sample.popSize()]*sample.popSize())
                sample.setIndInfo(rep_id, 'replicate')
                genotype_phenotype_data.extend(
                    np.asarray((sample.indInfo('ind_id'),
                                size_identifier,
                                sample.indInfo('generation'),
                                sample.indInfo('g'),
                                sample.indInfo('p'))).T)
        return np.asarray(genotype_phenotype_data)

    def collect_heterozygote_frequency_data(self, meta_population_library):
        """
        Collects minor allele frequency data of a multiple generation
        population library.
        """
        heterozygote_frequencies = []
        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                heterozygote_frequencies.append(np.asarray(
                    ([rep_id, sample.popSize()] + list(
                        sample.dvars().heteroFreq.values()))))
        return np.asarray(heterozygote_frequencies)

    def collect_genotype_frequency_data(self, meta_population_library, minor_alleles):
        """
        Collects the frequency of the minor allele homozygote data
        of a multiple replicate multiple sample population dictionary. The minor
        allele homozygote genotypes are constructed from the ``minor_alleles``
        parameter. Results are returned as a single array with columns

        + replicate
        + generation
        + locus1
        + locus2
        + so on and so forth

        """

        minor_allele_homozygotes = tuple(zip(minor_alleles, minor_alleles))
        minor_allele_homozygote_frequencies = []
        for rep_id, sample_list in meta_population_library.items():
            for sample in sample_list:
                minor_allele_homozygote_frequencies.append(
                    np.asarray(([rep_id, sample.popSize()] + list(
                    sample.dvars().genoFreq[locus][genotype]
                    for locus, genotype in
                    enumerate(minor_allele_homozygotes)))))
        return np.asarray(minor_allele_homozygote_frequencies)

    def store_allele_frequency_data(self, meta_population_library, minor_alleles,
                                    hdf_file_name):
        """
        Collects minor allele frequency data of a multiple generation
        population library. Stores the allele frequency data in an
        HDF5 file.

            af/replicate_id/generation_id

        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    hdf_file.create_dataset(
                        'af/' + str(rep_id) + '/' + str(sample.popSize()),
                        data=np.asarray(
                            list(sample.dvars().alleleFreq[locus][allele]
                                 for locus, allele in
                                 enumerate(minor_alleles))))

    def store_genotype_phenotype_data(self, meta_population_library, hdf_file_name):
        """
        Collects the genotype and phenotype data of a multiple replicate
        multiple sample population dictionary. Stores the results in
        an HDF5 file.

        Stored by
        geno_pheno/replicate_id/generation_id
        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    size_identifier = np.asarray([sample.popSize()]*sample.popSize())
                    sample.setIndInfo(rep_id, 'replicate')
                    g_and_p = np.asarray((sample.indInfo('ind_id'),
                                          sample.indInfo('replicate'),
                                          size_identifier,
                                          sample.indInfo('g'),
                                          sample.indInfo('p'))).T
                    hdf_file.create_dataset(
                        'geno_pheno/' + str(rep_id) + '/' + str(gen_id),
                        data=g_and_p)

    def store_heterozygote_frequency_data(self, meta_population_library,
                                          hdf_file_name):
        """
        Collects minor allele frequency data of a multiple generation
        population library.Stores the allele frequency data in an
        HDF5 file.

        hetf replicate_id/generation_id
        """
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    hdf_file.create_dataset(
                        'hetf/' + str(rep_id) + '/' + str(sample.popSize()),
                        data=np.asarray(
                            list(sample.dvars().heteroFreq.values())))

    def store_genotype_frequency_data(self, meta_population_library, minor_alleles,
                                      hdf_file_name):
        """
        Collects the frequency of the minor allele homozygote data
        of a multiple replicate multiple sample population dictionary. The minor
        allele genotypes are created using the ``minor_alleles`` parameter.
        Stores the results in an HDF5 file.

        Keyed by

            homf/replicate_id/generation_id
        """
        minor_homozygote_genotypes = tuple(zip(minor_alleles, minor_alleles))
        with h5py.File(hdf_file_name) as hdf_file:
            for rep_id, sample_list in meta_population_library.items():
                for sample in sample_list:
                    hdf_file.create_dataset(
                        'homf/' + str(rep_id) + '/' + str(sample.popSize()),
                        data=np.asarray(
                            tuple(sample.dvars().genoFreq[locus][genotype]
                                  for locus, genotype in
                                  enumerate(minor_homozygote_genotypes))))


def allele_data(pop, alleles, loci):
    """
    Determines the minor alleles, minor allele frequencies, major alleles and
    major allele frequencies.

    Example
    ~~~~~~~

        pop = sim.loadPopulation('magic1478.pop')
        loci = list(range(10, 200, 5))
        alleles = shelve.open('magic_1478_simulation_parameters')


    :param pop:
    :type pop:
    :param loci:
    :type loci:
    :return:
    :rtype:
    """

    reversed_allele_frequencies = {}
    allele_frq = {}

    for locus in loci:
        reversed_allele_frequencies[locus] = {}
        for allele in alleles[locus]:
            frequency = pop.dvars().alleleFreq[locus][allele]
            reversed_allele_frequencies[locus][frequency] = allele

    allele_frq['minor', 'alleles'] = col.OrderedDict()
    allele_frq['minor', 'frequency'] = col.OrderedDict()
    allele_frq['frequencies'] = col.OrderedDict()
    allele_frq['major', 'alleles'] = col.OrderedDict()
    allele_frq['major', 'frequency'] = col.OrderedDict()

    # Determine major/minor allele in aggregate population
    for locus in loci:
        temp_frq = []
        for allele in alleles[locus]:
#            np.array([pop.dvars().alleleFreq[locus][allele]])
            temp_frq.append(pop.dvars().alleleFreq[locus][allele])

        allele_frq['frequencies'][locus] = temp_frq
        minor_frequency = min(allele_frq['frequencies'][locus])
        major_frequency = max(allele_frq['frequencies'][locus])
        allele_frq['minor', 'alleles'][locus] = reversed_allele_frequencies[locus][minor_frequency]
        allele_frq["major", "alleles"][locus] = reversed_allele_frequencies[locus][major_frequency]

    # In some cases alleles may be at exactly equal frequencies i.e. both alleles
    # at 0.5. The result will be the inability to consistently distinguish
    # bewtween the major and minor allele.

    ties = np.array([locus for locus in loci
                    if allele_frq['minor', 'alleles'][locus] ==
                         allele_frq['major', 'alleles'][locus]])

    for tied_alleles in ties:
        allele_frq['major', 'alleles'][tied_alleles] = list(pop.dvars().alleleFreq[tied_alleles])[0]
        allele_frq['minor', 'alleles'][tied_alleles] = list(pop.dvars().alleleFreq[tied_alleles])[1]

    # Test to make sure the major/minor allele conflict has been resolved to
    # by making sure that the major and minor alleles are different at every
    # locus. Assert that the number of matches is equal to 0.

    major_minor_allele_conflicts = sum(
        np.equal(list(allele_frq['minor', 'alleles'].values()),
                 list(allele_frq['major', 'alleles'].values())))

    assert major_minor_allele_conflicts == 0, "At least one locus defines the major allele to be the same" \
                                              "as the minor allele.: {}".format(ties)

    minor_alleles = np.asarray([allele_frq['minor', 'alleles'][locus] for locus in loci], dtype=np.int8)
    minor_frequencies = [pop.dvars().alleleFreq[locus][minor_allele] for locus, minor_allele in enumerate(minor_alleles)]
    major_alleles = np.asarray([allele_frq['major', 'alleles'][locus] for locus in loci], dtype=np.int8)
    major_frequencies = [pop.dvars().alleleFreq[locus][major_allele] for locus, major_allele in enumerate(major_alleles)]

    allele_data_structure = \
        pd.DataFrame(np.array([minor_alleles, minor_frequencies, major_alleles, major_frequencies]),
                     index=['minor_allele', 'minor_frequency', 'major_allele', 'major_frequency'],
                     columns=loci).T

    return allele_data_structure


def allele_frq_table(pop, number_gens,
    allele_frq_data, recombination_rates, genetic_map):
    """
    Generates a large table which centralizes all allele frequency data.
    The data is inserted into a pandas DataFrame object.
    Useful for downstream analysis and insertion into a database.

    Allele frequency data is first built up in a regular *dict* object
    then inserted into a
    :param pop:
    :param number_gens:
    :param allele_frq_data:
    :param recombination_rates:
    :param genetic_map:
    """


    data_columns = ['abs_index', 'chrom', 'locus', 'major',
                    'minor', 'recom_rate', 'cM', 'v']
    generation_labels = ['G_'+str(i)
                         for i in range(0, number_gens+1, 2)]
    data_columns = data_columns + generation_labels + ['aggregate']
    data = {}
    breakpoints = col.OrderedDict()
    for locus in range(pop.totNumLoci()):
        if recombination_rates[locus] == 0.01:
            breakpoints[locus] = locus + 1

    diagram = ["|"]*pop.totNumLoci()

    for locus, point in breakpoints.items():
        try:
                diagram[point] = '*'
        except IndexError:
            pass

    qtl_diagram = ['o']*pop.totNumLoci()
    for locus in pop.dvars().triplet_qtl:
        qtl_diagram[locus] = 'x'

    chromosomes = []
    relative_loci = []
    for locus in range(pop.totNumLoci()):
        pair = pop.chromLocusPair(locus)
        chromosomes.append(pair[0]+1)
        relative_loci.append(pair[1])

    data['chrom'] = chromosomes
    data['locus'] = relative_loci


    data['recom_rate'] = recombination_rates
    data['v'] = diagram
    data['cM'] = genetic_map['cM_pos']
    data['qtl'] = qtl_diagram


    data['abs_index'] = [locus for locus in range(pop.totNumLoci())]
    data['minor'] = [allele_frq_data['minor', 'alleles'][locus] for locus in range(pop.totNumLoci())]
    data['major'] = [allele_frq_data['major', 'alleles'][locus] for locus in range(pop.totNumLoci())]
    for sp_idx, label in enumerate(generation_labels):
        data[label] = [allele_frq_data['minor', 'frequency', sp_idx][
                           locus] for locus in range(pop.totNumLoci())]
    data['aggregate'] = [allele_frq_data['minor', 'frequency'][locus] for locus in range(pop.totNumLoci())]
    af_table = pd.DataFrame(data, columns=data_columns)
    return af_table


def collect_haplotype_data(pop, allele_effects, quantitative_trait_loci):
    """

    :param pop:
    :type pop:
    :param allele_effects:
    :type allele_effects:
    :param quantitative_trait_loci:
    :type quantitative_trait_loci:
    :return:
    :rtype:
    """

    haplotypes = {}
    haplotypes['loci'] = {}
    for k, i in enumerate(range(0, len(quantitative_trait_loci), 3)):
        haplotypes['loci'][k] = (quantitative_trait_loci[i],
                             quantitative_trait_loci[i+1],
                             quantitative_trait_loci[i+2])

    haplotypes['alleles'] = {}
    haplotypes['effect'] = {}
    haplotypes['frequency'] = {}
    for loci in haplotypes['loci'].values():
        haplotypes['frequency'][loci] = {}
        for sp in range(pop.numSubPop()):
            haplotypes['frequency'][loci][sp] = {}

    sim.stat(pop, haploFreq=list(haplotypes['loci'].values()),
             vars=['haploFreq', 'haploFreq_sp'])

    for k, v in haplotypes['loci'].items():
        haplotypes['alleles'][v] = list(pop.dvars(0).haploFreq[v].keys())

    for sp in range(pop.numSubPop()):
        for loci, triplet in haplotypes['alleles'].items():
            for alleles in triplet:
                haplotypes['frequency'][loci][sp][alleles] = pop.dvars(
                    sp).haploFreq[loci][alleles]

    for htype, triplets in haplotypes['alleles'].items():
        haplotypes['effect'][htype] = {}
        for trip in triplets:
            htype_effect = allele_effects[htype[0]][trip[0]] +\
            allele_effects[htype[1]][trip[1]] +\
            allele_effects[htype[2]][trip[2]]
            haplotypes['effect'][htype][trip] = htype_effect

    return haplotypes


def generate_haplotype_data_table(pop, haplotype_data):
    """
    Generates a table for easy analysis and visualization of haplotypes,
    effects, frequencies and locations.


    :param pop:
    :type pop:
    :param haplotype_data:
    :type haplotype_data:
    :return:
    :rtype:
    """
    integer_to_snp = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '+', 5: '-'}
    haplotype_table = []
    data_columns = ['centered_on', 'relative_position', 'chromosome',
                    'haplotype', 'effect']
    generation_columns = ['G_'+str(i) for i in range(0, 2*(pop.numSubPop()),
                                                     2)]
    data_columns.extend(generation_columns)
    for locus in haplotype_data['loci'].values():
        for triplet in haplotype_data['alleles'][locus]:
            generational_frequencies = [haplotype_data['frequency'][locus][sp][triplet]
                                        for sp in range(pop.numSubPop())]
            effect = haplotype_data['effect'][locus][triplet]
            snp_triplet = integer_to_snp[triplet[0]] + \
                          integer_to_snp[triplet[1]] + \
                          integer_to_snp[triplet[2]]
            chromosome = pop.chromLocusPair(locus[1])[0] + 1
            relative_locus = pop.chromLocusPair(locus[1])[1]
            row = [locus[1]] + \
                  [relative_locus] +\
                  [chromosome] + \
                  [snp_triplet] + \
                  [effect] + \
                  generational_frequencies
            haplotype_table.append(row)
    return pd.DataFrame(haplotype_table, columns=data_columns)


class GWAS(object):
    """
    Class for performing principal component analyis on genotype matrices.
    Test for population structure significance tests the largest eigenvalue
    of the genotype covarience matrix. Details can be found in the paper:
    Population Structure and Eigenanalysis Patterson et al 2006.
    """
    def __init__(self, pop, loci, run_id):
        self.pop = pop
        self.loci = loci
        self.run_id = run_id
        self.int_to_snp_conversions = {0: 'A', 1: 'C',
                                       2: 'G', 3: 'T',
                                       4: '-', 5: '+'}


        self.individual_names = np.core.defchararray.add('I',
                    defchararray.array(np.asarray(self.pop.indInfo('ind_id'),
                        dtype=np.int_), copy=False, unicode=np.unicode_))

    def __str__(self):
        return "Conversions: \n "+ str(self.int_to_snp_conversions)


    # noinspection PyArgumentList
    def calculate_count_matrix(self, allele_subset, seg_loci,
                               count_matrix_file_name = None):
        """
        A function to calculate the copy numbers of either the minor or
        major allele for each individual at each locus.
        Minor or major
        alleles parameter is a single set of alleles which determines if the
        return is the minor or major allele count matrix.
        :param allele_subset: Allows user to choose a custom set of alleles to use i.e. minor vs major.
        :param count_matrix_filename: Output file name. If defined will write a file. Otherwise returns the count_matrix
        """
        comparison_array = np.asarray(allele_subset, dtype=np.int_)
        count_matrix = np.zeros((self.pop.popSize(), len(self.loci)), dtype=np.int_)
        for i, ind in enumerate(self.pop.individuals()):
            alpha_genotype = np.asarray(ind.genotype(ploidy=0))[list(seg_loci)]
            alpha_comparisons = np.equal(comparison_array, alpha_genotype,
                                         dtype=np.int8)
            omega_genotype = np.asarray(ind.genotype(ploidy=1))[list(seg_loci)]
            omega_comparisons = np.equal(comparison_array, omega_genotype,
                                         dtype=np.int8)
            counts = np.add(alpha_comparisons, omega_comparisons, dtype=np.int8)
            count_matrix[i, :] = counts

        if count_matrix_file_name is not None:
            np.savetxt(count_matrix_file_name, count_matrix, fmt="%d")

        return count_matrix


    def pop_struct_svd(self, count_matrix):
        """

        Follows procedure of Population Structure and Eigenanalysis
        Patterson et al 2006.
        Constructs a genotype matrix of bi-allelic loci where each entry is
        the number of copies of the major allele at each locus. The genotype
        matrix has dimensions (number_of_individuals)*(number_of_markers).
        :param pop:
        :param count_matrix:

        """
        shift = np.apply_along_axis(np.mean, axis=1, arr=count_matrix)
        p_vector = np.divide(shift, 2)
        scale = np.sqrt(np.multiply(p_vector, (1-p_vector)))

        shift_matrix = np.zeros((count_matrix.shape[0], count_matrix.shape[1]))
        scale_matrix = np.zeros((count_matrix.shape[0], count_matrix.shape[1]))
        for i in range(len(self.loci)):
            shift_matrix[:, i] = shift
            scale_matrix[:, i] = scale

        corrected_matrix = (count_matrix - shift_matrix)/scale_matrix
        # singular value decomposition using scipy linalg module
        eigenvectors, s, v = linalg.svd(corrected_matrix)
        eigenvalues = np.diagonal(
            np.square(
                linalg.diagsvd(s, count_matrix.shape[1], count_matrix.shape[0]))).T

        sum_of_eigenvalues = np.sum(eigenvalues)
        fraction_of_variance = np.divide(eigenvalues, sum_of_eigenvalues)
        eigen_data = {'vectors': eigenvectors, 'values': eigenvalues,
                      'fraction_variance': fraction_of_variance}
        return eigen_data

    def population_structure_formatter(self, eigen_data, pop_struct_file_name = None):
        """
        Writes the first two of the population structure matrix to a
        file. First column of the file is are names.
        :param structure_filename:
        :param eigen_data:
        """

        ordered_names = self.individual_names

        structure_matrix = pd.DataFrame([list(eigen_data['vectors'][:, 0].T),
                                         list(eigen_data['vectors'][:, 1].T)]).T
        structure_matrix.index = ordered_names

        if pop_struct_file_name is not None:
            header = "<Covariate>\t\t\n<Trait>\td1\td2\n"
            with open(pop_struct_file_name, 'w') as f:
                f.write(header)
                structure_matrix.to_csv(f, sep='\t', index=True, header=False)
#       structure_matrix.index = ordered_names
        return structure_matrix

    def test_statistic(self, eigenvalues):
        sum_of_eigenvalues = np.sum(eigenvalues)
        n_hat_numerator = (self.pop.popSize() + 1)*sum_of_eigenvalues
        n_hat_denom = (self.pop.popSize()-1)*sum_of_eigenvalues - sum_of_eigenvalues
        n_hat = n_hat_numerator/n_hat_denom
        lowercase_l = (self.pop.popSize() - 1)*eigenvalues[0]
        mu_hat = np.square((np.sqrt(n_hat - 1) +
                            np.sqrt(self.pop.popSize()))) / n_hat
        sigma_hat = ((np.sqrt(n_hat - 1) + np.sqrt(self.pop.popSize()))/n_hat) * \
                    (((1/np.sqrt(n_hat - 1)) + 1/np.sqrt(self.pop.popSize())) ** (
                        1 / 3.0))
        test_statistic = (lowercase_l - mu_hat) / sigma_hat
        return test_statistic

    def hapmap_formatter(self, segregating_loci,
                         alleles_column,
                             locus_names,
                             corresponding_chromosomes,
                             pos_column,
                             hapmap_file_name=None):

        # Need to guarantee that the column names are in same order as the
        # genoype data. Iterate over individuals in population to build up a
        #  list of names will guarantee that col names are in same order as
        # the hapmap_data


        hapmap_ordered_columns = ['rs', 'alleles',
                      'chrom', 'pos', 'strand', 'assembly',
                      'center', 'protLSID', 'assayLSID', 'panelLSID',
                                  'QCode'] + list(
            self.individual_names)

        segregating_loci_array = np.array(segregating_loci)
        hapmap_matrix = pd.DataFrame(columns=hapmap_ordered_columns)
        hapmap_matrix.rs = locus_names
        hapmap_matrix.alleles = alleles_column
        hapmap_matrix.chrom = corresponding_chromosomes
        hapmap_matrix.pos = pos_column
        hapmap_matrix.ix[:, 'strand':'QCode'] = np.core.defchararray.array([['NA']*len(locus_names)]*7).T

        for i, ind in enumerate(self.pop.individuals()):
            hapmap_matrix.ix[:, self.individual_names[i]] = [
                ''.join(sorted(self.int_to_snp_conversions[a]
                               + self.int_to_snp_conversions[b]))
                for a, b in zip(np.asarray(ind.genotype(ploidy=0))[list(segregating_loci)], np.asarray(ind.genotype(ploidy=1))[list(segregating_loci)])]

        if hapmap_file_name is not None:
            hapmap_matrix.to_csv(hapmap_file_name, sep='\t', index=False)

        return hapmap_matrix

    def trait_formatter(self, trait_file_name = None):
        """
        Simple function to automate the formatting of the phenotypic data.
        Because of the way the header must be written the file is opened in
        append mode. Rewriting to the same file many times could introduce an
        unsuspected bug.
        :param trait_filename:
        """

        trait_vector = pd.DataFrame(np.array([self.individual_names,
                               self.pop.indInfo('p')]).T)

        if trait_file_name is not None:
            header = "<Trait>\tsim\n"
            with open(trait_file_name, 'w') as f:
                f.write(header)
                trait_vector.to_csv(f, sep='\t', index=False, header=False)

        return trait_vector

    def calc_kinship_matrix(self, allele_count_matrix,
                            allele_frequencies,
                            kinship_matrix_file_name = None):
        """
        Calculates the kinship matrix according to VanRaden 2008:
        Efficient Methods to Compute Genomic Predictions and writes it to a
        file formatted for Tassel. The variable names try to be consistent
        with the variable names in the paper.

        The allele frequencies used for this function are with respect to
        the aggregate population: all individuals sampled during selection.

        Unsure if I should use the G0 allele or not.
        Will decide after TASSEL results.

        :param allele_count_matrix:
        :param allele_data:
        :param kinship_filename:
        :return:
        :rtype:
        """
        M = np.matrix(allele_count_matrix - 1)
        P = 2*(allele_frequencies - 0.5)
        Z = M - P
        scaling_terms = np.zeros((len(self.loci)))
        for idx, probability in enumerate(allele_frequencies):
            scaling_terms[idx] = 2*probability*(1 - probability)

        scaling_factor = sum(scaling_terms)

        G = Z*Z.T/scaling_factor

        annotated_G = pd.DataFrame(G, index=self.individual_names)

        if kinship_matrix_file_name is not None:
            header = "{}\n".format(self.pop.popSize())
            with open(kinship_matrix_file_name, 'w') as f:
                f.write(header)
                annotated_G.to_csv(f, sep='\t', index=True, header=False,
                                      float_format='%.3f')

        return annotated_G

    def replicate_tassel_gwas_configs(self, rep_id, sample_size,
                                          hapmap_file_name,
                                          kinship_file_name,
                                          phenotype_file_name,
                                          structure_file_name,
                                          output_file_prefix,
                                          config_file_template):


        """
        Creates an xml file to run TASSEL using a mixed linear model approach.
        Assumes use of hapmap, kinship, phenotype and population structure files.

        The TASSEL command line interface requires a considerable number of
        options to run GWAS. It is impractical to run the command line manually
        for the number of replications in a simulated study. The TASSEL command
        line interface allows the user to input a .xml file with the same
        information which is used in the terminal.

        :param input_directory: Directory path to send the input files.
        :param run_identifier_prefix: Identifier for single replicate of data
        :param config_file_templae: XML file already setup for running a
        specific kind of GWAS
        :return: XML file to run a single replicate of data using TASSEL
        """

        import xml.etree.ElementTree as ET
        import lxml.etree as etree

        tree = ET.parse(config_file_template)
        root = tree.getroot()
        lxml_tree = etree.fromstring(ET.tostring(root))
        lxml_root = lxml_tree.getroottree()

        lxml_root.find('fork1/h').text = hapmap_file_name
        lxml_root.find('fork2/t').text = phenotype_file_name
        lxml_root.find('fork3/q').text = structure_file_name
        lxml_root.find('fork4/k').text = kinship_file_name

        lxml_root.find('combine6/export').text = output_file_prefix

        lxml_root.write("C:\\tassel\\bin\\" + 'R' + rep_id + '_' + str(
            sample_size) + '_' + self.run_id + '_' + "_sim_gwas_pipeline.xml",
                        encoding="UTF-8",
                        method="xml", xml_declaration=True, standalone='',
                        pretty_print=True)


class GWASConfig(object):
    """
    Class to contain all the parameters associated with writing a "configFile"
    for TASSEL 5.0. A ``configFile`` is an xml document which instructs TASSEL
    what to do with a set of data.
    """

    def __init__(self, run_id, tassel_executable_path, tassel_input_path,
                 tassel_output_path, config_file_template, config_file_output):
        """

        This class is designed for the purpose of enabling platform independence.
        The user initalizes a GWASConfigFile object with the locations of where
        tassel is and where the user stores the input and output data.

        :param str tassel_executable_path: Location of tassel executable ``run_pipeline*``
        :param str tassel_input_path: Location of where the TASSEL input data is stored
        :param str tassel_outnput_path: Location of where TASSEL should store the output
        """
        self._run_id = run_id
        self._tassel_executable_path = tassel_executable_path
        self._tassel_output_path = tassel_output_path
        self._tassel_input_path = tassel_input_path
        self._config_file_template = config_file_template
        self._config_file_output_path = config_file_output_path

    def __str__(self):
        return "GWASConfig: {instance_name}\n" \
               "Location of TASSEL Executable: {tassel_executable}\n" \
               "Location of input: {tassel_input}\n" \
               "Destination of output: {tassel_output}".format(instance_name=self.__name__,
                         tassel_executable=_tassel_executable_path,
                         tassel_input=self._tassel_input_path,
                         tassel_output=self._tassel_output_path)



    @property
    def run_id(self):
        return self._run_id

    @property
    def tassel_executable_path(self):
        return self._tassel_executable_path

    @property
    def tassel_input_path(self):
        return self._tassel_input_path

    @property
    def tassel_output_path(self):
        return self._tassel_outnput_path

    @property
    def config_file_template(self):
        return self._config_file_template

    @property
    def config_file(self):
        return self.config_file

    def generate_config(self, sample_size, hapmap_file_name, kinship_file_name,
                            phenotype_file_name, structure_file_name,
                            output_file_prefix):
        """
        Creates an xml file to run TASSEL. The protocol performed by TASSEL
        is controlled by one of the xml tags. For example <mlm/> would run a
        moxed linear model using your data.


        The TASSEL command line interface requires a considerable number of
        options to run GWAS. It is impractical to run the command line manually
        for the number of replications in a simulated study. The TASSEL command
        line interface allows the user to input a .xml file with the same
        information which is used in the terminal.

        :param run_identifier_prefix: Identifier for single replicate of data
        :param config_file_templae: XML file already setup for running a
        specific kind of GWAS
        :return: XML file to run a single replicate of data using TASSEL
        """

        hapmap_file_path = os.join.path(self.tasel_input_path, hapmap_file_name)
        phenotype_file_path = os.join.path(self.tassel_input_path, phenotype_file_name)
        structure_file_path = os.join.path(self.tassel_input_path, structure_file_name)
        kinship_file_path = os.join.path(self.tassel_input_path, kinship_file_name)
        output_prefix = os.join.path(self.tassel_output_path, output_file_prefix)

        tree = ET.parse(self._config_file_template)
        root = tree.getroot()
        lxml_tree = etree.fromstring(ET.tostring(root))
        lxml_root = lxml_tree.getroottree()

        lxml_root.find('fork1/h').text = hapmap_file_path
        lxml_root.find('fork2/t').text = phenotype_file_path
        lxml_root.find('fork3/q').text = structure_file_path
        lxml_root.find('fork4/k').text = kinship_file_path

        lxml_root.find('combine6/export').text = output_file_prefix


        config_file_output_name = os.join(self._config_file_output_path, run_id, '_', str(sample_size),
                '_simulated_gwas_pipeline.xml')

        lxml_root.write(config_file_output_name,
                        encoding="UTF-8",
                       method="xml", xml_declaration=True, standalone='',
                        pretty_print=True)





def single_sample_analyzer(full_population, sample_size,
                               quantitative_trait_loci, alleles,
                                       allele_effects, heritability,
                                       run_id='infinite'):

    """
    A function to call all of the functions in order to get from population
    to TASSEL input and output.

    For the time being it is hard-coded in many aspects; however, this one
    function replaces an entire IPython notebook's worth of code.

    In future development this function will probably be separated into
    two or three separate methods under the same class.

    :warning: The paths for the location of the TasselIO are hard-coded for the time being.
    :warning: The path for the int_to_snp_map is also hard-coded.
    """

    syn_parameters = shelve.open('synthesis_parameters')
    int_to_snp_map = syn_parameters['integer_to_snp']
    syn_parameters.close()

    from . import parameters

    sample_population = sim.sampling.drawRandomSample(full_population,
                                                      sizes=sample_size)
    sim.stat(sample_population, alleleFreq=sim.ALL_AVAIL)
    sim.stat(sample_population, numOfSegSites=sim.ALL_AVAIL,
             vars=['segSites', 'numOfSegSites'])
    segregating_loci = sample_population.dvars().segSites
    operators.assign_additive_g(full_population, quantitative_trait_loci, allele_effects)
    operators.calculate_error_variance(sample_population, heritability)
    operators.calculate_p(sample_population)
    af = allele_data(sample_population, alleles,
                             range(sample_population.totNumLoci()))

    af.to_hdf(str(sample_size) + '_daoko_girl_af.hdf', 'af')


    gwas = GWAS(sample_population, segregating_loci,
                        np.array(af['minor_allele']), run_id)

    indir = "C:\\tassel\\input\\"
    ccm = gwas.calculate_count_matrix(indir + 'daoko_girl_MAC.txt')
    ps_svd = gwas.pop_struct_svd(ccm)
    ps_m = gwas.population_structure_formatter(ps_svd, indir + str(sample_size) +
                                               '_daoko_girl_structure_matrix.txt')
    hmap = gwas.hapmap_formatter(int_to_snp_map,
                                 indir + str(sample_size) + '_daoko_girl_simulated_hapmap.txt')
    phenos = gwas.trait_formatter(indir + str(sample_size) + '_daoko_girl_phenotype_vector.txt')
    ks_m = gwas.calc_kinship_matrix(ccm, af,
                                    indir + str(sample_size) + '_daoko_girl_kinship_matrix.txt')

    gwas.generate_tassel_gwas_configs(sample_size,
                  indir + str(sample_size) + '_daoko_girl_simulated_hapmap.txt',
                  indir + str(sample_size) + '_daoko_girl_kinship_matrix.txt',
                  indir + str(sample_size) + '_daoko_girl_phenotype_vector.txt',
                  indir + str(sample_size) + '_daoko_girl_structure_matrix.txt',
                  "C:\\tassel\\output\\" + str(sample_size) + "_daoko_girl_out_",
                  "C:\\Users\DoubleDanks\\BISB\\wisser\\code\\rjwlab-scripts\\"
              "saegus_project\\devel\\magic\\1478\\daoko_girl_gwas_pipeline.xml")
    return segregating_loci, aes_table

class Study(object):
    """
    Encapsulation of functions and properties to track information about a line
    of analysis or idea. Required parameter is a str run_id. Hence everything
    with the same run_id is related.

    """


    def __init__(self, run_id):
        """
        Studies are identified by their run_id. A study has its own set of
        input/output associated with it identifiable by the run_id in the file
        name.
        :param str run_id: Unique identifier to naturally group related pieces of data.
        """
        self.run_id = run_id
        self.integer_to_snp = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'D', 5: 'I'}
        self.snp_to_integer = {'A': 0, 'C': 1, 'D': 4, 'G': 2, 'I': 5, 'T': 3}


    def collect_samples(self, replicate_populations, sample_sizes):
        """
        Testing for concordance of segregating loci among samples requires that
        the samples be gathered in advance. Collects samples from replicate_populations

        :param replicate_populations: Multi-replicate population to analyze
        :param sample_sizes: Size of sample to gather.

        :note: :py:func:`len(sample_sizez)` == number of samples gathered from each replicate.

        :param str run_id: Identifier
        :return: List of populations
        """
        samples = {}
        for rep in replicate_populations.populations():
            samples[rep.dvars().rep] =\
                [sim.sampling.drawRandomSample(rep, sizes=sample_size)
                 for sample_size in sample_sizes]
        return samples


    def save_sample_populations(self, library_of_samples):
        for rep, sample_list in library_of_samples.items():
            for sample in sample_list:
                name = '_'.join([self.run_id, str(rep),
                                 str(sample.popSize())]) + '.pop'
                sample.save(os.path.join(os.getcwd(), name))


    def load_sample_library(self, number_of_replicates, sample_sizes, sub_run_id=''):
        reloaded_sample_library = {rep: [] for rep in range(number_of_replicates)}
        for rep in range(number_of_replicates):
            for size in sample_sizes:
                reloaded_sample_library[rep].append(sim.loadPopulation(self.run_id + '_' + str(rep) + '_' + str(size) + '.pop'))
        return reloaded_sample_library

    def load_particular_sample_library(self, replicates_to_load, sample_sizes, sub_run_id=''):
        particular_loaded_sample_library = {rep: [] for rep in replicates_to_load}
        for rep in replicates_to_load:
            for size in sample_sizes:
                particular_loaded_sample_library[rep].append(sim.loadPopulation(self.run_id+ '_' + str(rep) + '_' + str(size) + '.pop'))
        return particular_loaded_sample_library

    def calculate_power_fpr(self, panel_map, sample_sizes, number_of_replicates,
                                number_of_qtl):
        """
        Determines the power by calculating number of detected loci divided by
        the number of loci with effects.

        :param panel_map: Dictionary of dictionaries of pandas.DataFrames. Keyed by panel_map[size][rep] = pd.DataFrame
        :param sample_sizes: List of integers corresponding to how many individuals are sampled from each replicate.
        :param number_of_replicates: Number of replicates in the run
        :param number_of_qtl: Loci declared as QTL and assigned an effect
        :return: pd.DataFrame summarizing power and false positive rate across replicates and sample sizes.
        """


        true_positive_detected_loci = {}
        false_positive_detected_loci = {}

        analysis_columns = []
        for size in sample_sizes:
            analysis_columns.append('_'.join(['power', str(size)]))
            analysis_columns.append('_'.join(['fpr', str(size)]))
        power_fpr_results = pd.DataFrame(index=range(number_of_replicates),
                                         columns=analysis_columns)

        for size in sample_sizes:
            panel_of_tassel_results = pd.Panel(panel_map[size])
            potential_false_positives = len(panel_of_tassel_results[0])
            for rep in range(number_of_replicates):
                power_column = 'power_' + str(size)
                fpr_column = 'fpr_' + str(size)
                power_fpr_results.ix[rep, power_column] = len(
                    panel_of_tassel_results[rep][
                        (panel_of_tassel_results[rep].ix[:, 'q'] < 0.05)
                        & (panel_of_tassel_results[rep].ix[:,
                           'difference'] > 0.0)]) / number_of_qtl

                true_positive_detected_loci[size, rep] = panel_of_tassel_results[rep][
                    (panel_of_tassel_results[rep].ix[:, 'q'] < 0.05)
                    & (panel_of_tassel_results[rep].ix[:, 'difference'] > 0.0)]

                power_fpr_results.ix[rep, fpr_column] = len(
                    panel_of_tassel_results[rep][
                        (panel_of_tassel_results[rep].ix[:, 'q'] < 0.05)
                        & (
                        panel_of_tassel_results[rep].ix[:, 'difference'] == 0.0)]) / (
                                                        potential_false_positives - number_of_qtl)

                false_positive_detected_loci[size, rep] = panel_of_tassel_results[rep][
                    (panel_of_tassel_results[rep].ix[:, 'q'] < 0.05)
                    & (panel_of_tassel_results[rep].ix[:, 'difference'] == 0.0)]

        return power_fpr_results, true_positive_detected_loci, \
               false_positive_detected_loci

    def probability_of_detection(self, allele_effects_table, sample_sizes,
                                     number_of_replicates,
                                     true_positives_detected):
        """
        Calculates the probability that a locus with an effect is detected.
        Probability of detection is defined as the number of times a locus is detected
        divided by the total number of realizations

        Example
        -------
        If the number of realizations is 200 and a locus is detected in all 200 realizations
         then its probability of detection is 1.0

        :param allele_effects_table: Allele effects table given by generate_allele_effects_table
        :param sample_sizes: List of number of individuals sampled from each replicate
        :param number_of_replicates: Number of replicates in the run
        :param true_positives_detected: Dictionary of lists of loci with effects that were detected.
        :return: Modified version of allele effects table which includes the probability of detection column.
        """

        number_of_realizations = len(sample_sizes) * number_of_replicates
        aggregate_detected_loci = []
        for size in sample_sizes:
            for rep in range(number_of_replicates):
                aggregate_detected_loci.extend(list(true_positives_detected[size, rep].index))

        prob_detection_table = allele_effects_table.copy()
        prob_detection_table = pd.DataFrame(prob_detection_table, columns=list(
            prob_detection_table.columns) + ['detected'])
        prob_detection_table.fillna(0, inplace=True)
        prob_detection_table.index = list(prob_detection_table['tassel_locus'])
        prob_detection_table.drop('tassel_locus', axis=1, inplace=True)
        for locus, count_detected in col.Counter(aggregate_detected_loci).items():
            prob_detection_table.ix[
                locus, 'detected'] = count_detected / number_of_realizations

        return prob_detection_table

    def seg_loci_among_samples(self, sample_library):
        """
        Examines the segregating loci of all samples and counts how many
        times each set of segregating loci occurs among samples. If the final
        count is equal to number of replicates * length of sample sizes then
        all loci have the same segregating loci.

        :param sample_library:
        :return:
        """
        for rep in sample_library.values():
            for sample in rep:
                sim.stat(sample, numOfSegSites=sim.ALL_AVAIL, vars=['segSites', 'numOfSegSites'])
        seg_of_samples = (tuple(sample.dvars().segSites) for rep in
                          sample_library.values() for sample in rep)
        segregating_loci_counts = col.Counter(seg_of_samples)
        return segregating_loci_counts


    def store_allele_frequencies(self, library_of_samples, alleles):

        hdf_store = pd.HDFStore(self.run_id + '_allele_frequency_storage.h5')
        for rep_id, samples in library_of_samples.items():
            for sample in samples:
                af = allele_data(sample, alleles,
                                     range(sample.totNumLoci()))

                name = self.run_id + '/' + str(rep_id) + '/' + str(sample.popSize())
                hdf_store.put(name, af)

        hdf_store.close()



def write_multiple_sample_analyzer(library_of_samples, sample_size_list,
                             quantitative_trait_loci, alleles, allele_effects,
                         heritability, segregating_loci, run_id='infinite',
                         sub_run_id = '',
                         allele_frequency_hdf=''):


    syn_parameters = shelve.open('synthesis_parameters')
    int_to_snp_map = syn_parameters['integer_to_snp']
    syn_parameters.close()

    allele_frqs = pd.HDFStore(allele_frequency_hdf)

    for rep_id, sample_list in library_of_samples.items():
        for sample_population in sample_list:

            operators.assign_additive_g(sample_population, quantitative_trait_loci,
                                        allele_effects)
            operators.calculate_error_variance(sample_population, heritability)
            operators.calculate_p(sample_population)

            name = run_id + sub_run_id + '_' + str(rep_id) + '_' + str(sample_population.popSize())
            afrq_name = run_id + '/' + str(rep_id) + '/' + str(sample_population.popSize())
            minor_alleles = allele_frqs[afrq_name]['minor_allele']
            minor_allele_frequencies = np.array([allele_frqs[afrq_name]['minor_frequency'][locus] for locus in segregating_loci])

            gwas = GWAS(sample_population, segregating_loci, run_id)

            indir = "C:\\tassel\\input\\"

            ccm = gwas.calculate_count_matrix(minor_alleles)
            ps_svd = gwas.pop_struct_svd(ccm)
            gwas.population_structure_formatter(ps_svd  , indir + name + '_structure_matrix.txt')
            gwas.hapmap_formatter(int_to_snp_map,
                                         indir + name + '_simulated_hapmap.txt')
            gwas.trait_formatter(indir + name + '_phenotype_vector.txt')
            gwas.calc_kinship_matrix(ccm, minor_allele_frequencies,
                                            indir + name + '_kinship_matrix.txt')

            gwas.replicate_tassel_gwas_configs(str(rep_id), sample_population.popSize(),
                indir + name + '_simulated_hapmap.txt',
                indir + name + '_kinship_matrix.txt',
                indir + name + '_phenotype_vector.txt',
                indir + name + '_structure_matrix.txt',
                "C:\\tassel\\output\\" + name + '_out_',
                                              "C:\\Users\DoubleDanks\\BISB\\wisser\\code\\rjwlab-scripts\\"
                                              "saegus_project\\devel\\magic\\1478\\" + run_id + "_gwas_pipeline.xml")

    allele_frqs.close()


def generate_allele_effects_frequencies(sample_population,
                                        alleles,
                                        qtl,
                                        allele_effects):

    alpha_alleles = alleles[:, 0]
    beta_alleles = alleles[:, 1]

    alpha_allele_frequencies = np.array(
        [sample_population.dvars().alleleFreq[locus][alpha_allele] for
         locus, alpha_allele in enumerate(alpha_alleles)])
    beta_allele_frequencies = np.array(
        [sample_population.dvars().alleleFreq[locus][beta_allele] for
         locus, beta_allele in enumerate(beta_alleles)])

    alpha_effects, beta_effects = np.zeros(len(alpha_alleles)), np.zeros(
        len(beta_alleles))
    for locus in qtl:
        alpha_effects[locus] = allele_effects[locus][alpha_alleles[locus]]
        beta_effects[locus] = allele_effects[locus][beta_alleles[locus]]

    difference = np.abs(alpha_effects - beta_effects)
    loci = np.arange(len(alpha_alleles), dtype=np.int64)

    table_data = dict(locus=loci,
                      alpha_allele=alpha_alleles,
                      alpha_frequency=alpha_allele_frequencies,
                      alpha_effect=alpha_effects,
                      beta_allele=beta_alleles,
                      beta_frequency=beta_allele_frequencies,
                      beta_effect=beta_effects,
                      difference=difference)

    order_of_columns = ['locus', 'alpha_allele', 'alpha_frequency',
                        'alpha_effect',
                        'beta_allele', 'beta_frequency', 'beta_effect',
                        'difference']

    return pd.DataFrame(table_data, columns=order_of_columns)


def store_allele_effect_frequency_tables(sample_library, alleles,
                                         qtl,
                                         allele_effects,
                                         run_id, sub_run_id):

    alpha_alleles = alleles[:, 0]
    beta_alleles = alleles[:, 1]

    store_name = '_'.join([run_id, sub_run_id, 'allele_effects_and_frequencies.h5'])
    multi_aef_storage = pd.HDFStore(store_name)
    for rep_id, samples in sample_library.items():
        for sample in samples:
            expanded_allele_effects_and_frequencies = generate_allele_effects_frequencies(sample, alleles, qtl, allele_effects)
            name = '/'+ '/'.join([run_id, sub_run_id, str(rep_id), str(sample.popSize())])
            multi_aef_storage.put(name, expanded_allele_effects_and_frequencies)
    multi_aef_storage.close()


def generate_super_table(run_id,
                         rep_id,
                         sample_size,
                         segregating_loci,
                         sub_run_id=''):


    """
    Combines the TASSEL output with allele frequencies, allele effects and
    q values.


    """
    gwas_results_file_name = run_id +'_'+ sub_run_id+'_'+str(rep_id)+'_'+str(sample_size)+'_'+'out_2.txt'
    gwas_results = pd.read_csv(gwas_results_file_name, sep='\t')
    gwas_results.drop('Trait', axis=1, inplace=True)
    gwas_results.drop('Pos', axis=1, inplace=True)
    gwas_results.drop(0, axis=0, inplace=True)
    gwas_results = gwas_results.ix[:, 'Marker':'dom_p']
    gwas_results.index = gwas_results.index - 1
    gwas_results.drop('Marker', axis=1, inplace=True)

    q_values_file_name = run_id +'_'+sub_run_id +'_'+ str(rep_id) +'_'+ str(sample_size)+ '_' + 'qvalues.txt'
    qvalues = pd.read_csv(q_values_file_name, sep='\t')
    qvalues.columns = ['q']
    qvalues.index = qvalues.index - 1

    results = gwas_results.join(qvalues)

    allele_frequency_table = reload_allele_frequencies_table(run_id, rep_id, sample_size, sub_run_id)
    subsetted_af_table = remap_allele_frequency_table_loci(allele_frequency_table, segregating_loci)

    sub_results = results.join(subsetted_af_table)

    allele_effects_and_frequencies_table = reload_allele_effects_and_frequencies_table(run_id, sub_run_id, rep_id, sample_size)
    subsetted_aefrq_table = remap_allele_effect_and_frq_table_loci(allele_effects_and_frequencies_table, segregating_loci)

    super_results = sub_results.join(subsetted_aefrq_table)

    return super_results

def collect_power_analysis_data(run_id, sample_sizes, number_of_replicates,
                                    segregating_loci, sub_run_id):


    """
    Results will calculated using super tables. No longer reads a file,
    looks at a table in memory.


    :param sample_sizes:
    :param number_of_replicates:
    :param genome_wide_allele_effect_differences:
    :return:
    """
    panel_map = {}
    for size in sample_sizes:
        panel_map[size] = {}
        for rep in range(number_of_replicates):
            panel_map[size][rep] = generate_super_table(run_id, rep, size,
                                                        segregating_loci,
                                                        sub_run_id)
    return panel_map



def reload_allele_frequencies_table(run_id, rep_id, sample_size, sub_run_id=''):
    if sub_run_id != '':
        store_name = '_'.join([run_id, 'allele_frequency_storage.h5'])
    else:
        store_name = '_'.join([run_id, 'storage_diff.h5'])
    table_name = '/' + '/'.join([run_id, str(rep_id), str(sample_size)])
    reloaded_table = pd.read_hdf(store_name, key=table_name)
    return reloaded_table

def remap_allele_frequency_table_loci(allele_frequency_table,
                              segregating_loci):

    loci = list(allele_frequency_table.index)
    droppable_loci = [locus for locus in loci if locus not in segregating_loci]
    allele_frequency_table_subset = copy.deepcopy(allele_frequency_table)
    allele_frequency_table_subset.drop(droppable_loci, axis=0, inplace=True)
    new_index = list(range(len(segregating_loci)))
    allele_frequency_table_subset.index = new_index
    return allele_frequency_table_subset


def remap_allele_effect_and_frq_table_loci(
            allele_effect_and_frequency_table,
            segregating_loci):
    """
    Drops non-segregating loci from the allele effect - frequency hybrid table.

    :param allele_effect_and_frequency_table:
    :param saegus_to_tassel_loci:
    :return:
    """


    loci = list(allele_effect_and_frequency_table.index)
    droppable_loci = [locus for locus in loci if locus not in segregating_loci]
    allele_effect_and_frequency_table_subset = allele_effect_and_frequency_table.drop(
        droppable_loci,
        axis=0)
    new_index = list(range(len(segregating_loci)))
    allele_effect_and_frequency_table_subset.index = new_index
    return allele_effect_and_frequency_table_subset


def reload_allele_effects_and_frequencies_table(run_id, sub_run_id, rep_id, sample_size):
    store_name = '_'.join([run_id, sub_run_id, 'allele_effects_and_frequencies.h5'])
    table_name = '/' + '/'.join([run_id, sub_run_id, str(rep_id), str(sample_size)])
    reloaded_table = pd.read_hdf(store_name, key=table_name)
    return reloaded_table


def reconfigure_gwas_results(gwas_results_file, q_values_file,
                             expanded_allele_effects_table, delim="\t"):
    """
    Useful function which parses the output from TASSEL, collects the useful
    pieces of information such as p values, q values and allele effects
    into a useful table.
    """

    gwas_results = pd.read_csv(gwas_results_file, sep=delim)
    gwas_results.drop('Trait', axis=1, inplace=True)
    gwas_results.drop('Pos', axis=1, inplace=True)
    gwas_results.drop(0, axis=0, inplace=True)
    gwas_results = gwas_results.ix[:, 'Marker':'dom_p']
    gwas_results.index = gwas_results.index - 1
    gwas_results.drop('Marker', axis=1, inplace=True)
    qvalues = pd.read_csv(q_values_file, sep=delim)
    qvalues.columns = ['q']
    qvalues.index = qvalues.index - 1
    results = gwas_results.join(qvalues)
    greater_results = results.join(expanded_allele_effects_table)

    return greater_results


def calculate_power_fpr(panel_map, sample_sizes, number_of_replicates,
                            number_of_qtl):
    analysis_columns = []
    for size in sample_sizes:
        analysis_columns.append('_'.join(['power', str(size)]))
        analysis_columns.append('_'.join(['fpr', str(size)]))

    power_fpr_results = pd.DataFrame(index=range(number_of_replicates),
                                     columns=analysis_columns)
    for size in sample_sizes:
        panel_of_tassel_results = pd.Panel(panel_map[size])
        potential_false_positives = len(panel_of_tassel_results[0])
        for rep in range(number_of_replicates):
            power_column = 'power_' + str(size)
            fpr_column = 'fpr_' + str(size)
            power_fpr_results.ix[rep, power_column] = len(
                panel_of_tassel_results[rep][
                    (panel_of_tassel_results[rep].ix[:, 'q'] < 0.05)
                    & (panel_of_tassel_results[rep].ix[:,
                       'difference'] > 0.0)]) / number_of_qtl

            power_fpr_results.ix[rep, fpr_column] = len(
                panel_of_tassel_results[rep][
                    (panel_of_tassel_results[rep].ix[:, 'q'] < 0.05)
                    & (panel_of_tassel_results[rep].ix[:,
                       'difference'] == 0.0)]) / potential_false_positives - number_of_qtl
    return power_fpr_results



def replicate_gwas_results(gwas_results_file, q_values_file, expanded_allele_effects_table, delim="\t"):
    """
    Useful function which parses the output from TASSEL, collects the useful
    pieces of information such as p values, q values and allele effects
    into a useful table.
    """


    gwas_results = pd.read_csv(gwas_results_file, sep=delim)
    gwas_results.drop('Trait', axis=1, inplace=True)
    gwas_results.drop('Pos', axis=1, inplace=True)
    gwas_results.drop(0, axis=0, inplace=True)
    gwas_results = gwas_results.ix[:, 'Marker':'p']
    gwas_results.index = gwas_results.index - 1
    gwas_results.drop('Marker', axis=1, inplace=True)
    qvalues = pd.read_csv(q_values_file, sep=delim)
    qvalues.columns = ['q']
    qvalues.index = qvalues.index - 1
    results = gwas_results.join(qvalues)
    greater_results = results.join(expanded_allele_effects_table)

    return greater_results

class Synbreed(object):

    def __init__(self, genotypes, phenotypes, genetic_map):
        self.genotypes = genotypes
        self.phenotypes = phenotypes
        self.genetic_map = genetic_map

    def genotype_formatter(self, hap_map, syn_genotype_filename):
        """
        Modifies ``hap_map`` to obtain the format for Synbreed.
        Subtracts extraneous columns, transposes the result, adds 1 to each
        index to convert Python index to R index. Then prepends each locus
        with "M" indicating 'marker'.
        :param hap_map: A simulated hapmap file from hapmap_formatter.
        :param syn_genotype_filename: Name of the file which the modified
        hapmap is written to.

        """

        syn_geno_map = hap_map.drop(["rs", "alleles", "chrom", "pos", "strand",
                   "assembly", "center", "protLSID", "assayLSID",
                  "panelLSID", "QCcode"], axis=1).T
        syn_geno_map.columns = ["M" + str(idx + 1) for idx in syn_geno_map.columns]

        syn_geno_map.to_csv(syn_genotype_filename, sep="\t", index=True,
                            header=False)

        return syn_geno_map


class MetaData(object):
    """
    The wgs is extensively paramterized. Hence changing one parameter will potentially produce a significantly different
    result in the final population. Therefore, a set of replications is defined by a particular of parameterization.
    The parameterization will be described in a metadata document. The class MetaData is responsible for collecting
    the parameterization information and processing it into a writable file.
    """

    def __init__(self, prefounders, founders, population_sizes, allele_effect_information,
                 allele_effects_table, metadata_filename):
        """
        An instance of MetaData should have enough information to completely specify a population without using any
        external information.
        :param prefounders: Prefounder population of the 26 lines which were used to make the NAM population.
        :param founders: Subset of prefounders used to make a derived population.
        :param population_sizes: Size of the population during the F_one, F_two, 'mate-and-merge' phase and finally
        the selection phase.
        :param allele_effect_information: Information about the distribution of allele effects and the corresponding
        parameters, the random number generator package and random seed used to generate the allele effects.
            Ex: normal(0, 1), numpy.random, seed 1337.
        :param allele_effects_table: The actual tabular/dictionary representation of the realized allele effect values.
        """
        self.prefounders = prefounders
        self.founders = founders
        self.population_sizes = population_sizes
        self.allele_effect_information = allele_effect_information
        self.allele_effects_table = allele_effects_table
        self.metadata_filename = metadata_filename


        # A master function will use other functions to write the necessary information to file.

    @staticmethod
    def ascii_chromosome_representation(pop, reduction_factor, metadata_filename):
        """
        Writes a ascii representation of chromosomes with uninteresting loci
        as * and QTL as |. The representation is has scale 1 /
        reduction_factor to make it viable to put into a .txt document.
        :param pop:
        :param reduction_factor:
        :param metadata_filename:
        """
        reduced_chromosomes = [math.floor(chrom/reduction_factor) for chrom in list(pop.numLoci())]
        reduced_qtl = [math.floor(pop.chromLocusPair(locus)[1]/reduction_factor) for locus in pop.dvars().properQTL]
        chromosomes_of_qtl = [pop.chromLocusPair(qtl)[0] for qtl in pop.dvars().properQTL]
        aster_chroms = [['*']*chrom_len for chrom_len in reduced_chromosomes]
        for red_qtl, chrom_of_qtl in zip(reduced_qtl, chromosomes_of_qtl):
            aster_chroms[chrom_of_qtl][red_qtl] = '|'
        with open(metadata_filename, 'a') as chrom_file:
            chrom_file.write('Scale: 1/%d\n' % reduction_factor)
            for idx, chrom in enumerate(aster_chroms):
                idx += 1
                chrom_file.write('Chromosome: %d\n' % idx)
                chrom_file.write(''.join(chrom) + '\n')

    @staticmethod
    def coefficient_of_dispersion(pop):
        """
        Mean to variance ratio of pairwise sequential distances of quantitative trait loci.
        Note that this statistic contributes nothing if there is only one qtl on a chromosome.
        :param pop:
        """
        chrom_loc_pairs = [pop.chromLocusPair(pop.dvars().properQTL[i]) for i in range(len(pop.dvars().properQTL))]
        chromosomes = [chrom_loc_pairs[i][0] for i in range(len(chrom_loc_pairs))]
        diffs = []
        for i in range(len(chrom_loc_pairs)):
            if chromosomes[i-1] == chromosomes[i]:
                diffs.append(math.fabs(chrom_loc_pairs[i-1][1] - chrom_loc_pairs[i][1]))
        diffs = np.array(diffs)
        mean = np.mean(diffs)
        var = np.var(diffs)
        var_to_mean_ratio = var/mean
        return var_to_mean_ratio

    @staticmethod
    def genomic_dispersal(pop):
        """
        Genomic dispersal is a novel statistics which measures the spread of loci over a genome.z All loci of a chromosome
        are compared to the center of the genetic map (in cMs) and weighted by the length of that chromosome.
        :param pop: Population used for recurrent selection
        :return: Dimensionless constant describing the parameterization
        """
        chrom_loc_pairs = [pop.chromLocusPair(pop.dvars().properQTL[i]) for i in range(len(pop.dvars().properQTL))]
        chromosomes = [chrom_loc_pairs[i][0] for i in range(len(chrom_loc_pairs))]
        qtl_positions = [(chrom_loc_pairs[i][1]) for i in range(len(chrom_loc_pairs))]
        chromosome_midpoints = {i: (pop.numLoci()[i]/2) for i in range(len(pop.numLoci()))}
        diffs = []
        for pos, chrom, i in zip(qtl_positions, chromosomes, range(len(chrom_loc_pairs))):
            diffs.append(pos - chromosome_midpoints[chrom])
        squared_diffs = np.square(np.array(diffs))
        root_squared_diffs = np.sqrt(squared_diffs)
        denominator_lengths = np.array(list(pop.numLoci()))
        pre_genetic_dispersal = np.divide(root_squared_diffs, denominator_lengths)
        genomic_dispersal = sum(pre_genetic_dispersal)
        return genomic_dispersal