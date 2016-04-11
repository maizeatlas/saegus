# -*- coding: utf-8 -*-
import simuPOP as sim
import math
import numpy as np
import pandas as pd
import collections as col
import os
import random
import shelve
from scipy import linalg
from . import operators, parameters

import parameters


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
    sim.stat(pop, alleleFreq=loci, vars=['alleleFreq'])

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
            np.array([pop.dvars().alleleFreq[locus][allele]])
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

    minor_alleles = [allele_frq['minor', 'alleles'][locus] for locus in loci]
    minor_frequencies = [pop.dvars().alleleFreq[locus][minor_allele] for locus, minor_allele in enumerate(minor_alleles)]
    major_alleles = [allele_frq['major', 'alleles'][locus] for locus in loci]
    major_frequencies = [pop.dvars().alleleFreq[locus][major_allele] for locus, major_allele in enumerate(major_alleles)]

    allele_data_structure = \
        pd.DataFrame(np.array([minor_alleles, minor_frequencies, major_alleles, major_frequencies]),
                     index=['minor_allele', 'minor_frequency', 'major_allele', 'major_frequency'],
                     columns=loci).T

    return allele_data_structure




def rank_allele_effects(self, pop, loci, alleles,
                       allele_effects):
    """
    Collects information about alleles at quantitative trait loci into a
    dictionary. Determines favorable/unfavorable allele and corresponding
    frequency. Keys of quantitative_trait_alleles have similar hierarchy
    for both the alleles and their frequencies.
    :param pop:
    :param loci:
    :param alleles:
    :param allele_effects:


    """
    quantitative_trait_alleles = {}
    quantitative_trait_alleles['effects'] = col.OrderedDict()
    quantitative_trait_alleles['alleles'] = col.OrderedDict()
    quantitative_trait_alleles['alleles']['favorable'] = col.OrderedDict()
    quantitative_trait_alleles['alleles']['unfavorable'] = col.OrderedDict()
    quantitative_trait_alleles['frequency'] = col.OrderedDict()
    quantitative_trait_alleles['frequency']['favorable'] = col.OrderedDict()
    quantitative_trait_alleles['frequency']['unfavorable'] = col.OrderedDict()
    for locus in loci:
        temp_effects = []
        for allele in alleles[locus]:
            temp_effects.append(allele_effects[locus][allele])
        quantitative_trait_alleles['effects'][locus] = temp_effects

    for locus in loci:
        for allele in alleles[locus]:
            if allele_effects[locus][allele] == max(
                    quantitative_trait_alleles['effects'][locus]):
                quantitative_trait_alleles['alleles']['favorable'][locus] =\
                    allele
            if allele_effects[locus][allele] == min(
                    quantitative_trait_alleles['effects'][locus]):
                quantitative_trait_alleles['alleles']['unfavorable'][locus] =\
                    allele

    for locus, allele in quantitative_trait_alleles['alleles']['favorable'].items():
        quantitative_trait_alleles['frequency']['favorable'][locus] = \
            pop.dvars().alleleFreq[locus][allele]
    for locus, allele in quantitative_trait_alleles['alleles']['unfavorable'].items():
        quantitative_trait_alleles['frequency']['unfavorable'][locus] =\
            pop.dvars().alleleFreq[locus][allele]

    return quantitative_trait_alleles

def allele_frq_table(self, pop, number_gens,
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


def generate_allele_effects_table(qtl, alleles, allele_effects):
    """
    Creates a simple pd.DataFrame for allele effects. Hard-coded
    for bi-allelic case.

    :parameter list qtl: List of loci declared as QTL
    :parameter np.array alleles: Array of alleles at each locus
    :parameter dict allele_effects: Mapping of effects for alleles at each QTLocus
    """
    ae_table = {
        'locus': [],
        'alpha_allele': [],
        'alpha_effect': [],
        'beta_allele': [],
        'beta_effect': [],
    }

    for locus in qtl:
        ae_table['locus'].append(locus)
        alpha_allele, beta_allele = alleles[locus]
        ae_table['alpha_allele'].append(alpha_allele)
        ae_table['beta_allele'].append(beta_allele)
        alpha_effect = allele_effects[locus][alpha_allele]
        ae_table['alpha_effect'].append(alpha_effect)
        beta_effect = allele_effects[locus][beta_allele]
        ae_table['beta_effect'].append(beta_effect)
    order_of_columns = ['locus', 'alpha_allele', 'alpha_effect', 'beta_allele', 'beta_effect']
    allele_effect_frame = pd.DataFrame(ae_table, columns=order_of_columns)
    return allele_effect_frame


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


def plot_frequency_vs_effect(pop, haplotype_table, plot_title,
                             plot_file_name,
                             color_map='Dark2'):
    """
    Uses the haplotype data table to arrange data into a chromosome
    color coded multiple generation plot which shows the change in
    haplotype frequency over time. Haplotypes are dots with fixed
    x-position which shows their effect. Their motion along the y-axis
    which is frequency shows changes over time.
    :param plot_title:
    :param plot_file_name:
    :param color_map:
    :param pop:
    :param haplotype_table:
    """

    plt.style.use('ggplot')

    distinct_chromosomes = list(set(haplotype_table['chromosome']))
    number_of_different_colors = len(distinct_chromosomes)
    generation_labels = ['G_' + '{' + str(i) + '}' for i in
                          range(0, 2*(pop.numSubPop()), 2)]
    generations = ['G_' + str(i) for i in range(0, 2*(pop.numSubPop()), 2)]

    c_map = plt.get_cmap(color_map)

    colors = c_map(np.linspace(0, 1, number_of_different_colors))

    chromosome_colors = {distinct_chromosomes[i]: colors[i] for i in
                         range(number_of_different_colors)}

    effect_frq_by_chromosome = {}

    for sp in range(pop.numSubPop()):
        effect_frq_by_chromosome[sp] = {}
        for chrom in distinct_chromosomes:
            haplotype_frequencies = np.array(
                haplotype_table.loc[
                    haplotype_table['chromosome'] == chrom][generations[sp]])

            haplotype_effects = np.array(
                haplotype_table.loc[
                    haplotype_table['chromosome'] == chrom]['effect'])

            effect_frq_by_chromosome[sp][chrom] = np.array([
                haplotype_frequencies, haplotype_effects])

    # Figure parameters
    maximum_haplotype_effect = max(haplotype_table['effect'])

    generations = ['G_'+str(i) for i in range(0, 2*(pop.numSubPop()) + 1, 2)]

    f, ax = plt.subplots(pop.numSubPop(), 1, figsize=(15, 40))
    for sp in range(pop.numSubPop()):
        ax[sp].set_xlim(-0.5, maximum_haplotype_effect+4)
        ax[sp].set_ylim(-0.1, 1.1)
        for chrom in distinct_chromosomes:
            ax[sp].plot(effect_frq_by_chromosome[sp][chrom][1],
                    effect_frq_by_chromosome[sp][chrom][0],
                    markersize=8, linewidth=0.0, marker='*',
                    color=chromosome_colors[chrom],
                        label="Chrom {}".format(chrom))
        #handles, labels = ax[sp].get_legend_handles_labels()
        ax[sp].set_xlabel("Effect")
        ax[sp].set_ylabel("Frequency")
        ax[sp].set_title(r'${gen}$'.format(gen=generation_labels[sp]),
                         fontsize=12)
        ax[sp].legend(loc='best')
    f.suptitle(plot_title,
               fontsize=24)

    f.savefig(plot_file_name, dpi=300)

    return effect_frq_by_chromosome







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


class GWAS(object):
    """
    Class for performing principal component analyis on genotype matrices.
    Test for population structure significance tests the largest eigenvalue
    of the genotype covarience matrix. Details can be found in the paper:
    Population Structure and Eigenanalysis Patterson et al 2006.
    """
    def __init__(self, pop, loci, allele_subset, run_id):
        self.pop = pop
        self.loci = loci
        self.allele_subset = allele_subset
        self.run_id = run_id

        self.individual_names = ['I' + str(ind.ind_id)[:-2]
                                 for ind in pop.individuals()]

        self.locus_names = list(range(len(loci)))
        self.pos_names = list(range(len(loci)))

    def calculate_count_matrix(self, count_matrix_filename):
        """
        A function to calculate the copy numbers of either the minor or
        major allele for each individual at each locus. Minor or major
        alleles parameter is a single set of alleles which determines if the
        return is the minor or major allele count matrix.
        :param pop:
        :param allele_subset:
        :param count_matrix_filename:
        """
        comparison_array = np.array([self.allele_subset[locus]
                                     for locus in self.loci], dtype=np.int8)

        count_matrix = np.zeros((self.pop.popSize(), len(self.loci)))


        for i, ind in enumerate(self.pop.individuals()):
            alpha_genotype = np.array([ind.genotype(ploidy=0)[locus]
                              for locus in self.loci])

            alpha_comparisons = np.equal(comparison_array, alpha_genotype,
                                         dtype=np.int8)


            beta_genotype = [ind.genotype(ploidy=1)[locus]
                             for locus in self.loci]

            beta_comparisons = np.equal(comparison_array, beta_genotype,
                                         dtype=np.int8)

            counts = np.add(alpha_comparisons, beta_comparisons, dtype=np.int8)
            count_matrix[i, :] = counts

        np.savetxt(count_matrix_filename, count_matrix, fmt="%d")
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
        eigen_data = {}
        eigen_data['vectors'] = eigenvectors
        eigen_data['values'] = eigenvalues
        eigen_data['fraction_variance'] = fraction_of_variance
        return eigen_data

    def population_structure_formatter(self, eigen_data, structure_filename):
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

        header = "<Covariate>\t\t\n<Trait>\td1\td2\n"

        cwd = os.getcwd()
        file_out_path = os.path.join(cwd, structure_filename)

        if os.path.exists(file_out_path):
            os.remove(file_out_path)
        with open(structure_filename, 'w') as f:
            f.write(header)
            structure_matrix.to_csv(f, sep='\t', index=True, header=False)

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

    def hapmap_formatter(self, int_to_snp_conversions, hapmap_matrix_filename):
        """
        Converts genotype data from sim.Population object to HapMap file format
        in expectation to be used in TASSEL for GWAS. At present the column
        names will be hardcoded as will some of the values.
        ``hapmap_matrix_filename`` is the name of the file the formatted
        matrix will be written to.
        :param int_to_snp_conversions:
        :param hapmap_matrix_filename:
        :return:
        :rtype:
        """
        hapmap_data = {}
        hapmap_data['rs'] = self.locus_names
        hapmap_data['alleles'] = ['NA']*len(self.loci)
        hapmap_data['chrom'] = [self.pop.chromLocusPair(locus)[0]+1 for
                                locus in self.loci]
        hapmap_data['pos'] = self.pos_names

        # Several columns which are set to 'NA'.
        extraneous_columns = ['strand', 'assembly', 'center', 'protLSID',
                              'assayLSID', 'panelLSID', 'QCcode']
        for column in extraneous_columns:
            hapmap_data[column] = ['NA']*len(self.loci)

        # Each individual has a column. Simulated individuals will have names
        # reflecting some information about them. 'RS' recurrent selection,
        # 'R' rep, 'G' generation and 'I' ind_id

        # Finally each individual's genotype must be converted to HapMap format


        for i, ind in enumerate(self.pop.individuals()):
            alpha_genotype = np.array([ind.genotype(ploidy=0)[locus]
                                       for locus in self.loci])

            beta_genotype = [ind.genotype(ploidy=1)[locus]
                             for locus in self.loci]

            hapmap_data[self.individual_names[i]] = [
                ''.join(sorted(int_to_snp_conversions[a]
                               + int_to_snp_conversions[b]))
                                 for a, b in zip(alpha_genotype, beta_genotype)]


        # Need to guarantee that the column names are in same order as the
        # genoype data. Iterate over individuals in population to build up a
        #  list of names will guarantee that col names are in same order as
        # the hapmap_data

        hapmap_ordered_columns = ['rs', 'alleles', 'chrom', 'pos', 'strand',
                           'assembly', 'center', 'protLSID', 'assayLSID',
                               'panelLSID', 'QCcode'] + self.individual_names

        hapmap_matrix = pd.DataFrame(columns=hapmap_ordered_columns)
        for k, v in hapmap_data.items():
            hapmap_matrix[k] = v

        hapmap_matrix.to_csv(hapmap_matrix_filename, sep='\t',
                             index=False)

        return hapmap_matrix

    def trait_formatter(self, trait_filename):
        """
        Simple function to automate the formatting of the phenotypic data.
        Because of the way the header must be written the file is opened in
        append mode. Rewriting to the same file many times could introduce an
        unsuspected bug.
        :param trait_filename:
        """
        header = "<Trait>\tsim\n"

        trait_vector = pd.DataFrame(np.array([self.individual_names,
                               self.pop.indInfo('p')]).T)

        cwd = os.getcwd()
        file_out_path = os.path.join(cwd, trait_filename)

        if os.path.exists(file_out_path):
            os.remove(file_out_path)
        with open(trait_filename, 'w') as f:
            f.write(header)
            trait_vector.to_csv(f, sep='\t', index=False, header=False)

        return trait_vector

    def calc_kinship_matrix(self, allele_count_matrix, allele_data,
                            kinship_filename):
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

        # Second allele in the unselected, un-inbred base population.
        # Refers to major allele in G_0
        allele_frequencies = np.array([allele_data['minor_frequency'][locus] for locus in self.loci])

        P = 2*(allele_frequencies - 0.5)

        Z = M - P

        scaling_terms = np.zeros((len(self.loci)))
        for idx, probability in enumerate(allele_frequencies):
            scaling_terms[idx] = 2*probability*(1 - probability)

        scaling_factor = sum(scaling_terms)

        G = Z*Z.T/scaling_factor

        annotated_G = pd.DataFrame(G, index=self.individual_names)

        # Tassel example has number of individuals in the header of the G
        # matrix file
        header = "{}\n".format(self.pop.popSize())

        cwd = os.getcwd()
        file_out_path = os.path.join(cwd, kinship_filename)

        if os.path.exists(file_out_path):
            os.remove(file_out_path)
        with open(kinship_filename, 'w') as f:
            f.write(header)
            annotated_G.to_csv(f, sep='\t', index=True, header=False,
                               float_format='%.3f')

        return annotated_G

    def generate_tassel_gwas_configs(self, sample_size,
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


        lxml_root.write("C:\\tassel\\bin\\daoko_girl_sim_gwas_pipeline.xml", encoding="UTF-8",
                       method="xml", xml_declaration=True, standalone='',
                        pretty_print=True)

def population_sample_analyzer(full_population, sample_size,
                               quantitative_trait_loci, alleles,
                                       allele_effects, heritability,
                                       run_id='daoko_girl'):

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
  #  quantitative_trait_loci = sorted(
   #     random.sample(sample_population.dvars().segSites, number_of_qtl))
   # add_trait = parameters.Trait()

    #aes = add_trait.assign_allele_effects(alleles, quantitative_trait_loci,
     #                                     dist_function,
      #                                    *dist_func_parameters,
       #                                   multiplicity=multiplicity)
    aes_table = generate_allele_effects_table(quantitative_trait_loci,
                                                      alleles, allele_effects)
    operators.assign_additive_g(full_population, quantitative_trait_loci, allele_effects)
    operators.calculate_error_variance(sample_population, heritability)
    operators.phenotypic_effect_calculator(sample_population)
    af = allele_data(sample_population, alleles,
                             range(sample_population.totNumLoci()))

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
                  "C:\\tassel\\bin\\" + str(sample_size) + "_daoko_girl_out_",
                  "C:\\Users\DoubleDanks\\BISB\\wisser\\code\\rjwlab-scripts\\"
              "saegus_project\\devel\\magic\\1478\\daoko_girl_gwas_pipeline.xml")
    return segregating_loci, aes_table

def remap_ae_table_loci(allele_effect_table, saegus_to_tassel_loci):
    """
    Converts the absolute indices of saegus population to the indices of the truncated
    set of segregating loci. Allows for appropriate comparisons of loci.
    """
    remapped_loci = [saegus_to_tassel_loci[locus]
                     for locus in allele_effect_table['locus']]
    allele_effect_table['locus'] = remapped_loci

    allele_effect_table['difference'] = np.abs(allele_effect_table['alpha_effect'] -
                                          allele_effect_table['beta_effect'])

    allele_effect_table.index = allele_effect_table.locus

    expanded_ae_table = pd.DataFrame(np.zeros((len(saegus_to_tassel_loci))).T,
                                              columns=['difference'])

    for qtlocus in allele_effect_table.locus:
        expanded_ae_table.ix[qtlocus, 'difference'] =  allele_effect_table.ix[qtlocus, 'difference']

    return expanded_ae_table

def reconfigure_gwas_results(gwas_results_file, q_values_file, expanded_allele_effects_table, delim="\t"):
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
