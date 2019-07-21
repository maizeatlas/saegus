__author__ = 'John J. Dougherty III'
__project__ = 'wgs'
# -*- coding: utf-8 -*-
import simuPOP as sim
import math
import numpy as np
import pandas as pd
import collections as col
import csv
from scipy import linalg
import matplotlib.pyplot as plt
plt.ioff()


def allele_effects_writer(pop, filename):
    ae = pop.dvars().alleleEffects
    with open(filename, 'w') as ae_file:
        ae_writer = csv.writer(ae_file, delimiter=',')
        for k, v in ae.items():
            ae_writer.writerow([k, v])


def qtl_file_writer(pop, filename):
    ae = pop.dvars().alleleEffects
    with open(filename, 'w') as qtl_file:
        qtl_writer = csv.writer(qtl_file, delimiter=',')
        for k in ae.keys():
            qtl_writer.writerow([k[0]])


def ae_file_reader(pop, filename):
    ae_dict = col.OrderedDict()
    with open(filename, newline='\n') as ae_file:
        ae_reader = csv.reader(ae_file, delimiter=',', quotechar='"')
        for row in ae_reader:
            splitrows = row[0].split(",")
            stripped_sublocus = splitrows[0].lstrip('(')
            stripped_nucleotide = splitrows[1].rstrip(')')
            sl = int(stripped_sublocus)
            nt = int(stripped_nucleotide)
            effect = float(row[1])
            ae_dict[sl, nt] = effect
    pop.dvars().alleleEffects = ae_dict


def population_information_writer(pop: sim.Population, info_field_list: list, filename: str):
    """
    Writes information equivalent to a simuPOP.utils Exporter; however, this function is much more
     flexible. I am uncertain of what information will be necessary in the future.
    :param pop: sim.Population
    :param info_field_list: List of information fields to be included in output file.
    :param filename: Name of output file.
    :return: CSV File 'filename'
    """
    with open(filename, 'w') as info:
        info_writer = csv.writer(info, delimiter=',')
        info_writer.writerow(info_field_list)
        for ind in pop.individuals():
            info_writer.writerow([ind.ind_id, ind.mother_id, ind.father_id, ind.ge, ind.pe,
                                  ind.fitness, ind.genotype()])


def meta_population_information_writer(meta_pop: sim.Population, info_field_list: list, filename: str):
    """
    Writes information equivalent to a simuPOP.utils Exporter; however, this function is much more
     flexible. I am uncertain of what information will be necessary in the future.
    :param info_field_list: List of information fields to be included in output file.
    :param filename: Name of output file.
    :return: CSV File 'filename'
    """

    for i in range(meta_pop.numSubPop()):
        full_filename = filename + "_" + str(i) + ".txt"
        with open(full_filename, 'w') as info:
            info_writer = csv.writer(info, delimiter=',')
            info_writer.writerow(info_field_list)
            for ind in meta_pop.individuals(i):
                info_writer.writerow([ind.ind_id, ind.mother_id, ind.father_id, ind.ge, ind.pe])



def generate_lineage_of_qtl(pop):
    qtl = pop.dvars().qtl
    if pop.numSubPop() == 1:
        pop.splitSubPop(0, sizes=[pop.popSize()-1, 1])
    linray = lin_callback(pop.individual(0, 1), qtl)
    for ind in pop.individuals(0):
        linray = np.vstack((linray, lin_callback(ind, qtl)))
    sim.mergeSubPops(pop)
    return linray


def generate_lineage_of_all(pop):
    qtl = pop.dvars().qtl
    if pop.numSubPop() == 1:
        pop.splitSubPop(0, sizes=[pop.popSize()-1, 1])
    linray = lin_callback(pop.individual(0, 1), qtl)
    for ind in pop.individuals(0):
        linray = np.vstack((linray, lin_callback(ind, qtl)))
    sim.mergeSubPops(pop)
    return linray


def lin_callback(ind: sim.Individual, qtl: list):
    """
    Callback function which allows for vectorized computation of frequency data of each individual.
    :param ind: Individual of population
    :param qtl: List of absolute indicies of quantitative trait loci
    :return:
    """
    return np.array([ind.lineage()[qtl_idx] for qtl_idx in qtl])


def triplet_calc_along_axis(lineage_array, pop):
    """
    Usage
    -----
    Uses numpy vectorized computations to create a list of triplet effects.

    Parameters
    ----------
    :param lineage_array: numpy.array of the lineage values which are triplets
    :param pop: simuPOP.Population object
    :return: numpy.array of equal axial dimension as input array.
    """
    return np.apply_along_axis(triplet_effect_list, 1, lineage_array, pop)


def triplet_effect_list(triplets, pop):
    triplet_effects = [pop.dvars().alleleEffects.ix[idx-1, str(trip)[1]] +
                       pop.dvars().alleleEffects.ix[idx, str(trip)[2]] +
                       pop.dvars().alleleEffects.ix[idx+1, str(trip)[3]]
                       for idx, trip in zip(pop.dvars().qtl, triplets)]
    return triplet_effects


def calc_triplet_frequencies(pop):
    arr = np.split(pop.lineage(), 2, 1)
    alpha = arr[0]
    beta = arr[1]
    splitlet_frequency_dict = col.defaultdict(int)
    for locus in range(3*pop.totNumLoci()):
        alpha_slice = alpha[:, locus]
        beta_slice = beta[:, locus]
        combined = np.append(alpha_slice, beta_slice)
        splitlet_frequency_dict[locus] = col.defaultdict(int, col.Counter(combined))
        for splitlet in iter\
                    (splitlet_frequency_dict[locus].keys()):
            splitlet_frequency_dict[locus][splitlet] = splitlet_frequency_dict[locus][splitlet]/800
    return splitlet_frequency_dict


def calc_af_changes(pop):
    af_changes = pd.DataFrame(index=pop.dvars().properQTL, columns=list(range(pop.dvars().gen)))
    offspring_afrqs = pd.Panel(pop.dvars().current_frqs).fillna(0)
    parental_afrqs = pd.Panel(pop.dvars().parental_frqs).fillna(0)
    for i in range(pop.dvars().gen):
        delta = (parental_afrqs[i]-offspring_afrqs[i]).max()
        af_changes.ix[:, i] = delta
    return af_changes


def chromosome_abs_and_rel(pop):
    chromosomes = col.OrderedDict()
    relative_indexes = col.OrderedDict()
    for locus in pop.dvars().properQTL:
        chromosomes[locus] = pop.chromLocusPair(locus)[0]
        relative_indexes[locus] = pop.chromLocusPair(locus)[1]
    return chromosomes, relative_indexes


def index_to_cM(index):
    if index == 0:
        return -1*0.2
    if index == 1:
        return 0.0
    if index == 2:
        return 0.2


def arbitrary_ordering_of_founders(prefounders, pop, ordered_founder_indices):
    """
    Allows user to enter pairs of founders regardless of their absolute index.
    """
    assert tuple([float(id) for id in list(range(1, 27))]) == prefounders.indInfo('ind_id'), "Prefounders must have ind_ids 1 through 26."
    pop.addInfoFields('custom')
    complement = [j for j in range(prefounders.popSize()+1) if j not in ordered_founder_indices]
    prefounders.removeIndividuals(complement)
    for ord_idx, i in zip(ordered_founder_indices, range(len(ordered_founder_indices))):
        pop.indByID(ord_idx).custom = i
    pop.sortIndividuals(infoFields='custom')


def table_for_development(prefounders: sim.Population, founders: list, pop: sim.Population, replicate_type: str):
    """
    Function to generate a summary data table for population subjected to either drift or selection.

    :param prefounders: Population from which all populations are derived
    :param founders: List of integers specifying which prefounders to use in making a population
    :param pop: Population subjected to either drift or selection
    :param replicate_type: String either 'drift' or 'selection' reflecting whether 'pop' is subjected to drift or
    selection
    :return: hierarchically indexed pandas.DataFrame with frequency and effect data for all qtl
    """
    data = []
    columns = ['QTL', 'Chromosome', 'cM_Position', 'Singlet', 'SingletEffect', 'SingletFrequency', 'Triplet',
               'TripletEffect', 'NumberOfFounders']
    generational_columns = [pop.dvars().generations[i] for i in range(pop.dvars().gen)]
    columns.extend(generational_columns)
    chromosomes, relative_positions = chromosome_abs_and_rel(pop)
    qtls = []
    stf = test_subloci_calc_freq(prefounders, pop)
    for i in range(len(pop.dvars().properQTL)):
        qtls.append(prefounders.dvars().properQTL[i]-1)
        qtls.append(prefounders.dvars().properQTL[i])
        qtls.append(prefounders.dvars().properQTL[i]+1)
    for qtl in pop.dvars().properQTL:
        for trip in sorted(list(pop.vars()['triplet_freq'][0][qtl].keys())):
            for i in range(3):
                datarow = [qtl,  # locus
                           chromosomes[qtl],  # chromosome of locus
                           relative_positions[qtl]+index_to_cM(i), # position in cM
                           trip[i],  # singlet
                           pop.dvars().alleleEffects[qtl+i-1, int(trip[i])],  # singlet effect
                           stf[qtl, qtl+i-1][trip[i]],  # singlet frequency
                           trip,  # triplet
                           sum([pop.dvars().alleleEffects[qtl-1, int(trip[0])],  # triplet effect
                           pop.dvars().alleleEffects[qtl, int(trip[1])],
                           pop.dvars().alleleEffects[qtl+1, int(trip[2])]]),
                           len(founders)]  # number of founders selected from prefounders
                frq_per_gen = [pop.vars()['triplet_freq'][j][qtl][trip] for j in range(8)]
                rep_frq_per_gen = [pop.vars()[replicate_type + '_triplet_freq'][k][qtl][trip] for k in range(8, pop.dvars().gen)]
                frq_per_gen.extend(rep_frq_per_gen)
                datarow.extend(frq_per_gen)
                data.append(datarow)
    s = pd.DataFrame(data, columns=columns)
    for i in range(8, 13):
        for term in pop.dvars().mav_terms:
            s.set_value(term, pop.dvars().generations[i], pop.dvars().mavs[i, term])
    return s.fillna(0)


def nucleotide_translator(nucleotide):
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


def all_snps_setter(pop, map_file_name):
    """
    Updated version of genotype_setter. This function sets all SNPs involve in the lineage triplets as regular alleles
    of the pop genotype. Using all_snps_setter requires the user to specify t
    :param pop:
    :param map_file_name:
    :return:
    """
    genotype_table = pd.read_table(map_file_name, index_col=0)
    founder_names = genotype_table.columns[3:]

    central_snps = [i for i in range(len(genotype_table))
                    if float(genotype_table['cM_pos'][i]) == float(int(genotype_table['cM_pos'][i]))]
    left_flank = [i - 1 for i in central_snps]
    right_flank = [i + 1 for i in central_snps]
    combined = []
    combined.extend(central_snps)
    combined.extend(left_flank)
    combined.extend(right_flank)
    combined = sorted(combined)

    all_genotypes = [[nucleotide_translator(genotype_table.ix[idx, name][0])for idx in combined] +
                    [nucleotide_translator(genotype_table.ix[idx, name][1]) for idx in combined]
                     for name in founder_names]

    for ind, gt in zip(pop.individuals(), all_genotypes):
        ind.setGenotype(gt)


def partial_stringer(lin):
    yield list(map(lambda triplet: str(triplet), lin))


def stringer(pop):
    reformatted_lineage = []
    complete = []
    for ind in pop.individuals():
        incomplete_lineage = list(map(lambda triplet: str(triplet), ind.lineage()))
        reformatted_lineage.append(incomplete_lineage)
    for i in range(len(reformatted_lineage)):
        temp = reformatted_lineage[i]
        complete.append(list(map(lambda triplet: triplet[1:4], temp)))
    return complete

def iterator_stringer(pop):
    reformatted_lineage = []
    complete = []
    for ind in pop.individuals():
        incomplete_lineage = list(map(lambda triplet: str(triplet), ind.lineage()))
        reformatted_lineage.append(incomplete_lineage)
    for i in iter(reformatted_lineage):
        complete.append(list(map(lambda triplet: triplet[1:4], i)))
    return complete


def yielder_stringer(ind):
    yield list(map(lambda triplet: str(triplet), ind.lineage()))


def triplet_splitter(trip):
    return [trip[0], trip[1], trip[2]]


def master_splitter(reformatted_lineage_array):
    complete_splitlets = []
    full_length = len(reformatted_lineage_array)
    for i in range(full_length):
        current_lineage = reformatted_lineage_array.pop()
        splitlet_generator = map(triplet_splitter, current_lineage)
        singlets = []
        for j in range(2956):
            singlets.extend(next(splitlet_generator))
        complete_splitlets.append(singlets)
    return complete_splitlets


def lineage_triplets_to_splitlets(pop):
    reformed_lineage = stringer(pop)
    complete_splitlets = master_splitter(reformed_lineage)
    return np.flipud(np.array(complete_splitlets))


def splitlet_frequencies(pop, singlet_lineage):
    arr = np.split(singlet_lineage, 2, 1)
    alpha = arr[0]
    beta = arr[1]
    splitlet_frequency_dict = col.defaultdict(int)
    for locus in range(3*pop.totNumLoci()):
        alpha_slice = alpha[:, locus]
        beta_slice = beta[:, locus]
        combined = np.append(alpha_slice, beta_slice)
        splitlet_frequency_dict[locus] = col.defaultdict(int, col.Counter(combined))
        for splitlet in iter\
                    (splitlet_frequency_dict[locus].keys()):
            splitlet_frequency_dict[locus][splitlet] = splitlet_frequency_dict[locus][splitlet]/800
    return splitlet_frequency_dict


def find_minor_alleles(pop: sim.Population):
    minor_allele_frequencies = col.OrderedDict()
    minor_alleles = np.empty(pop.totNumLoci())
    for locus in range(pop.totNumLoci()):
        if len(list(pop.dvars().alleleFreq[locus].keys())) >= 2:
            minor_allele_frequencies[locus] = min(pop.dvars().alleleFreq[locus].values())
            for k in pop.dvars().alleleFreq[locus].keys():
                if pop.dvars().alleleFreq[locus][k] == minor_allele_frequencies[locus]:
                    minor_alleles[locus] = k
        elif len(list(pop.dvars().alleleFreq[locus].keys())) == 1:
            # If an allele is at fixation then we want to show a zero for the MAC format.
            minor_alleles[locus] = 7
    return minor_alleles


def find_major_alleles(pop):
    """
    Finds major allele of each locus and builds a dictionary of locus, major allele pairs.
    :param allele_frequencies:
    :return:
    """
    major_allele_frequencies = col.OrderedDict()
    major_alleles = col.OrderedDict()
    for locus in range(pop.totNumLoci()):
        major_allele_frequencies[locus] = max(pop.dvars().alleleFreq[locus].values())
        for k in pop.dvars().alleleFreq[locus].keys():
            if pop.dvars().alleleFreq[locus][k] == major_allele_frequencies[locus]:
                major_alleles[locus] = k
    return major_alleles


def minor_allele_genotype_matrix(pop, minor_allele_dict, mac_genotype_filename):
    print("Constructing minor allele genotype matrix.")
    minor_allele_list = np.array(list(minor_allele_dict.values()), dtype=int)
    minor_allele_matrix = pd.DataFrame(np.zeros((pop.popSize(), pop.totNumLoci())))
    for i, ind in enumerate(pop.individuals()):
        cmp_one = np.equal(minor_allele_list, ind.genotype(ploidy=0), dtype=int)
        cmp_two = np.equal(minor_allele_list, ind.genotype(ploidy=1), dtype=int)
        minor_allele_count = np.add(cmp_one, cmp_two, dtype=int)
        minor_allele_matrix.ix[i, :] = minor_allele_count
    index = list(map(int, pop.indInfo('ind_id')))
    minor_allele_matrix.index = index
    genframe = pd.DataFrame(list(pop.indInfo('generation')), index=index, columns=['popdata'])
    minor_allele_matrix = pd.concat([genframe, minor_allele_matrix], axis=1)
    print("Writing minor allele genotype matrix to filename: %s" % mac_genotype_filename)
    minor_allele_matrix.to_csv(mac_genotype_filename, float_format='%d', sep='\t')
    print("Finished.")
    return minor_allele_matrix

class Synbreed(object):

    def genotype_writer(self, pop, minor_allele_dict, singlet_lineage, genotype_filename):
        minor_allele_list = np.array(list(minor_allele_dict.values()), dtype=int)
        # option to add the index of 'ind_id' -- may be more suitable for databases
        minor_allele_matrix = pd.DataFrame(np.zeros((pop.popSize(), 3*pop.totNumLoci())))
        arr = np.split(singlet_lineage, 2, axis=1)
        alpha = np.array(arr[0], dtype=int)
        beta = np.array(arr[1], dtype=int)
        for i in range(pop.popSize()):
            cmp_one = np.equal(minor_allele_list, alpha[i, :], dtype=int)
            cmp_two = np.equal(minor_allele_list, beta[i, :], dtype=int)
            minor_allele_count = np.add(cmp_one, cmp_two, dtype=int)
            minor_allele_matrix.ix[i, :] = minor_allele_count
        minor_allele_matrix.index = ['ID' + str(int(ind.ind_id)) for ind in pop.individuals()]
        minor_allele_matrix.columns = ['M' + str(i) for i in range(4434)]
        minor_allele_matrix.to_csv(genotype_filename, index=True, header=True, sep='\t')

    def phenotype_writer(self, pop, phenotype_filename):
        """
        Writes information equivalent to a simuPOP.utils Exporter; however, this function is much more
         flexible. I am uncertain of what information will be necessary in the future.
        :param filename: Name of output file.
        :return: CSV File 'filename'
        """
        with open(phenotype_filename, 'w') as info:
            info_writer = csv.writer(info, delimiter='\t')
            info_writer.writerow(['ind_id', 'pe'])
            for ind in pop.individuals():
                info_writer.writerow(['ID' + str(int(ind.ind_id)), ind.pe])

    def sex_as_int_to_string_converter(self, sex_as_int):
        if sex_as_int == 1:
            return 'M'
        elif sex_as_int == 2:
            return 'F'

    def pedigree_writer(self, pop, pedigree_filename):
        with open(pedigree_filename, 'w') as pedigree:
            pedigree_writer = csv.writer(pedigree, delimiter='\t')
            pedigree_writer.writerow(['ind_id', 'mother_id', 'father_id', 'generation'])
            for ind in pop.individuals():
                pedigree_writer.writerow(['ID'+str(int(ind.ind_id)), 'ID' + str(int(ind.mother_id)),
                                          'ID' + str(int(ind.father_id)), ind.generation])

    def marker_map(self, pop, marker_map_filename):
        splitlet_loci_per_chromosome = [3*number_loci for number_loci in pop.numLoci()]
        centi_morgan_positions = []
        for i in range(pop.totNumLoci()):
            centi_morgan_positions.append(pop.locusPos(i) - 0.2)
            centi_morgan_positions.append(pop.locusPos(i))
            centi_morgan_positions.append(pop.locusPos(i) + 0.2)
        corresponding_splitlet_chromosome_list = []
        for chrom_index, number_loci in zip(range(1, 11), splitlet_loci_per_chromosome):
            corresponding_splitlet_chromosome_list.extend([chrom_index]*number_loci)
        marker_prefixed_index = ["M" + str(i) for i in range(3*pop.totNumLoci())]
        syn_map = pd.DataFrame([list(np.array(corresponding_splitlet_chromosome_list, dtype=int)), centi_morgan_positions]).T
        syn_map.columns = ['chr', 'pos']
        syn_map.index = marker_prefixed_index
        syn_map.to_csv(marker_map_filename, sep="\t", index=True, names=['chr', 'pos'], header=True)

    def write_population_data(self, pop, singlet_lineage, minor_allele_dict, genotype_filename='syn_genotype.txt',
                                    phenotype_filename='syn_phenotype.txt', pedigree_filename='syn_pedigree.txt'):
        self.syn_genotype_writer(pop, minor_allele_dict, singlet_lineage, genotype_filename)
        self.syn_phenotype_writer(pop, phenotype_filename)
        self.syn_pedigree_writer(pop, pedigree_filename)



def write_info_fields_to_file(pop:sim.Population, info_fields: list, output_filename: str, fmt_string):
    """
    Utilize numpy's simple numpy.savetxt function to write various types of pop.IndInfo(info_string) to files.
    Superior solution to using csv.writer's for simple .txt output.
    + Example
    ----------
    pop.infoFields() = ('ind_id', 'ge', 'pe')
    info_fields = ['ind_id', 'ge', 'pe']
    info_array = np.array([pop.indInfo(info) for info in info_fields])
    np.savetxt(output_filename, info_array, delimiter='\t', newline='\n', header=">Generation: pop.dvars().gen",
        comments='')
    """
    info_array = np.array([pop.indInfo(info) for info in info_fields]).T
    np.savetxt(output_filename, info_array, fmt='%d\t%.2f\t%.2f', delimiter='\t', newline='\n', header=">Generation: 0\nID\tG\tP",
               comments='')



# Matplotlib code to shade an interval.
"""
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cbook as cbook

# load up some sample financial data
datafile = cbook.get_sample_data('goog.npy')
r = np.load(datafile).view(np.recarray)

# create two subplots with the shared x and y axes
fig, (ax1, ax2) = plt.subplots(1,2, sharex=True, sharey=True)

pricemin = r.close.min()

ax1.plot(r.date, r.close, lw=2)
ax2.fill_between(r.date, pricemin, r.close, facecolor='blue', alpha=0.5)

for ax in ax1, ax2:
    ax.grid(True)

ax1.set_ylabel('price')
for label in ax2.get_yticklabels():
    label.set_visible(False)

fig.suptitle('Google (GOOG) daily closing price')
fig.autofmt_xdate()
"""


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
        Writes a ascii representation of chromosomes with uninteresting loci as * and QTL as |.
        The representation is has scale 1 / reduction_factor to make it viable to put into a .txt document.
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


class PCA(object):
    """
    Class for performing principal component analyis on genotype matrices. Test for population structure significance
    tests the largest eigenvalue of the genotype covarience matrix. Details can be found in the paper:
    Population Structure and Eigenanalysis Patterson et al 2006.
    """
    def __init__(self, population_size, number_of_markers):
        self.population_size = population_size
        self.number_of_markers = number_of_markers

    def svd(self, major_allele_matrix):
        """
        Follows procedure of Population Structure and Eigenanalysis Patterson et al 2006.
        Constructs a genotype matrix of bi-allelic loci where each entry is the number of copies of the major allele at
        each locus. The genotype matrix has dimensions (number_of_individuals)*(number_of_markers)
        """
        if 'popdata' in list(major_allele_matrix.columns):
            major_allele_matrix.drop('popdata', axis=1)
        shifting_factor = np.apply_along_axis(np.mean, axis=1, arr=major_allele_matrix)
        p_vector = np.divide(shifting_factor, 2)
        scaling_factor = np.sqrt(np.multiply(p_vector, (1-p_vector)))
        corrected_matrix = np.array(pd.DataFrame(list(map(lambda i:
                                                 (major_allele_matrix.ix[major_allele_matrix.index[i], :]-shifting_factor[i])/scaling_factor[i],
                                                 range(self.population_size)))))
        # singular value decomposition using scipy linalg module
        eigenvectors, s, v = linalg.svd(corrected_matrix)
        eigenvalues = np.diagonal(np.square(linalg.diagsvd(s, self.population_size, self.number_of_markers))).T
        sum_of_eigenvalues = np.sum(eigenvalues)
        fraction_of_variance = np.divide(eigenvalues, sum_of_eigenvalues)
        eigenvalues = np.vstack((eigenvalues, fraction_of_variance))
        return eigenvectors, eigenvalues

    def shift_and_scale(self, major_allele_matrix):
        shifting_factor = np.apply_along_axis(np.mean, axis=1, arr=major_allele_matrix)
        p_vector = np.divide(shifting_factor, 2)
        scaling_factor = np.sqrt(np.multiply(p_vector, (1-p_vector)))
        corrected_matrix = np.array(pd.DataFrame(list(map(lambda i:
                                                 (major_allele_matrix.ix[i, :]-shifting_factor[i])/scaling_factor[i],
                                                 range(self.population_size)))))
        return corrected_matrix

    @staticmethod
    def triplet_genotype_matrix(pop, major_allele_dict, singlet_lineage):
        major_allele_list = np.array(list(major_allele_dict.values()), dtype=int)
        # option to add the index of 'ind_id' -- may be more suitable for databases
        major_allele_matrix = pd.DataFrame(np.zeros((pop.popSize(), 3*pop.totNumLoci())))
        arr = np.split(singlet_lineage, 2, axis=1)
        alpha = np.array(arr[0], dtype=int)
        beta = np.array(arr[1], dtype=int)
        for i in range(pop.popSize()):
            cmp_one = np.equal(major_allele_list, alpha[i, :], dtype=int)
            cmp_two = np.equal(major_allele_list, beta[i, :], dtype=int)
            major_allele_count = np.add(cmp_one, cmp_two, dtype=int)
            major_allele_matrix.ix[i, :] = major_allele_count
        return major_allele_matrix

    @staticmethod
    def major_allele_genotype_matrix(pop, major_allele_dict, mac_genotype_filename):
        print("Constructing major allele genotype matrix.")
        major_allele_list = np.array(list(major_allele_dict.values()))
        major_allele_matrix = pd.DataFrame(np.zeros((pop.popSize(), pop.totNumLoci())))
        for i, ind in enumerate(pop.individuals()):
            cmp_one = np.equal(major_allele_list, ind.genotype(ploidy=0), dtype=np.int)
            cmp_two = np.equal(major_allele_list, ind.genotype(ploidy=1), dtype=int)
            major_allele_count = np.add(cmp_one, cmp_two, dtype=int)
            major_allele_matrix.ix[i, :] = major_allele_count
        index = list(map(int, pop.indInfo('ind_id')))
        major_allele_matrix.index = index
        genframe = pd.DataFrame(list(pop.indInfo('generation')), index=index, columns=['popdata'])
        major_allele_matrix = pd.concat([genframe, major_allele_matrix], axis=1)
        print("Writing major allele genotype matrix to filename: %s" % mac_genotype_filename)
        major_allele_matrix.to_csv(mac_genotype_filename, float_format='%d', sep='\t')
        print("Finished.")
        return major_allele_matrix

    @staticmethod
    def test_statistic(eigenvalues, m_individuals):
        sum_of_eigenvalues = np.sum(eigenvalues[:, 0])
        n_hat_numerator = (m_individuals + 1)*sum_of_eigenvalues
        n_hat_denom = (m_individuals - 1)*sum_of_eigenvalues - sum_of_eigenvalues
        n_hat = n_hat_numerator/n_hat_denom
        lowercase_l = (m_individuals - 1)*eigenvalues[0, 0]
        mu_hat = np.square((np.sqrt(n_hat - 1) + np.sqrt(m_individuals))) / n_hat
        sigma_hat = ((np.sqrt(n_hat - 1) + np.sqrt(m_individuals))/n_hat) * \
                    (((1/np.sqrt(n_hat - 1)) + 1/np.sqrt(m_individuals)) ** (1 / 3.0))
        test_statistic = (lowercase_l - mu_hat) / sigma_hat
        return test_statistic

    @staticmethod
    def genotype_writer(pop, major_allele_genotype_matrix, major_allele_genotype_filename):
        major_allele_genotype_matrix.index = ['ID' + str(int(ind.ind_id)) for ind in pop.individuals()]
        major_allele_genotype_matrix.columns = ['M' + str(i) for i in range(3*pop.totNumLoci())]
        major_allele_genotype_matrix.to_csv(major_allele_genotype_filename, index=True, header=True, sep='\t')

    @staticmethod
    def write_eigenanalysis_to_file(eigenvectors, eigenvalues, eigenvector_filename, eigenvalue_filename):
        np.savetxt(eigenvector_filename, eigenvectors, fmt='%3f', delimiter='\t')
        np.savetxt(eigenvalue_filename, eigenvalues.T, fmt='%3f', delimiter='\t')

    @staticmethod
    def plot_principal_components(eigenvectors, output_filename):
        f, ax = plt.subplots(1, 1, figsize=(8, 8))
        f.suptitle('Population Structure: PCA')
        ax.set_xlabel('First Eigenvector')
        ax.set_ylabel('Second Eigenvector')
        #generational_eigenvectors = np.split(eigenvectors, number_of_generations, axis=1)
        #for i in range(number_of_generations):
        ax.plot(eigenvectors[:100, 0], eigenvectors[:100, 1], 'r*', linewidth=0.0, markersize=5, alpha=0.5)
        ax.plot(eigenvectors[100:200, 0], eigenvectors[100:200, 1], 'c^', linewidth=0.0, markersize=5, alpha=0.5)
        ax.plot(eigenvectors[200:300, 0], eigenvectors[200:300, 1], 'go', linewidth=0.0, markersize=5, alpha=0.5)
        ax.plot(eigenvectors[300:400, 0], eigenvectors[300:400, 1], 'y8', linewidth=0.0, markersize=5, alpha=0.5)
        ax.plot(eigenvectors[400:, 0], eigenvectors[400:, 1], 'mv', linewidth=0.0, markersize=5, alpha=0.5)
        plt.savefig(output_filename, dpi=300)

    @staticmethod
    def plot_multiple_replicate_population_structure(eigenvectors, output_filename):
        """
        Generates a plot of population stucture of multiple replicates of the meta-population. Parameter
        eigenvectors is a list of (pop_size)*(2) numpy.arrays because only the first two eigenvectors
        are considered. Generations are given different markers and colors for purposes of comparison.
        """
        import matplotlib.lines as mlines
        f, axes = plt.subplots(2, 2, figsize=(16, 16))
        f.suptitle("Population Structure Under Five Different Scenarios")

        itervecs = iter(eigenvectors)
        for i in range(2):
            for j in range(2):
                e = next(itervecs)
                axes[i, j].plot(e[:100, 0], e[:100, 1], 'r*',
                                label='Generation One', linewidth=0.0, markersize=5, alpha=0.5)
                axes[i, j].plot(e[100:200, 0], e[100:200, 1], 'c^',
                                label='Generation Two', linewidth=0.0, markersize=5, alpha=0.5)
                axes[i, j].plot(e[200:300, 0], e[200:300, 1], 'go',
                                label='Generation Three', linewidth=0.0, markersize=5, alpha=0.5)
                axes[i, j].plot(e[300:400, 0], e[300:400, 1], 'y8',
                                label='Generation Four', linewidth=0.0, markersize=5, alpha=0.5)
                axes[i, j].plot(e[400:, 0], e[400:, 1], 'mv',
                                label='Generation Five', linewidth=0.0, markersize=5, alpha=0.5)
        axes[0, 0].set_xlabel('First Eigenvector')
        axes[0, 0].set_ylabel('Second Eigenvector')

        red_stars = mlines.Line2D([], [], color='red', marker='*', markersize=8, label='Generation One')
        cyan_triangles = mlines.Line2D([], [], color='cyan', marker='^', markersize=8, label='Generation Two')
        green_circles = mlines.Line2D([], [], color='green', marker='o', markersize=8, label='Generation Three')
        yellow_octagons = mlines.Line2D([], [], color='yellow', marker='8', markersize=8, label='Generation Four')
        inverted_triangles = mlines.Line2D([], [], color='m', marker='v', markersize=8, label='Generation Five')

        labels = ['Generation One', 'Generation Two', 'Generation Three', 'Generation Four', 'Generation Five']
        plt.figlegend(handles=[red_stars, cyan_triangles,
                               green_circles, yellow_octagons,
                               inverted_triangles], labels=labels, loc='upper right')

        axes[0, 0].set_title('Selection with Random Mating')
        axes[0, 1].set_title('Drift with Random Mating, Largest')
        axes[1, 0].set_title('Selection without Random Mating')
        axes[1, 1].set_title('Drift without Random Mating')
        plt.savefig(output_filename, dpi=300)


def extract_genotype_matrix(pop, geno_filename):
    """
    Extracts genotypes of all individuals and stores them as a
    ploidy x population_size x number_loci array. Creates number of files equal to ploidy of the genotype
    information and returns a local copy.

    simuPOP stores the genotype information for the entire population as a single large custom
    data type called a ``carray``. Granted simuPOP comes with functions to access individual
    genotypes this can be cumbersome. This function solves that problem by providing a numpy array
    of the genotype values.
    """
    pop_ploidy = pop.ploidy()
    pop_size = pop.popSize()
    genotype_array = np.empty((pop_ploidy, pop.popSize(), pop.totNumLoci()))
    base_genotypes = np.array(pop.genotype())
    split_genotypes = np.split(base_genotypes, pop_ploidy)
    for i in range(pop_ploidy):
        genotype_array[i] = np.array(np.split(split_genotypes[i], pop_size))
        filename_with_ploidy_prefix = 'ploidy_' + str(i) + '_' + geno_filename
        np.savetxt(filename_with_ploidy_prefix, genotype_array[i], fmt='%d', delimiter='\t')
    return genotype_array




