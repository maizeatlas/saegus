==============
:mod:`analyze`
==============



.. py:function:: allele_frequencies(pop, alleles, loci)

   Determine major and minor alleles in each generation the aggregate
   population. Generations in a meta-populations correspond to
   sub-populations.

   :param pop:
   :param loci:

.. py:function:: rank_allele_effects(pop, loci, alleles, allele_effects)

   Collects information about alleles at quantitative trait loci into a
   dictionary. Determines favorable/unfavorable allele and corresponding
   frequency. Keys of quantitative_trait_alleles have similar hierarchy
   for both the alleles and their frequencies.

   :param pop:
   :param loci:
   :param alleles:
   :param allele_effects:

.. py:function:: allele_frq_table(pop, number_gens, allele_frq_data, recombination_rates, genetic_map)

   Tabulates useful information about each locus and allele frequency

   :param pop: Population with multiple sub-populations. Usually represents multiple generations of recurrent selection or drift.
   :param int number_gens: Number of generations of selection or drift
   :param dict allele_frq_data: Allele frequency data and the major/minor alleles at each locus.
   :param list recombination_rates: Recombination rates for each locus in order.
   :param pd.DataFrame genetic_map: Chromosome:cM position correspondence.


   Columns
   -------

   + abs_idex: Absolute index of locus
   + chrom: Chromosome
   + locus: Relative index of locus
   + major: Major allele
   + minor: Minor allele
   + recom_rate: Probability of recombination immediately *after* locus
   + cM: centiMorgan position on genetic map
   + v: miniature diagram where `*` is a breakpoint and `|` is a non-recombining locus
   + generation_labels: Generation number prefixed by `G`, frequency of minor allele at locus


   Example
   -------

   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+
   | abs_index | chrom | locus | major | minor | recom_rate |  cM  | v   | G_0 | G_2 | G_4 |
   +===========+=======+=======+=======+=======+============+======+=====+=====+=====+=====+
   |    0      |   1   |   0   |   0   |   2   |     0.00   | -4.8 | `|` | 0.0 | 0.0 | 0.0 |
   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+
   |    1      |   1   |   1   |   2   |   0   |     0.01   | -4.6 | `|` | 0.0 | 0.0 | 0.0 |
   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+
   |    2      |   1   |   2   |   2   |   0   |     0.00   | -4.4 | `*` | 0.0 | 0.0 | 0.0 |
   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+
   |    3      |   1   |   3   |   1   |   3   |     0.00   | -4.2 | `|` | 0.1 | 0.2 | 0.3 |
   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+
   |    4      |   1   |   4   |   1   |   2   |     0.00   | -4.0 | `|` | 0.0 | 0.0 | 0.0 |
   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+
   |    5      |   1   |   5   |   5   |   4   |     0.00   | -3.8 | `|` | 0.0 | 0.0 | 0.0 |
   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+
   |    6      |   1   |   6   |   2   |   0   |     0.01   | -3.6 | `|` | 0.0 | 0.0 | 0.0 |
   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+
   |    7      |   1   |   7   |   1   |   3   |     0.00   | -3.4 | `*` | 0.0 | 0.0 | 0.0 |
   +-----------+-------+-------+-------+-------+------------+------+-----+-----+-----+-----+

   Usage
   -----

   .. code:: python

      allele_data = analyze.Frq(pop, triplet_qtloci[0], alleles_by_locus, qt_allele_effects[0])
      allele_frequencies = allele_data.allele_frequencies(pop, range(pop.totNumLoci())
      allele_frequency_table = selection_qtd.allele_frq_table(pop, 10, allele_frq_data, recombination_rates,
                                                            genetic_map)






.. py:function:: qt_allele_table(self, qt_alleles, allele_effects):
   Generates a pd.DataFrame object of data relevant to quantitative
   trait alleles across all generations.
   :param qt_alleles:
   :type qt_alleles:
   :param allele_effects:
   :type allele_effects:
   :return:
   :rtype:

.. py:function:: collect_haplotype_data(pop, allele_effects, quantitative_trait_loci)

    :param pop:
    :param allele_effects:
    :param quantitative_trait_loci:



.. py:function:: generate_haplotype_data_table(pop, haplotype_data)

    Generates a table for easy analysis and visualization of haplotypes,
    effects, frequencies and locations.


    :param pop:
    :type pop:
    :param haplotype_data:
    :type haplotype_data:
    :return:
    :rtype:

.. py:function:: plot_frequency_vs_effect(pop, haplotype_table, plot_title,
                             plot_file_name,
                             color_map='Dark2')

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

.. py:class:: MetaData(object)

    The wgs is extensively paramterized. Hence changing one parameter will potentially produce a significantly different
    result in the final population. Therefore, a set of replications is defined by a particular of parameterization.
    The parameterization will be described in a metadata document. The class MetaData is responsible for collecting
    the parameterization information and processing it into a writable file.

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


.. py:class::PCA

    Class for performing principal component analyis on genotype matrices.
    Test for population structure significance tests the largest eigenvalue
    of the genotype covarience matrix. Details can be found in the paper:
    Population Structure and Eigenanalysis Patterson et al 2006.

    def __init__(self, pop, loci, qt_data):
        self.pop = pop
        self.loci = loci
        self.qt_data = qt_data

    def __str__(self):
        return "Parameters: PopSize {}, Number of Loci: {}, " \
               "Keys of Data: {}.".format(self.pop.popSize(), len(self.loci),
                                         self.qt_data.keys())

    def calculate_count_matrix(self, pop, alleles, count_matrix_filename):
        """
        A function to calculate the copy numbers of either the minor or
        major allele for each individual at each locus. Minor or major
        alleles parameter is a single set of alleles which determines if the
        return is the minor or major allele count matrix.
        :param pop:
        :param alleles:
        :param count_matrix_filename:
        """
        comparison_array = [alleles[locus] for locus in range(pop.totNumLoci())]
        count_matrix = np.zeros((pop.popSize(), len(alleles)))
        for i, ind in enumerate(pop.individuals()):
            alpha = np.equal(np.array(comparison_array), ind.genotype(
                ploidy=0), dtype=np.int8)
            beta = np.equal(np.array(comparison_array), ind.genotype(ploidy=1),
                            dtype=np.int8)
            counts = np.add(alpha, beta, dtype=np.int8)
            count_matrix[i, :] = counts
        np.savetxt(count_matrix_filename, count_matrix, fmt="%d")
        return count_matrix

    def svd(self, pop, count_matrix):
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

        shift_matrix = np.zeros((pop.popSize(), pop.totNumLoci()))
        scale_matrix = np.zeros((pop.popSize(), pop.totNumLoci()))
        for i in range(pop.totNumLoci()):
            shift_matrix[:, i] = shift
            scale_matrix[:, i] = scale

        corrected_matrix = (count_matrix - shift_matrix)/scale_matrix
        # singular value decomposition using scipy linalg module
        eigenvectors, s, v = linalg.svd(corrected_matrix)
        eigenvalues = np.diagonal(np.square(linalg.diagsvd(s, pop.popSize(),
                                                           pop.totNumLoci()))).T
        sum_of_eigenvalues = np.sum(eigenvalues)
        fraction_of_variance = np.divide(eigenvalues, sum_of_eigenvalues)
        eigen_data = {}
        eigen_data['vectors'] = eigenvectors
        eigen_data['values'] = eigenvalues
        eigen_data['fraction_variance'] = fraction_of_variance
        return eigen_data

    def test_statistic(self, pop, eigenvalues):
        sum_of_eigenvalues = np.sum(eigenvalues)
        n_hat_numerator = (pop.popSize() + 1)*sum_of_eigenvalues
        n_hat_denom = (pop.popSize()-1)*sum_of_eigenvalues - sum_of_eigenvalues
        n_hat = n_hat_numerator/n_hat_denom
        lowercase_l = (pop.popSize() - 1)*eigenvalues[0]
        mu_hat = np.square((np.sqrt(n_hat - 1) +
                            np.sqrt(pop.popSize()))) / n_hat
        sigma_hat = ((np.sqrt(n_hat - 1) + np.sqrt(pop.popSize()))/n_hat) * \
                    (((1/np.sqrt(n_hat - 1)) + 1/np.sqrt(pop.popSize())) ** (
                        1 / 3.0))
        test_statistic = (lowercase_l - mu_hat) / sigma_hat
        return test_statistic

.. py:class::GWAS

    A class to collect and format all data in preparation for GWAS using TASSEL.

    def __init__(self, pop, individual_names, locus_names, positions, *args,
                 **kwargs):
        self.pop = pop
        self.individual_names = individual_names
        self.locus_names = locus_names
        self.positions = positions


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
        hapmap_data['alleles'] = ['NA']*self.pop.totNumLoci()
        hapmap_data['chrom'] = [self.pop.chromLocusPair(locus)[0]+1 for
                                locus in
                                range(self.pop.totNumLoci())]
        hapmap_data['pos'] = self.positions

        # Several columns which are set to 'NA'.
        extraneous_columns = ['strand', 'assembly', 'center', 'protLSID',
                              'assayLSID', 'panelLSID', 'QCcode']
        for column in extraneous_columns:
            hapmap_data[column] = ['NA']*self.pop.totNumLoci()

        # Each individual has a column. Simulated individuals will have names
        # reflecting some information about them. 'RS' recurrent selection,
        # 'R' rep, 'G' generation and 'I' ind_id

        # Finally each individual's genotype must be converted to HapMap format

        for ind in self.pop.individuals():
            alpha, beta = ind.genotype(ploidy=0), ind.genotype(ploidy=1)
            hapmap_data[self.individual_names[ind.ind_id]] = [
                int_to_snp_conversions[a]+int_to_snp_conversions[b]
                                 for a, b in zip(alpha, beta)]


        # Need to guarantee that the column names are in same order as the
        # genoype data. Iterate over individuals in population to build up a
        #  list of names will guarantee that col names are in same order as
        # the hapmap_data
        ordered_names = [self.individual_names[ind.ind_id] for ind in
                         self.pop.individuals()]

        hapmap_ordered_columns = ['rs', 'alleles', 'chrom', 'pos', 'strand',
                           'assembly', 'center', 'protLSID', 'assayLSID',
                               'panelLSID', 'QCcode'] + ordered_names

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
        header = "<Trait> sim\n"

        # Ensure phenotype and name are coming from the same individual


        phenotypes = []
        ind_names = []
        for ind in self.pop.individuals():
            phenotypes.append(ind.p)
            ind_names.append(self.individual_names[ind.ind_id])

        trait_vector = pd.DataFrame([ind_names, phenotypes]).T

        cwd = os.getcwd()
        file_out_path = os.path.join(cwd, trait_filename)

        if os.path.exists(file_out_path):
            os.remove(file_out_path)
        with open(trait_filename, 'a') as f:
            f.write(header)
            trait_vector.to_csv(f, sep=' ', index=False, header=False)

        return trait_vector

    def population_structure_formatter(self, eigen_data, structure_filename):
        """
        Writes the first two of the population structure matrix to a
        file. First column of the file is are names.
        :param structure_filename:
        :param eigen_data:
        """

        ordered_names = [self.individual_names[ind.ind_id] for ind in
                         self.pop.individuals()]

        structure_matrix = pd.DataFrame([list(eigen_data['vectors'][:, 0].T),
                                     list(eigen_data['vectors'][:, 1].T)]).T

        structure_matrix.index = ordered_names

        header = "<Covariate>\t\t\n<Trait>\td1\td2\n"

        cwd = os.getcwd()
        file_out_path = os.path.join(cwd, structure_filename)

        if os.path.exists(file_out_path):
            os.remove(file_out_path)
        with open(structure_filename, 'a') as f:
            f.write(header)
            structure_matrix.to_csv(f, sep='\t', index=True, header=False)

        return structure_matrix

    def calc_kinship_matrix(self, allele_counts, allele_data,
                            kinship_filename):
        """
        Calculates the kinship matrix according to VanRaden 2008:
        Efficient Methods to Compute Genomic Predictions and writes it to a
        file formatted for Tassel. The variable names try to be consistent
        with the variable names in the paper.

        The allele frequencies used for this function are with respect to
        the base population or G0: after random mating and before selection.
        :param allele_counts:
        :param allele_data:
        :param kinship_filename:
        :return:
        :rtype:
        """

        M = np.matrix(allele_counts - 1)

        major_allele_frequencies = \
            np.array([allele_data['major', 'frequency', 0][locus]
                      for locus in range(self.pop.totNumLoci())])

        P = 2*(major_allele_frequencies - 0.5)

        Z = M - P

        scaling_factor = sum(2*P*(1 - P))

        G = (Z*Z.T)/scaling_factor

        annotated_G = pd.DataFrame(G, index=[self.individual_names[ind.ind_id]
                                             for ind in
                                             self.pop.individuals()])

        # Tassel example has number of individuals in the header of the G
        # matrix file
        header = "{}\n".format(self.pop.popSize())

        cwd = os.getcwd()
        file_out_path = os.path.join(cwd, kinship_filename)

        if os.path.exists(file_out_path):
            os.remove(file_out_path)
        with open(kinship_filename, 'a') as f:
            f.write(header)
            annotated_G.to_csv(f, sep=' ', index=True, header=False)

        return annotated_G


def generate_tassel_gwas_configs(input_directory_prefix,
                                 out_directory_prefix,
                                 config_file_prefix,
                                 run_identifier_prefix,
                                 xml_pipeline_template):
    """
    Creates an xml file to run TASSEL using a mixed linear model approach.
    Assumes use of hapmap, kinship, phenotype and population structure files.




    The TASSEL command line interface requires a considerable number of
    options to run GWAS. It is impractical to run the command line manually
    for the number of replications in a simulated study. The TASSEL command
    line interface allows the user to input a .xml file with the same
    information which is used in the terminal.

    :param input_directory_prefix: Directory path to send the input files.
    :param run_identifier_prefix: Identifier for single replicate of data
    :param xml_pipeline_template: XML file already setup for running a
    specific kind of GWAS
    :return: XML file to run a single replicate of data using TASSEL
    """


    import xml.etree.ElementTree as ET
    import lxml.etree as etree

    tree = ET.parse(xml_pipeline_template)
    root = tree.getroot()
    lxml_tree = etree.fromstring(ET.tostring(root))
    lxml_root = lxml_tree.getroottree()

    lxml_root.find('fork1/h').text = os.path.join(input_directory_prefix,
                                                  run_identifier_prefix +
                                                  'simulated_hapmap.txt')

    lxml_root.find('fork2/t').text = os.path.join(input_directory_prefix,
                                                  run_identifier_prefix +
                                                  'phenotype_vector.txt')
    lxml_root.find('fork3/q').text = os.path.join(input_directory_prefix,
                                                  run_identifier_prefix +
                                                  'structure_matrix.txt')
    lxml_root.find('fork4/k').text = os.path.join(input_directory_prefix,
                                                  run_identifier_prefix +
                                                  'kinship_matrix.txt')

    lxml_root.find('combine6/export').text = os.path.join(
        out_directory_prefix, run_identifier_prefix +'gwas_out_')

    xml_file_name = run_identifier_prefix + 'sim_gwas_pipeline.xml'
    config_file_out = os.path.join(config_file_prefix, run_identifier_prefix
                                   + 'sim_gwas_pipeline.xml')

    lxml_root.write(config_file_out, encoding="UTF-8",
                   method="xml", xml_declaration=True, standalone='',
                    pretty_print=True)


def parameter_set_writer(directory_prefix, run_prefix, mating,
                         quantitative, effects,
                         genetic_structure):
    """
    Simulation parameters are collected in separate dictionary objects.
    This function writes all parameter information into a set of human
    readable .yaml files.

    :param directory_prefix:
    :param run_prefix: Identifier for a set of simulated data
    :param mating: Parameters which specifying mating
    :param quantitative: Dictionary of qtl for each replicate
    :param effects:
    :param genetic_structure:
    :return:
    """

    import yaml


    file_names = {}

    file_names[run_prefix + 'mating.yaml'] = mating
    file_names[run_prefix + 'qtl.yaml'] = quantitative
    file_names[run_prefix + 'allele_effects.yaml'] = effects
    file_names[run_prefix + 'genetic_structure.yaml'] = genetic_structure

    for name, param_set in file_names.items():
        with open(name, 'w') as p_stream:
            yaml.dump(param_set, p_stream)


def parameter_set_reader(parameter_filename):
    """
    Reads a file of .yaml parameters for an easy way to parameterize a
    simulation. Alternately the user would have to derive a great deal of
    information from raw files.
    :param parameter_filename:
    :return:
    """

    pass
