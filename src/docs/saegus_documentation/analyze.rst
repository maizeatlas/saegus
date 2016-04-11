==============
:mod:`analyze`
==============



.. py:function:: allele_data(pop, alleles, loci)

   Determines the minor alleles, minor allele frequencies, major alleles and
   major allele frequencies.


   :parameter pop: Population intended for GWAS analysis
   :parameter list loci: All loci of a population or a subset
   :parameter dict alleles: Dictionary of alleles present at each locus

   This function is used to find the major/minor alleles of a Population
   ``pop`` given a list of ``alleles`` at each locus given in ``loci``.
   The output is intended to be used in other functions to determine the
   kinship matrix and population structure.

   Additionally this function will also resolve ties between the
   major and minor alleles which result when both alleles have exactly equal
   frequency i.e. 0.50.


.. code-block:: python

   pop = sim.loadPopulation('magic1478.pop')
   loci = list(range(pop.totNumLoci()))
   alleles = shelve.open('magic_1478_simulation_parameters')
   alleles

   {0: [1, 2],
    1: [1, 3],
    2: [3, 1],
    3: [0, 2],
    4: [2, 0],
    5: [0, 2],
    6: [0, 2],
    7: [3, 1],
    8: [0, 2],
    ...}

    af = analyze.allele_data(magic1478_2718, alleles, list(range(1478)))
    af

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>minor_allele</th>
          <th>minor_frequency</th>
          <th>major_allele</th>
          <th>major_frequency</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>2</td>
          <td>0.00000</td>
          <td>1</td>
          <td>1.00000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>3</td>
          <td>0.13275</td>
          <td>1</td>
          <td>0.86725</td>
        </tr>
        <tr>
          <th>2</th>
          <td>1</td>
          <td>0.06575</td>
          <td>3</td>
          <td>0.93425</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2</td>
          <td>0.00000</td>
          <td>0</td>
          <td>1.00000</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0</td>
          <td>0.05675</td>
          <td>2</td>
          <td>0.94325</td>
        </tr>
        <tr>
          <th>5</th>
          <td>2</td>
          <td>0.24875</td>
          <td>0</td>
          <td>0.75125</td>
        </tr>
        <tr>
          <th>6</th>
          <td>2</td>
          <td>0.12300</td>
          <td>0</td>
          <td>0.87700</td>
        </tr>
        <tr>
          <th>7</th>
          <td>1</td>
          <td>0.00000</td>
          <td>3</td>
          <td>1.00000</td>
        </tr>
        <tr>
          <th>8</th>
          <td>2</td>
          <td>0.24000</td>
          <td>0</td>
          <td>0.76000</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
      </tbody>
    </table>
    <p>1478 rows × 4 columns</p>
    </div>

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
   :param genetic_map: Chromosome:cM position correspondence.

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


.. py:function:: generate_allele_effects_table(qtl, founder_alleles, allele_effects)

    Creates a simple pd.DataFrame for allele effects. Hard-coded
    for bi-allelic case.

    :parameter list qtl: List of loci declared as QTL
    :parameter np.array alleles: Array of alleles at each locus
    :parameter dict allele_effects: Mapping of effects for alleles at each QTLocus

**Example**

.. code-block:: python

   alleles
   array([[1, 2],
        [1, 3],
        [3, 1],
        ...,
        [1, 0],
        [3, 0],
        [3, 1]], dtype=int64)

   qtl
   [44, 103, 168, 340, 488, 639, 737, 819, 981, 1065]

   allele_effects
   {44: {0: 5.629446187924926, 2: 1.8962727055819322},
   103: {0: 1.3097813991257303, 2: 6.14070564290979},
   168: {2: 6.718096248082958, 3: 4.697238579652859},
   340: {1: 1.521689147484636, 2: 2.2131077852927032},
   488: {1: 2.512286137462885, 3: 2.486777318327935},
   639: {0: 1.1268072986309254, 3: 1.3391282487711016},
   737: {0: 1.4879865577936147, 1: 1.607534785598338},
   819: {1: 2.2153417608326986, 3: 0.20077940947200731},
   981: {0: 3.9513501430851568, 3: 1.78843909724396},
   1065: {0: 0.998194377898828, 2: 1.5139052352904945}}

    aeframe

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>locus</th>
          <th>alpha_allele</th>
          <th>alpha_effect</th>
          <th>beta_allele</th>
          <th>beta_effect</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>44</td>
          <td>0</td>
          <td>5.629446</td>
          <td>2</td>
          <td>1.896273</td>
        </tr>
        <tr>
          <th>1</th>
          <td>103</td>
          <td>0</td>
          <td>1.309781</td>
          <td>2</td>
          <td>6.140706</td>
        </tr>
        <tr>
          <th>2</th>
          <td>168</td>
          <td>2</td>
          <td>6.718096</td>
          <td>3</td>
          <td>4.697239</td>
        </tr>
        <tr>
          <th>3</th>
          <td>340</td>
          <td>2</td>
          <td>2.213108</td>
          <td>1</td>
          <td>1.521689</td>
        </tr>
        <tr>
          <th>4</th>
          <td>488</td>
          <td>3</td>
          <td>2.486777</td>
          <td>1</td>
          <td>2.512286</td>
        </tr>
        <tr>
          <th>5</th>
          <td>639</td>
          <td>0</td>
          <td>1.126807</td>
          <td>3</td>
          <td>1.339128</td>
        </tr>
        <tr>
          <th>6</th>
          <td>737</td>
          <td>1</td>
          <td>1.607535</td>
          <td>0</td>
          <td>1.487987</td>
        </tr>
        <tr>
          <th>7</th>
          <td>819</td>
          <td>1</td>
          <td>2.215342</td>
          <td>3</td>
          <td>0.200779</td>
        </tr>
        <tr>
          <th>8</th>
          <td>981</td>
          <td>0</td>
          <td>3.951350</td>
          <td>3</td>
          <td>1.788439</td>
        </tr>
        <tr>
          <th>9</th>
          <td>1065</td>
          <td>2</td>
          <td>1.513905</td>
          <td>0</td>
          <td>0.998194</td>
        </tr>
      </tbody>
    </table>
    </div>





.. py:function:: collect_haplotype_data(pop, allele_effects, quantitative_trait_loci)

    :param pop:
    :param allele_effects:
    :param quantitative_trait_loci:

.. py:function:: generate_haplotype_data_table(pop, haplotype_data)

    Generates a table for easy analysis and visualization of haplotypes,
    effects, frequencies and locations.

    :param pop:
    :param haplotype_data:

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

.. py:class:: GWAS(pop, loci, allele_subset, run_id, individual_names, locus_names, pos_names)

   A class to collect all the functions necessary to prepare viable input for TASSEL's
   mixed linear model method. Most of the methods contribute to generating one
   of the files required to run MLM in TASSEL.

   :parameter pop: Population under analysis. Can be full sized or sample.
   :parameter loci: Either all loci in a population or a subset. Usually the segregating loci.
   :parameter allele_subset: A 1D list of alleles at each locus. Typically a list of minor alleles.
   :parameter run_id: A string prefixed to files and various objects to track the corresponding meta-data
   :parameter individual_names: Stringified version of the individual IDs of :param:`pop`. Can be used to attach meta-data to each ID
   :parameter locus_names: List of loci indices for the ``hapmap`` required by TASSEL.
   :parameter pos_names: Dummy column for our purposes. However, it is properly used to represent physical position on a genome.


   .. py:method:: hapmap_formatter(self, int_to_snp_conversions, hapmap_matrix_filename)

      :parameter dict int_to_snp_conversions: Simple dictionary which encodes integer alleles as strings
      :parameter hapmap_matrix_filename: Output filename

      Converts genotype data from sim.Population object to HapMap file format
      in expectation to be used in TASSEL for GWAS. At present the column
      names will be hardcoded as will some of the values. ``hapmap_matrix_filename``
      is the name of the file the formatted will be written to.

   .. py:method:: trait_formatter(self, trait_filename)

      :parameter str trait_filename: Output for the ID-phenotype vector

      Simple function to automate the formatting of the phenotype data.
      TASSEL requires a specific header in the trait file. This function
      simply write the header first and then adds the rest of the data in using
      :function:`np.savetxt`

      **Example**

      .. code-block:: python

         header = "<Trait> sim\n"


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

