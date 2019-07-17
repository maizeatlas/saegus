
.. _parameters:

==========
Parameters
==========


.. _population_structure:

Population Structure
====================


.. py:class:: PopulationStructure(pop)

   :parameter pop: simuPOP.Population

   This class was developed using the Tuson population and its corresponding raw
   data. The Tuson prefounder population structure was estimated from Tuson
   G\:sub:`0` sample of 105 individuals. To recapitulate that population structure
   in the simulated populations we use the estimated proportions of inheritance
   for each individual. For example:


   sample_id |    1      |    2     |    3    |    4    |    5   |   6
   1         |    0.0    |   0.96   |  0.0    |    0.0  |  0.04  |   0.0
   2    	    |    0.0110 |   0.0015 |  0.0004 | 0.1047  |  0.0   |  0.8824

   Individual 1 derives 96 percent of its genome from prefounder population 2 and
   4 percent of its genome from prefouner population 5. Individual 2 derives
   1 percent of its genome from prefounder population 1, 0.15 percent from population
   2, 0.04 percent from 3, 10 percent from 4 and 88 percent from 6. Individual 1's
   primary population is 2 and individual 2's primary population is 6. This is
   assigned to the ``primary`` infoField. Hence if individual 1 is picked to
   for a mating event the probability it will mate with an individual from
   prefounder population 2 is 0.96. The propability it will mate with an
   individual from population 5 is 0.04. So on and so forth for individual 2.

   .. py:method:: parse_and_remap_population_structure(population_structure_matrix_file_name)

      :parameter str population_structure_matrix_file_name: Input pop strct file

      Parses a population structure matrix from a file and converts it into
      Python dictionary. The :file:`population_structure_matrix.xlsx` file is
      not in the same order as :file:`genotype_matrix.txt`. This function
      remaps the inheritance proportions to match the genotype matrix of the
      population.

   .. py:method:: generate_mating_pmfs(population_structure_dict)

      :parameter population_structure_dict: Dictionary of lists keyed by integer corresponding to ``ind_id`` infoField

      Converts a dictionary of lists of probabilities into scipy.stats.rv_discrete
      customized probability mass functions.

   .. py:method:: assign_primary_subpopulation(pop, struct_mating_probabilities)

      :param struct_mating_probabilities: Dict of lists of probabilities

      Assigns the primary subpopulation to each individual according to
      ``ind_id``. Primary subpopulation is the population from which the
      individual derives most of its genome.

   .. code-block:: py
      :caption: Example use-case of :py:class:`PopulationStructure`

      >>> tuson = sim.loadPopulation('artemis_tuson.pop')
      >>> artemis_popst = parameters.PopulationStructure(tuson)
      >>> struct_mating_probs = artemis_popst.parse_and_remap('population_structure_matrix.xlsx')
      {
      1: [0.0, 0.99960000000000004, 0.0, 0.00040000000000000002, 0.0, 0.0],
      2: [0.010999999999999999, 0.0015, 0.00040000000000000002, 0.1047, 0.0,
         0.88239999999999996],
      3: [0.0, 0.38319999999999999, 0.0, 0.61680000000000001, 0.0, 0.0],
      ...,
      }
      >>> mating_pmfs = artemis_popst.generate_mating_pmfs(struct_mating_probs)
      >>> mating_pmfs
      {
      1: <scipy.stats._distn_infrastructure.rv_discrete at 0xf5a1358>,
      2: <scipy.stats._distn_infrastructure.rv_discrete at 0xf64da58>,
      3: <scipy.stats._distn_infrastructure.rv_discrete at 0xf595c88>,
      ...,
      }
      >>> artemis_popst.assign_primary_subpopulations(struct_mating_probabilities)


.. _missing_genotype_data:

.. py:class:: MissingGenotypeData(raw_genotype_array, number_of_individuals, number_of_markers, missing_genotype_token)

   :parameter raw_genotype_array: Array of genotypes
   :parameter int number_of_individuals: Individuals in a population
   :parameter int number_of_markers:
   :parameter str missing_genotype_token: Missing genotype token i.e. "./." or "NA"

   .. py:method:: genotype_counter()

   .. py:method:: determine_individuals_missing_data(loci_missing_data, missing_genotype_token)

      :parameter loci_missing_data:
      :parameter str missing_genotype_token:

   .. py:method:: convert_counts_to_frq(genotype_counts):

      :parameter genotype_counts:

   .. py:method:: generate_pmf_mappings(empirical_pmf_mappings)

      :parameter empirical_pmf_mappings:

   .. py:method:: replace_missing_genotypes()


.. _trait:

Trait
-----

.. py:class:: Trait


   .. py:method:: construct_allele_effects_table(pop, qtl, distribution_function, *distribution_function_parameters)








   .. _load_alleles:

   .. py:method:: load_alleles(allele_file_name)

      :parameter str allele_file_name: HDF File name containing alleles at each locus

      .. code-block:: python
         :caption: load_alleles_example

         >>> alleles = load_alleles('parameters\\alleles_at_1478_loci.hdf')
         >>> alleles
          array([[1, 2],
          [1, 3],
          [3, 1],
          ...,
          [1, 0],
          [3, 0],
          [3, 1]], dtype=int64)



   .. _assign_ae:

   .. py:method:: assign_allele_effects(alleles, qtl, distribution_function, *distribution_function_args, multiplicity=3)

      :parameter list alleles: np.array or list of lists of alleles at each locus
      :parameter list qtl: loci designated as contributing to a quantitative trait
      :parameter distribution_function: function such as random.expovariate
      :parameter distribution_function_args: arguments necessary for the distribution function
      :parameter int multiplicity: Number of random draws to take from the distribution


.. code-block:: python
   :caption: Assigning allele effects to a population

   >>> qtl = [85, 94, 378, 417, 431, 730, 935, 1108, 1348, 1355]
   >>> additive_trait = parameters.Trait()
   >>> allele_effects = additive_trait.assign_allele_effects(alleles, qtl, random.expovariate, 1, multiplicity=3)
   {85: {1: 0.7639459962395068, 3: 1.1275557092940487},
   94: {0: 0.8082841215038653, 2: 1.8820116489441723},
   378: {0: 7.048513796426754, 2: 1.4224519757549239},
   417: {1: 6.714168847163591, 3: 1.268012923400879},
   431: {1: 2.6270165938652026, 3: 4.909446892623217},
   730: {1: 3.378195420119303, 3: 3.752044147848409},
   935: {0: 3.1937192305039086, 2: 4.8342880250866},
   1108: {0: 3.214484353047612, 1: 5.40893005938693},
   1348: {1: 5.138900439370714, 3: 4.188077952052308},
   1355: {0: 3.323581565680311, 3: 5.605738561429297}}

