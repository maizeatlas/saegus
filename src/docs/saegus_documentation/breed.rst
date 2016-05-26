================
The Breed Module
================


A module which has functions and classes related to customized mating schemes
implemented in simuPOP. ``saegus`` is unique largely because of its support for
complex mating schemes.

.. py:class:: MAGIC(population, recombination_rates)

   :parameter population: simuPOP.Population to undergo the MAGIC mating protocol. Also works with simuPOP.Simulator objects.
   :parameter recombination_rates: Probability of recombination between sequential loci.

   .. py:method:: generate_f_one(parental_id_pairs, offspring_per_pair)

      :parameter parental_id_pairs: Nested list of lists. Inner lists are pairs of individuals intended to be crossed.
      :parameter offspring_per_pair: Family size per cross.

      Designed to be used with an initial population of prefounders. In cases
      where more control is required over which individuals are crossed see the
      class :py:class:`SecondOrderPairwiseIDChooser`.

      .. _example_generate_f_one_single_population:

      **Example**: With a single population:

      .. code-block:: python

         >>> prefounders = sim.loadPopulation('prefounders1478.pop')
         >>> sim.tagID(prefounders, reset=27)
         >>> magic = breed.MAGIC(prefounders, [0.01]*1478)

         >>> founders = [[1, 2], [3, 4], [5, 6], [7, 8]]
         >>> os_per_pair = 500

         >>> magic.generate_f_one(founders, os_per_pair)
         Generation: 0

      .. _example_generate_f_one_replicates:

      **Example**: With replicates of a population:

      .. code-block:: python

         >>> prefounders = sim.loadPopulation('prefounders1478.pop')
         >>> sim.tagID(prefounders, reset=27)
         >>> multi_prefounders = simuPOP.Simulator(prefounders, 5, stealPops=False)
         >>> multi_magic = breed.MAGIC(multi_prefounders, [0.01]*1478)

         >>> founders = [[1, 2], [3, 4], [5, 6], [7, 8]]
         >>> os_per_pair = 500

         >>> multi_magic.generate_f_one(founders, os_per_pair)
         Generation: 0
         Generation: 0
         Generation: 0
         Generation: 0
         Generation: 0

   .. py:method:: random_mating(generations_of_random_mating, pop_size)

      :parameter int generations_of_random_mating: Number of generations to randomly mate the population(s).
      :parameter int pop_size: Population size which can be different or the same as starting population(s).

      .. _example_random_mating:

      **Example**: Random mating with some tagging information.

      .. code-block:: python

         >>> base_population = simuPOP.loadPopulation('prefounders1478.pop')
         >>> recombination_rates = [0.01] * 1478
         >>> magic = breed.MAGIC(base_population, recombination_rates)
         >>> magic.random_mating(3, 2000)
         Initiating random mating for 3 generations.
         Generation: 0
         Generation: 1
         Generation: 2



MAGIC: Multi-parent Advanced Generation Inter-crosses

MAGIC begins with ``founders`` arranged into pairs. Pairs of parents
are crossed with each other to make hybrid offspring. The hybrid offspring
are crossed with a different group of hybrid offspring to make
double-hybrid offspring. This process continues until only a single
sub-population of individuals remains. Given the ``founders`` MAGIC
merges successive generations. The process is demonstrated visually with
a graph.

   **Example**

   .. code-block:: python

      founders = [[1, 2], [3, 4], [5, 6], [7, 8]]

   .. figure:: reformed_graph.png
      :scale: 50%

      Successive rounds of mating merges together the genomes of distinct founders.



   The result of these three rounds of crossing is a population with high
   genetic diversity. High genetic diversity stemming from the combination of
   many different parental lines greatly enhances mapping resolution.

   .. class::(pop, recombination_rates)

   :param pop: A simuPOP.Population subjected to MAGIC mating scheme.
   :param recombination_rates: List of recombination rates at each locus.
   :note: Recombination rates refer to the simuPOP definition.


Predicting Number of Rounds of Mating
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given :math:`n` founders to start with:

If :math:`n` is the :math:`m` a power of two, :math:`m`, then there are *m* rounds of mating.
That is if:

.. math::

   2^m = n

The number of rounds of mating is the largest power of :math:`2` which can
be subtracted from :math:`n` and then one more generation to combine the final pair
of individuals.

.. _pairwise_id_chooser:

.. py:class:: PairwiseIDChooser(pairs_of_parents, offspring_per_pair=1)

   :parameter pairs_of_parents: A list of lists or tuples. Contents of inner lists or tuples are a pair of IDs
   :parameter int offspring_per_pair: Number of offspring to generate from each pair of individuals

   .. py:method:: by_id_pairs(population)

      :parameter population: simuPOP.Population with infoField ``ind_id`` defined.

      .. _example_pairwise_id_chooser:

         **Example** Nested list of prefounder IDs




.. _second_order_id_chooser:

.. py:class:: SecondOrderPairIDChooser(female_parent_ids, male_parent_ids, offspring_per_parental_pair=1)

   :parameter list female_parent_ids: List of ind_id of individuals chosen as females for mating
   :parameter list male_parent_ids: List of ind_id of individuals chosen as males.
   :parameter int offspring_parental_pair: Family size

.. _multi_random_cross:

.. py:class:: MultiRandomCross(multi_replicate_population, number_sub_pops, sub_pop_size)

   A class intended for use with simulations using multiple replicates which
   require predictable mating among sub-populations of each replicate.

   :parameter multi_replicate_population: simuPOP.Simulator
   :parameter int number_sub_pops: Each replicate is split into ``number_sub_pops`` sub-populations. Parent IDs are sampled from each sub-population.
   :parameter sub_pop_size: Size of each sub-population.

   .. py:method:: determine_random_cross()

      Creates separate dictionaries for IDs of mothers and fathers respectively.
      Entries are keyed corresponding to ``rep`` of the replicate the
      IDs are taken from.

      .. _example_determine_random_cross:

      **Example**: Creating arrays of randomly chosen individuals for each replicate.

      .. code-block:: python

      >>> mrc = breed.MultiRandomCross(multi_prefounders, 4, 500)
      >>> mothers, fathers = mrc.determine_random_cross()
      >>> mothers[0]
      array([  525.,   482.,   294., ...,  1128.,  1405.,  1297.])
      >>> fathers[0]
      array([  904.,   825.,   751., ...,  1582.,  1911.,  1562.])


