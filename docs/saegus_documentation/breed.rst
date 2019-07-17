
.. _breed_docs:

================
The Breed Module
================

A module which has functions and classes related to customized mating schemes
implemented in simuPOP. ``saegus`` is unique largely because of its support for
complex mating schemes.

Mating Schemes
==============

A list of crosses that I frequently use. Most of these are wrappers around
the ``Population.evolve`` function with some assumed ``Taggers`` arguments.

.. py:class:: SelfCross(recombination_rates)

   .. py:method:: create_top_crosses(population, offspring_per_individual)

      :param population: Whatever

      Creates a population by selfing each individual of ``population`` without replacement.
      In other words each individual self-mates a single time.
      Each individual will generate ``offspring_per_individual`` offspring giving
      a population of ``offspring_per_individual*population`` individuals.

.. code-block:: python
   :caption: Example of create_self_crosses

   >>> pop = sim.Population(100, [10, 10])
   >>> self_crosser = SelfCross([0.01]*20)
   >>> self_crosser.create_self_crosses(pop, 10)
   (1,)
   >>> pop.popSize()
   1000

.. py:class:: MAGIC(population, recombination_rates)

   :parameter population: simuPOP.Population to undergo the MAGIC mating protocol. Also works with simuPOP.Simulator objects.
   :parameter recombination_rates: Probability of recombination between sequential loci.

   .. py:method:: generate_f_one(parental_id_pairs, offspring_per_pair)

      :parameter parental_id_pairs: Nested list of lists. Inner lists are pairs of individuals intended to be crossed.
      :parameter offspring_per_pair: Family size per cross.

      Designed to be used with an initial population of prefounders. In cases
      where more control is required over which individuals are crossed see the
      class :py:class:`SecondOrderPairwiseIDChooser`.

   .. py:method:: random_mating(generations_of_random_mating, pop_size)

      :parameter int generations_of_random_mating: Number of generations to randomly mate the population(s).
      :parameter int pop_size: Population size which can be different or the same as starting population(s).


.. _example_generate_f_one_single_population:

.. code-block:: python
   :caption: Example using a single population

   >>> prefounders = sim.loadPopulation('prefounders1478.pop')
   >>> sim.tagID(prefounders, reset=27)
   >>> magic = breed.MAGIC(prefounders, [0.01]*1478)

   >>> founders = [[1, 2], [3, 4], [5, 6], [7, 8]]
   >>> os_per_pair = 500

   >>> magic.generate_f_one(founders, os_per_pair)
   Generation: 0


.. _example_generate_f_one_replicates:

.. code-block:: python
   :caption: With replicates of a population

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


.. _example_random_mating:

.. code-block:: python
   :caption: example of random mating.

   >>> base_population = simuPOP.loadPopulation('prefounders1478.pop')
   >>> recombination_rates = [0.01] * 1478
   >>> magic = breed.MAGIC(base_population, recombination_rates)
   >>> magic.random_mating(3, 2000)
   Initiating random mating for 3 generations.
   Generation: 0
   Generation: 1
   Generation: 2

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

.. code-block:: python
   :caption: Example of determining a random cross

   >>> mrc = breed.MultiRandomCross(multi_prefounders, 4, 500)
   >>> mothers, fathers = mrc.determine_random_cross()
   >>> mothers[0]
   array([  525.,   482.,   294., ...,  1128.,  1405.,  1297.])
   >>> fathers[0]
   array([  904.,   825.,   751., ...,  1582.,  1911.,  1562.])



Parent Choosers
===============

This is a list of the classes which are utilized to implement customized mating
schemes. A customized mating scheme requires that individuals are passed to the
``OffspringGenerator`` in a specific order. The classes here handle the choosing
and ordering of parents.

.. py:class:: PairwiseIDChooser(pairs_of_parents, offspring_per_parental_pair)

   :parameter pairs_of_parents: Nested lists of list, each sub-list contains a single pair of parents.
   :parameter int offspring_per_parental_pair: Family size

.. py:function:: by_id_pairs(pop)

   :parameter pop: Population to which this parent chooser is applied

.. _second_order_id_chooser:

.. py:class:: SecondOrderPairIDChooser(female_parent_ids, male_parent_ids, offspring_per_parental_pair=1)

   :parameter list female_parent_ids: List of ind_id of individuals chosen as females for mating
   :parameter list male_parent_ids: List of ind_id of individuals chosen as males.
   :parameter int offspring_parental_pair: Family size

.. _multi_second_order_id_chooser:

.. py:class:: MultiSecondOrderPairIDChooser(multi_mother_ids, multi_father_ids, offspring_per_parental_pair=1)

   :parameter multi_mother_ids: Dictionary keyed by replicate of lists of individual IDs (selfing allowed)
   :parameter multi_father_ids: Dictionary keyed by replicate of lists of individual IDs (selfing allowed)
   :parameter offspring_per_parental_pair: The number of offspring to generate per pair of parents.

.. _enforced_population_structure:

.. py:class:: ForcedPopulationStructureParentChooser(expanded_population_size, mating_probabilities)

   :parameter expanded_population_size: The number of offspring to make from the mating of the individuals
   :parameter mating_probabilities: Dictionary keyed by individual ID of probability mass functions implemented in scipy.random

.. _half_sib_bulk_balance:

.. py:class:: HalfSibBulkBalanceChooser(inds_per_breeding_subpop, offspring_per_female)

   :parameter inds_per_breeding_subpop: Number of individuals in a selected breeding group
   :parameter offspring_per_female: Number of times to randomly mate a given female.

   .. py:method:: recursive_pairwise_parent_chooser(pop, subPop)

      :parameter pop: Population to choose parents for
      :parameter subPop: Sub-population to iterate through making the appropriate number of offspring for each breeding group
