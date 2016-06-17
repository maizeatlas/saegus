
.. _simulate:

===========
Simulations
===========

.. module:: simulate


.. py:class:: Truncation(generations_of_random_mating=1, generations_of_selection=1, operating_population_size=2000, proportion_of_individuals_saved=0.05, overshoot_as_proportion=0.50, heritability=0.7, meta_pop_sample_sizes=100, number_of_replicates=1, prefounders_file_name='', ae_file_name='')

   :parameter int generations_of_random_mating: Number of generations to perform random mating before selection begins.
   :parameter int generations_of_selection: Number of generations to perform recurrent selection (occurs after random mating)
   :parameter int operating_population_size: Size of the population you want to do selection on. Must be divisible by the number of pairs of founders.
   :parameter float proportion_of_individuals_saved: The percentage of the population at the extreme to save. Currently hardcoded for single direction truncation selection. Used as a proxy for selection intensity
   :parameter float overshoot_as_proportion: The multiple of individuals to create over and beyond :py:data:`operating_population_size`. The operating population is sampled from the overshoot.
   :parameter float heritability: Narrow sense heritability. Used to define the amount of error in the phenotype calculations
   :parameter meta_pop_sample_sizes: Either a list or a single number. A single number implies a single sample is taken. Multiple numbers imply that multiple samples are taken.
   :paramter str prefounders_file_name: File name to load the prefounder .pop file.
   :parameter str ae_file_name: File name to load allele effects from.

   ``Truncation`` is a class which encapsulates all parameters and
   functions to perform recurrent truncation selection on a quantitative
   trait. The trait is assumed to have a simple additive basis so all
   phenotypes are calculated by adding contribution of :math:`n` loci plus error.

   .. py:method:: recombinatorial_convergence(pop, recombination_rates):

      :parameter pop: simuPOP.Population object
      :parameter recombination_rates: A list of probabilities of recombination between sequential loci

   .. py:method:: expand_by_selfing(pop, recombination_rates)

      :parameter pop: simuPOP.Population object
      :parameter recombination_rates: A list of probabilities of recombination between sequential loci

   .. py:method:: interim_random_mating(pop, recombination_rates)

      :parameter pop: simuPOP.Population object
      :parameter recombination_rates: A list of probabilities of recombination between sequential loci

   .. py:method:: recurrent_truncation_selection(pop, recombination_rates)

      :parameter pop: simuPOP.Population object
      :parameter recombination_rates: A list of probabilities of recombination between sequential loci

   .. py:method:: replicate_selection(pop, recombination_rates)

      :parameter pop: simuPOP.Population object
      :parameter recombination_rates: A list of probabilities of recombination between sequential loci
