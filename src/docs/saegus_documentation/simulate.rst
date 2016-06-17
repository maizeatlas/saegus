
.. _simulate:

===========
Simulations
===========

.. module:: simulate


.. py:class:: Truncation(generations_of_random_mating=1, generations_of_selection=1, operating_population_size=2000, proportion_of_individuals_saved=0.05, overshoot_as_proportion=0.50, heritability=0.7, meta_pop_sample_sizes=100, number_of_replicates=1, prefounders_file_name='', ae_file_name='')

   :parameter int generations_of_random_mating: Number of generations to perform random mating before selection begins.
   :parameter int generations_of_selection: Number of generations to perform recurrent selection (occurs after random mating)
   :parameter int operating_population_size: Size of the population you want to do selection on. Must be divisible by the number of pairs of founders.

``Truncation`` is a class which encapsulates all parameters and
functions to perform recurrent truncation selection on a quantitative
trait. The trait is assumed to have a simple additive basis so all
phenotypes are calculated by adding contribution of :math:`n` loci plus error.




