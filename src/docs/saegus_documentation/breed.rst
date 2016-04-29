============
:mod:`breed`
============


A module which has functions and classes related to customized mating schemes
implemented in simuPOP. ``saegus`` is unique largely because of its support for
complex mating schemes.

.. py:class:: MAGIC(population, recombination_rates)

   .. py:method:: generate_f_one(parental_id_pairs, offspring_per_pair)

   .. py:method:: interim_random_mating(generations_of_random_mating, pop_size)

      :parameter int generations_of_random_mating:
      :parameter int pop_size:


      random_mater = breed.MAGIC(base_population, recombination_rates)
      random_mater.interim_random_mating(3, 2000)



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
