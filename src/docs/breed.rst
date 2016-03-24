============
:mod:`breed`
============


A module which has functions and classes related to customized mating schemes
implemented in simuPOP. ``saegus`` is unique largely because of its support for
complex mating schemes.



:class:`MAGIC`
==============

   MAGIC: Multi-parent Advanced Generation Inter-crosses
   MAGIC begins with founders arbitrarily arranged into pairs. Pairs of parents
   are crossed with each other to make hybrid offspring. The hybrid offspring
   are crossed with a different group of hybrid offspring to make
   double-hybrid offspring. This process continues until only a single
   sub-population of individuals remains.

   **Example**

      founders = [[1, 2], [3, 4], [5, 6], [7, 8]]

   First Round of Crosses
   ----------------------

      1 x 2 | 3 x 4 | 5 x 6 | 7 x 8
       1/2  |  3/4  |  5/6  | 7/8

   Second Round of Crosses
   -----------------------

      1/2 x 3/4 | 5/6 x 7/8
       1/2//3/4 | 5/6//7/8

   Third (Final) Round of Crosses
   ------------------------------

      1/2 // 3/4 x 5/6 // 7/8
       1/2//3/4 /// 5/6//7/8

   The result of these three rounds of crossing is a population with high
   genetic diversity. High genetic diversity stemming from the combination of
   many different parental lines greatly enhances mapping resolution.

   .. py:method:: __init__(pop, recombination_rates)

   :param pop: A simuPOP.Population subjected to MAGIC mating scheme.

   :param recombination_rates: List of recombination rates at each locus.

      :note: Recombination rates refer to the simuPOP definition.