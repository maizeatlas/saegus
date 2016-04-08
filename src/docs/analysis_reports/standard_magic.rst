=========================================
Standardization of Development Population
=========================================

   I have decided to create a standard population to use for all further development.
   Instead of creating a brand new population every single time I change something I
   have created a population called ``standard_magic.pop``. This document explains
   how I created the population. ``standard_magic.pop`` was tested at every single
   generation for an expected pedigree. Below I explain the implementation of
   MAGIC in ``saegus`` and how I tested its correctness.

Creating a Population with Saegus
=================================

   The general scenario in saegus is to begin with a relatively small number of
   individuals called ``prefounders`` and use a customized mating procedure
   to generate a full sized ``operating_population``. ``saegus`` is coded
   on top of the excellent general purpose and highly memory optimized package
   simuPOP_. However, ``simuPOP`` emphasizes quantitative models of disease
   in humans. In saegus we are primarily concerned with quantitative traits.

   .. _simuPOP: http://simupop.sourceforge.net/

Overview of MAGIC in saegus
===========================

   For simplicity I used 8 distinct individuals from our pool of 26 ``prefounders``.
   At minimum every individual in saegus has a unique identifier called ``ind_id``,
   ``mother_id`` and ``father_id``. We can easily examine, store or manipulate full
   pedigrees. The ``standard_magic.pop`` is created from ``founders`` 1 through 8. The
   ``f_one`` population is specified by pairs of individuals as a list. Hence 1 is
   crossed with 2, 3 is crossed with 4 so on and so forth.

   Some parameters we will see throughout this example
   ---------------------------------------------------

   * ``operating_population_size``: size of our population remains fixed starting with the f_one generation.
   * ``pop``: Our workhorse simuPOP.Population to carry out all evolutionary processes.
   * ``ind_id``: A unique floating point number corresponding to an individual. ``simuPOP`` gives us to use  ``ind_id`` to arbitrarily manipulate populations.
   * ``mother_id``: Albeit *Zea mays* can be considered both female and male it is useful at times to distinguish  ``mother`` from ``father``.
   * ``father_id``: Defined for the same reason as ``mother_id``. In general we follow the genetics convention of writing the *female* on the left of the cross.

   :note: All of the parameters are stored in using the Python shelve_ module.

   .. _shelve: https://docs.python.org/3.4/library/shelve.html



   .. code-block::

      founders = [[1, 2], [3, 4], [5, 6], [7, 8]]
       operating_population_size = 2000
       offspring_per_cross = operating_population_size / len(founders)



Implementation of the MAGIC Mating Procedure
============================================

   The major change from previous version of the MAGIC implementation in saegus is how
   I choose which individuals mate at each generation. In previous versions beyond ``f_one``
   I would simply specify which sub-populations *could* mate with one another. However, in
   this version I specify every single mating interaction before the ``evolve`` process.

Creating the F :sub:`1`
~~~~~~~~~~~~~~~~~~~~~~~

   The F :sub:`1` is graphical depiction of the MAGIC mating scheme applied
   to our choice of 8 founders.

      .. image:: reformed_graph.png