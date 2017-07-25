
.. _operators_of_saegus:

#############################
Operators of :py:mod:`saegus`
#############################

All of the classes in this module are derived from the
:class:`simuPOP.PyOperator` class. Each class performs some task inside of the
evolutionary process. :py:mod:`simuPOP` offers a standard library of operators
for common population genetics processes. The operators defined in this module
perform operations which would either be impossible or difficult to implement
using standard :py:mod:`simuPOP` operators.

.. _operators:

Operators
#########

.. _geno_additive_array:

:py:class:`GenoAdditiveArray`
=============================

.. py:class:: GenoAdditiveArray(qtl, allele_effects)

   :param list qtl: List of int
   :param numpy.array allele_effects: Array rows of locus, columns by allele state

   .. py:method:: additive_model(pop)

      :param simuPOP.Population: Diploid population with ``g`` defined

   Operator form used during evolutionary scenario to sum the allele effects
   and assign the result to the infoField ``g``.

.. _calculate_error_variance:

:py:class:`CalculateErrorVariance`
==================================

.. py:class:: CalculateErrorVariance(heritability)

   An operator to calculate the variance of the experimental error distribution.
   We assume that there is some degree of error when measuring phenotypes in
   an actual experiment. Measurement error is represented as a random draw
   from a normal distribution with mean zero and variance ``epsilon`` where

.. math::

   e = V_g * (1/h^2 - 1)

``epsilon`` is assigned as a population variable. This operator is typically
called once in the initOps phase of an evolutionary process. At present
:class:`CalculateErrorVariance` is hard coded to calculate
``genotypic_variance`` as the sample variance of the infoField ``g``.
Population must have infoField ``g``, 0 < ``heritability`` < 1.
:class:`GenoAdditive` must be called before :class:`CalculateErrorVariance` or
values of ``g`` must be assigned to each individual.

.. _pheno_additive:

:py:class:`PhenoAdditive`
=========================

.. py:class:: PhenoAdditive()

Mainly an operator to add 'noise' or 'error' to the value from genotypic
effects. Requires ``epsilon`` to be defined as in
py:class:`CalculateErrorVariance`


:py:class:`HDF5AlleleFrequencies`
=================================

.. py:class:: HDF5AlleleFrequencies(allele_frequency_group, allele_data)

   :param h5py.Group allele_frequency_group: group for allele frequencies
   :param numpy.array allele_data: Array of allele states

Operator to store allele frequencies during :py:func:`evolve` process.
``allele_data`` is gathered using :py:func:`analyze.gather_allele_data`. See
the entry in the user guide for collecting and storing data

.. todo:: Show examples of each HDF5 operator

:py:class:`HDF5GenotypeFrequencies`
===================================

.. py:class:: HDF5GenotypeFrequencies(genotype_frequency_group)

   :param h5py.Group genotype_frequency_group: group for genotype frequencies

Operator to store genotype frequencies during :py:func:`evolve` process.
Results are stored in a 3d :py:class:`numpy.array`. The axes are
locus x alpha_allele x omega_allele. The genotypes are interpretted
as coordinates for the purpose of easy storage and access.

:py:class:`HDF5Trait`
=====================

.. py:class:: HDF5Trait(trait_information_field, trait_group)

   :param trait_information_field str: string corresponding to trait
   :param h5py.Group trait group: group for trait information

Operator to store the data from :py:class:`Population.indInfo(trait_information_field)`
such as ``g`` and ``p`` in the User Guide examples. Traits are stored as
generation + '/' trait_information_field.

:py:class:`HDF5Close`
=====================

.. py:class:: HDF5Close(hdf5_file)

Closes the HDF5 file. Meant to be used in ``finalOps`` after an evolutionary
process.


.. py:class:: CullPopulation()

.. py:class:: Sorter()

.. py:class:: MetaPopulation()

.. py:class:: ReplicateMetaPopulation()

.. py:class:: SaveMetaPopulation()

.. py:class:: RandomlyAssignFemaleFitness()

.. py:class:: RandomlyAssignMaleFitness()

.. py:class:: DiscardRandomOffspring()

.. py:class:: SaveMetaPopulations()


.. _function_forms_of_operators:

Function Forms of Operators
###########################

.. _assign_additive_g_function:

:py:func:`assign_additive_g`
============================

.. py:function:: assign_additive_g(pop, qtl, allele_effects)

   :parameter pop: simuPOP.Population
   :parameter qtl: Loci assigned allele effects
   :parameter allele_effects: Dictionary keyed by locus and allele of the effect of an alelle

.. warning::

   :func:`assign_additive_g` assumes that the population has infoField ``g`` defined.

.. _calculate_g:

:py:func:`calculate_g`
======================

.. :py:func:: calculate_g(pop, allele_effects_array)

   :param simuPOP.Population pop: Diploid population with ``g`` defined
   :param allele_effects_array: Array with rows of loci and columns as allele states

.. _calculate_error_variance_function:

:py:func:`calculate_error_variance`
===================================

.. :py:func:: calculate_error_variance(pop, heritability)

   :parameter pop: simuPOP.Population with a quantitative trait
   :parameter heritability: Parameter determining how much noise exists between genotype and phenotype

.. :py:func:: phenotypic_effect_calculator(pop)

   :parameter pop: simuPOP.Population with quantitative trait

   This function only takes into to account additive phenotypic effects for the time being.
   Calculates ``p`` by adding an error term ``epsilon`` to the additive genotypic effect ``g``.
   The error term ``epsilon`` is defined as a random draw from a normal distribution with
   mean :math:`0` and variance :math:`1 - V_g(1/h^2 - 1)`.

   :math:`V_g` is defined as the variance of genotypic effect ``g``.

.. warning::

   :func:`phenotypic_effect_calculator` assumes that the population has infoField ``p`` defined.

.. _calculate_p:

:py:func:`calculate_p`
======================

.. :py:func:: calculate_p(pop)

   Adds error term to each individual's ``g`` value drawn from a normal
   distribution with mean ``0`` and variance as defined in
   :py:func:`calculate_error_variance`