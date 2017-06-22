
.. _operators:

===================================
Operators Added By :py:mod:`saegus`
===================================

All of the classes in this module are derived from the simuPOP.PyOperator
class. Each class performs some task inside of the evolutionary process.
simuPOP offers a standard library of operators for common population genetics
processes. The operators defined in this module perform operations which
would either be impossible or difficult to implement using standard simuPOP
operators.

.. _function_forms_of_operators:

.. py:function:: assign_additive_g(pop, qtl, allele_effects)

   :parameter pop: simuPOP.Population
   :parameter qtl: Loci assigned allele effects
   :parameter allele_effects: Dictionary keyed by locus and allele of the effect of an alelle

.. warning::

   :func:`assign_additive_g` assumes that the population has infoField ``g`` defined.


.. py:function:: calculate_error_variance(pop, heritability)

   :parameter pop: simuPOP.Population with a quantitative trait
   :parameter heritability: Parameter determining how much noise exists between genotype and phenotype

.. py:function:: phenotypic_effect_calculator(pop)

   :parameter pop: simuPOP.Population with quantitative trait

   This function only takes into to account additive phenotypic effects for the time being.
   Calculates ``p`` by adding an error term ``epsilon`` to the additive genotypic effect ``g``.
   The error term ``epsilon`` is defined as a random draw from a normal distribution with
   mean :math:`0` and variance :math:`1 - V_g(1/h^2 - 1)`.

   .. note::

      :math:`V_g` is defined as the variance of genotypic effect ``g``.

.. warning::

   :func:`phenotypic_effect_calculator` assumes that the population has infoField ``p`` defined.



.. py:class:: GenoAdditive()

.. py:class:: PhenotypeCalculator()

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
Population must have infoField ``g``, 0 < heritability < 1.
:class:`GenoAdditive` must be called before :class:`CalculateErrorVariance` or
values of ``g`` must be assigned to each individual.

.. py:class:: CullPopulation()

.. py:class:: Sorter()

.. py:class:: MetaPopulation()

.. py:class:: ReplicateMetaPopulation()

.. py:class:: SaveMetaPopulation()

.. py:class:: RandomlyAssignFemaleFitness()

.. py:class:: RandomlyAssignMaleFitness()

.. py:class:: DiscardRandomOffspring()

.. py:class:: SaveMetaPopulations()

.. py:class:: InfoAndGenotypeWriter()
