:mod:`operators` --- Operators: Functions Applied During Evolutionary Processes
===============================================================================



 All of the classes in this module are derived from the simuPOP.PyOperator
 class. Each class performs some task inside of the evolutionary process.
 simuPOP offers a standard library of operators for common population genetics
 processes. These operators are defined for the purpose of investigating
 the genetic and statistical properties of recurrently selected populations.



class CalculateErrorVariance
----------------------------


class GenotypicEffectCalculator
-------------------------------

.. class:: GenotypicEffectCalculator

   In our simulations an individual's genotype determines phenotype. We are
   primarily interested in investigating quantitative trait loci (QTL). The
   user supplies a set of QTL in the form of a list of integers. Each integer
   corresponds to a locus. The effect for each allele is a random draw
   from a common distribution. At present each of the effects is drawn from
   an exponential distribution.

   .. method:: GenotypicEffectCalculator(absolute_index_qtl, allele_effects,
begin=0, end=-1, step=1, at=[], subPops=ALL_AVAIL, infoFields=[])

      Creates an operator which simply sums the contribution of each locus by
      examining the genotype at each QTL and retrieving the allele effect
      from ``allele_effects``. Parameter ``absolute_index_qtl`` specifies
      which loci are considered QTL. As the name suggests the user must
      define each locus by its absolute index.

    .. function:: genotypic_contribution_calculator(self, ind)




class PhenotypeCalculator
--------------------------------

.. class:: PhenotypicEffectCalculator

   Our simulation attempts to reflect reality as closely as possible.
   Therefore, to reflect error which occurs in actual experiments
   an error term ``epsilon`` is added to an individual's genotypic
   contribution ``g`` to yield phenotype ``p``. ``proportion_selected`` refers
   to the individuals which are in the top ``proportion_selected`` by genotype
   of the population. The population is sorted in descending order by ``g``
   and phenotypes are calculated

   .. method:: PhenotypicEffectCalculator(proportion_selected, begin=0, end=-1,
step=1, at=[], subPops=ALL_AVAIL, infoFields=[])

      An operator to calculate the phenotype of an individual by adding an
      error term to the contribution from genotype. The error term ``epsilon``
      is a random draw from a normal distribution with mean 0 and variance
      variance(1/heritability - 1).
