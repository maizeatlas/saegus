================
:mod:`operators`
================

All of the classes in this module are derived from the simuPOP.PyOperator
class. Each class performs some task inside of the evolutionary process.
simuPOP offers a standard library of operators for common population genetics
processes. The operators defined in this module perform operations which
would either be impossible or difficult to implement using standard simuPOP
operators.

Classes
=======

:class:`GenoAdditive`


:class:`CalculateErrorVariance`

An operator to calculate the variance of the experimental error distribution.
We assume that there is some degree of error when measuring phenotypes in
an actual experiment. Measurement error is represented as a random draw
from a normal distribution with mean zero and variance ``epsilon`` where

   ``epsilon`` = ``genotypic_variance`` * (1/``heritability`` - 1)

``epsilon`` is assigned as a population variable. This operator is typically
called once in the initOps phase of an evolutionary process. At present
:class:`CalculateErrorVariance` is hard coded to calculate
``genotypic_variance`` as the sample variance of the infoField ``g``.

Requires
--------

Population must have infoField ``g``, 0 < heritability < 1.
:class:`GenoAdditive` must be called before class:`CalculateErrorVariance` or
values of ``g`` must be assigned to each individual.

:param float heritability: Floating point value greater than zero and less than 1.