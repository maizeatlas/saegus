.. _analysis_using_tassel:

#####################################################
Using Standalone Tassel for Multi-Locus Trait Mapping
#####################################################

This example will show how to perform association mapping of a quantitative
trait using TASSEL-5-standalone. TASSEL requires a set of input for its mixed
linear model.

* Kinship Matrix
* Trait

:py:mod:`saegus` has functions to calculate
the (kinship) relationship matrix using marker data [VanRaden2008]_. Population
structure is computed using [Patterson2006]_.





.. [Patterson2006] Patterson, N, Price, A, Reich, D. (2006). Population Structure and Eigenanalysis

.. [VanRaden2008] VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. Journal of Dairy Science, 91(11), 4414–23. doi:10.3168/jds.2007-0980
