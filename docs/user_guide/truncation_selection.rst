

.. _truncation_selection_guide:

##################################
Truncation Selection in ``saegus``
##################################

This guide will walk through the way of performing truncation or directional
selection. At present the only form of selection implemented in ``saegus`` is
directional selection. Other forms of selection can be implemented with not
too much difficulty. In this case we will perform ten generations of selection.
50 individuals will be sampled at each generation. Selection will be according
to an additive trait given by 30 QTL with an exponential allele effects
distribution.

Setting up the Population
#########################

We will add the necessary information fields and the necessary groups in the
HDF5 file.