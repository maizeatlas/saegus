NAM Genotype Data
=================

This page documents the process of extracting the information from
*hapmap3.txt* to make the nam_prefounders.pop file.

   Parse the *hapmap3.txt* file using a function from the pandas_ package.
   .. _pandas: http://pandas.pydata.org/
   The genotypes are given in columns of string characters. The genetic map,
    genotypes and names of each founder are all contained within this file.
    Recombination rates are tailored to fit the unique triplet scheme.

    =====

