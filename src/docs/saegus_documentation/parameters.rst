
.. _parameters:

==========
Parameters
==========


.. py:class:: Trait


   .. _load_alleles:

   .. py:method:: load_alleles(allele_file_name)

      :parameter str allele_file_name: HDF File name containing alleles at each locus

      .. code-block:: python
         :caption: load_alleles_example

         >>> alleles = load_alleles('parameters\\alleles_at_1478_loci.hdf')
         >>> alleles
          array([[1, 2],
          [1, 3],
          [3, 1],
          ...,
          [1, 0],
          [3, 0],
          [3, 1]], dtype=int64)



   .. _assign_ae:

   .. py:method:: assign_allele_effects(alleles, qtl, distribution_function, *distribution_function_args, multiplicity=3)

      :parameter list alleles: np.array or list of lists of alleles at each locus
      :parameter list qtl: loci designated as contributing to a quantitative trait
      :parameter distribution_function: function such as random.expovariate
      :parameter distribution_function_args: arguments necessary for the distribution function
      :parameter int multiplicity: Number of random draws to take from the distribution

   Example: Assigning allele effects to an additive trait.

   .. code-block:: python

      >>> qtl = [85, 94, 378, 417, 431, 730, 935, 1108, 1348, 1355]
      >>> additive_trait = parameters.Trait()
      >>> allele_effects = additive_trait.assign_allele_effects(alleles, qtl, random.expovariate, 1, multiplicity=3)
      {85: {1: 0.7639459962395068, 3: 1.1275557092940487},
      94: {0: 0.8082841215038653, 2: 1.8820116489441723},
      378: {0: 7.048513796426754, 2: 1.4224519757549239},
      417: {1: 6.714168847163591, 3: 1.268012923400879},
      431: {1: 2.6270165938652026, 3: 4.909446892623217},
      730: {1: 3.378195420119303, 3: 3.752044147848409},
      935: {0: 3.1937192305039086, 2: 4.8342880250866},
      1108: {0: 3.214484353047612, 1: 5.40893005938693},
      1348: {1: 5.138900439370714, 3: 4.188077952052308},
      1355: {0: 3.323581565680311, 3: 5.605738561429297}}
