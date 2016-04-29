=================
:mod:`parameters`
=================


.. py:class:: Trait

   .. _assign_ae::

   .. py:method:: assign_allele_effects(alleles, qtl, distribution_function, *distribution_function_args, multiplicity=3)

      :parameter list alleles:
      :parameter list qtl:
      :parameter distribution_function:
      :parameter *distribution_function_args:
      :parameter int multiplicity:

   **Example**: *Assigning allele effects to an additive trait.*

   additive_trait = parameters.Trait()
   allele_effects = additive_trait.assign_allele_effects(alleles, qtl, random.expovariate, 1, multiplicity=3)