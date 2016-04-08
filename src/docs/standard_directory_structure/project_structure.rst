****************************************
Project Skeleton for Using :mod:`saegus`
****************************************

Henceforth I will use a standard directory structure whenever I am working
with ``saegus``. Assuming I am working on the ``devel`` directory a
``saegus`` project will be organized as:


Population Type
===============

<RunIdentifier>
---------------

parameters
~~~~~~~~~~

general
'''''''

   * genotype_source_filename
   * founders
   * operating_population_size
   * number_of_replicates
   * integer_to_snp
   * snp_to_integer

trait
'''''

   * number_of_qtl
   * qtl
   * allele_effect_distribution
   * allele_effect_distribution_parameters
   * allele_effects

genetic_map
'''''''''''

   * loci
   * chromosome
   * cM_position
   * loci_per_chromosome