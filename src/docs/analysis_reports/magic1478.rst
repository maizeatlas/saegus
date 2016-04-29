======================
MAGIC1478 vs MAGIC7386
======================

The ``standard_magic`` population uses all 7386 loci from :file:`hapmap3.txt`. However,
if I use all 7386 loci and the triplet scheme this alters the GWAS results. I need
to create a population which uses the loci which occupy integer-valued positions on the
genetic map i.e. 1.0 cM, 2.0 cM. ``standard_magic`` has the loci at 0.6 cM, 0.8 cM, 1.0 cM, 1.2 cM
so on and so forth.

Allele Effects and Triplets of Non-Recombining Loci
===================================================

The ``standard_magic`` population assigns effects to an integer-valued position
and the positions immediately upstream and downstream. When I use the subset of 1478 loci
I will assign the alleles at those loci effects as **three independent** draws from
an exponential distribution. The goal is to make ``magic1478`` somewhat comparable to ``standard_magic``.

Strategy of Creating MAGIC1478
==============================

I am going to use the ``nam_prefounders`` and withdraw the integer-valued positions of each individual.
This saves me the work of building up the population from the raw files.


Creating the simuPOP Population Object for MAGIC1478
----------------------------------------------------

A simuPOP ``Population`` requires the number of chromosomes and the number of loci per chromosome.
Luckily I saved the absolute indexes of the integer-valued positions of the genetic map in a
:mod:`shelve` instance.

.. code-block:: python

   misc_gmap = shelve.open('misc_gmap')
   integral_valued_loci = misc_gmap['integral_valued_loci']
   relative_integral_valued_loci = misc_gmap['relative_integral_valued_loci']
   print(integral_valued_loci[:10])
   [4, 9, 14, 19, 24, 29, 34, 39, 44, 49]

``relative_integral_valued_loci`` is a ``dict`` which is keyed by absolute locus index. We will extract the chromosome
of each absolute valued index and then use ``collections.Counter`` to count how many loci are on each chromosome.

.. code-block:: python


   print(relative_integral_valued_loci[4])
   (1, -4.0)

   integer_loci_per_chromosome = [relative_integral_valued_loci[abs_locus][0] for abs_locus in integral_valued_loci]

   import collections as col

   integer_valued_loci_per_chromosome = col.Counter(loci_per_chromosome)
   print(integer_valued_loci_per_chromosome)

   Counter({1.0: 210, 3.0: 164, 2.0: 161, 5.0: 157, 4.0: 152, 7.0: 139, 8.0: 138, 9.0: 132, 10.0: 113, 6.0: 112})

   loci_per_chromosome = tuple(integer_valued_loci_per_chromosome.values())

   print(loci_per_chromosome)
   (210, 161, 164, 152, 157, 112, 139, 138, 132, 113)

So then I can simply do:

.. code-block:: python

   prefounders_1478 = sim.Population(size=26, ploidy=2, loci=loci_per_chromosome)

This initializes the population and now I can copy over genotypes at the integer valued positions from
the prefounders7386 population.

.. code-block:: python

   prefounders7386 = sim.loadPopulation('nam_prefounders.pop')
   for ind_7386, ind_1478 in zip(prefounders_7386.individuals(), prefounders_1478.individuals()):
      sub_alpha = [ind_7386.genotype(ploidy=0)[integer_position]
                  for integer_position in integer_subset]
      sub_beta = [ind_7386.genotype(ploidy=1)[integer_position]
                 for integer_position in integer_subset]
      genotype_1478 = sub_alpha + sub_beta
      ind_1478.setGenotype(genotype_1478)

.. _making-prefounders-1478:

Generating Standard MAGIC1478 from Prefounders1478
==================================================

I followed the same mating and testing strategy to make a ``standard_magic_1478.pop`` file. I used
the same founders as in ``standard_magic`` i.e. founders 1 through 8.

.. code-block:: python

   prefounders_1478 = sim.loadPopulation('prefounders_1478.pop')


Determine the Mating Pairs of Each Generation
---------------------------------------------

I created a standardized MAGIC1478 population as I did with MAGIC7386. At each step
I pre-determine the mating pairs and record them in ``lists`` which have the title
``expected_x_mother_ids`` and ``expected_x_father_ids``. The expected parental id pairs are
mated in order. The offspring have infoFields which record the ID of their mother and ID father.


After mating the offspring parental IDs are compared with the expected parental IDs.
Below is an example of this mating-testing cycle.

.. code-block:: python

   first_sp_mothers = [random.choice(pop.indInfo('ind_id', 0)) for i in range(1000)]
   second_sp_fathers = [random.choice(pop.indInfo('ind_id', 1)) for i in range(1000)]
   third_sp_mothers = [random.choice(pop.indInfo('ind_id', 2)) for i in range(1000)]
   fourth_sp_fathers = [random.choice(pop.indInfo('ind_id', 3)) for i in range(1000)]

   expected_f_two_mother_ids = first_sp_mothers + third_sp_mothers
   expected_f_two_father_ids = second_sp_fathers + fourth_sp_fathers

The expected parental IDs are written to disk using a ``shelve`` for post-comparison should it be necessary.

.. code-block:: python

   breeding_parameters['expected_f_two_mother_ids'] = expected_f_two_mother_ids
   breeding_parameters['expected_f_two_father_ids'] = expected_f_two_father_ids

A ``parent_chooser`` is initiated which determines how offspring are created.

.. code-block:: python

   second_order_pc = breed.SecondOrderPairIDChooser(expected_f_two_mother_ids, expected_f_two_father_ids)

Then mating occurs:

.. code-block:: python

   pop.evolve(
       matingScheme=sim.HomoMating(
           sim.PyParentsChooser(second_order_pc.snd_ord_id_pairs),
           sim.OffspringGenerator(ops=[
               sim.IdTagger(),
               sim.ParentsTagger(),
               sim.PedigreeTagger(),
               sim.Recombinator(rates=0.01)
           ],
               numOffspring=1),
           subPopSize=[2000],
       ),
       gen=1,
   )

We organize mother and father IDs into ``observed`` lists and compare them to the expected IDs.
We count the number of matches between expected and observed mother IDs and expected and observed father IDs.
The number should be equal to the population size.


.. code-block:: python

   breeding_parameters['expected_f_two_mother_ids'] = expected_f_two_mother_ids
   breeding_parameters['expected_f_two_father_ids'] = expected_f_two_father_ids

   breeding_parameters['observed_f_two_mother_ids'] = observed_f_two_mother_ids
   breeding_parameters['observed_f_two_father_ids'] = observed_f_two_father_ids

   breeding_parameters['number_of_matches_f_two_mother_ids'] = sum(np.equal(expected_f_two_mother_ids, observed_f_two_mother_ids))
   breeding_parameters['number_of_matches_f_two_father_ids'] = sum(np.equal(expected_f_two_father_ids, observed_f_two_father_ids))

   assert breeding_parameters['number_of_matches_f_two_father_ids'] == 2000, "Incorrect father IDs."

   breeding_parameters['number_of_matches_f_two_mother_ids'] == 2000, "Incorrect mother IDs."

The function :func:`np.equal` checks if the IDs match by location so the order is preserved.
Otherwise the script will crash via an ``AssertionError``.


Parameter Set Stored Using :mod:`shelve`
========================================

I used the ``shelve`` module to store the parameters and entire mating history of ``standard_magic``.
I did the same thing for ``magic_1478``.

.. code-block:: python

   m1478_sim_parameters = shelve.open('magic_1478_simulation_parameters')
   m1478_sim_parameters['founders'] = founders
   m1478_sim_parameters['number_of_replicates'] = 5
   m1478_sim_parameters['prefounder_file_name'] = 'prefounders_1478.pop'
   m1478_sim_parameters['alleles'] = magic1478_alleles
   m1478_sim_parameters['operating_population_size'] = 2000
   m1478_sim_parameters['recombination_rates'] = [0.01]*1478
   m1478_sim_parameters.close()

   m1478_trait_parameters = shelve.open('magic_1478_trait_parameters')
   m1478_trait_parameters['number_of_qtl'] = 10
   m1478_trait_parameters['allele_effect_parameters'] = 1
   m1478_trait_parameters.close()

