.. _multiple_generations_data_storage:

#########################################################
Integrated Example: Multiple Generations and Data Storage
#########################################################

This example will combine the what we have learned in prior examples. We will
create a multi-parental population, perform several generations of random
mating, calculate an additive trait at each generation, store the
allele frequencies, genotype frequencies and additive trait.


.. code-block:: python
   :caption: Module imports

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='short', numThreads=4, quiet=True)
   >>> import simuPOP as sim
   >>> import pandas as pd
   >>> import numpy as np
   >>> import random
   >>> from saegus import analyze, breed, operators, parameters, parse
   >>> np.set_printoptions(precision=3, suppress=True)

.. _creating_multi_parental_population:

Creating the Multi-Parental Population
######################################

.. _load_population_and_genetic_map:

Load Population and Genetic Map
===============================

As usual we start by re-loading the example population we have been working
with in all other tutorials.

.. code-block:: python
   :caption: Loading the population and genetic map

   >>> example_pop = sim.loadPopulation('example_pop.pop')
   >>> example_pop.addInfoFields(['ind_id', 'mother_id', 'father_id', 'g', 'p'])
   >>> sim.tagID(example_pop)
   >>> tf = parse.TusonFounders()
   >>> recom_map = tf.parse_recombination_rates('genetic_map.txt')

.. _achieving_recombinatorial_convergence:

Achieving Recombinatorial Convergence
=====================================


Creating the F\ :sub:`1`
------------------------

As in the :ref:`creating_multi_parental_populations` tutorial we will create
the multi-parental population from individuals ''1'' through ''8'' of
:py:class:`example_pop`. This requires three generations of mating to achieve
recombinatorial convergence. We will create a population of size ``1000``.

.. code-block:: python
   :caption: Repeating the MAGIC-like protocol

   >>> founders = [[1, 2], [3, 4], [5, 6], [7, 8]]
   >>> offspring_per_pair = 250
   >>> magic = breed.MAGIC(example_pop, recom_map)
   >>> magic.generate_f_one(founders, offspring_per_pair)

First Converging Cross
----------------------

.. code-block:: python
   :caption: Second step of recombinatorial convergence

   >>> first_random_cross = breed.RandomCross(example_pop, 4, 250)
   >>> first_mothers, first_fathers = first_random_cross.converging_random_cross()
   >>> first_parent_chooser = breed.SecondOrderPairIDChooser(first_mothers, first_fathers)
   >>> example_pop.evolve(
   ...  matingScheme=sim.HomoMating(
   ...      sim.PyParentsChooser(first_parent_chooser.snd_ord_id_pairs),
   ...      sim.OffspringGenerator(ops=[
   ...          sim.IdTagger(),
   ...          sim.PedigreeTagger(),
   ...          sim.Recombinator(rates=recom_map)],
   ...       numOffspring=1),
   ...      subPopSize=1000),
   ...   gen=1,
   ...  )
   1



Second Converging Cross
-----------------------

.. code-block:: python
   :caption: Second step of recombinatorial convergence

   >>> final_random_cross = breed.RandomCross(example_pop, 2, 500)
   >>> final_mothers, final_fathers = final_random_cross.converging_random_cross()
   >>> final_parent_chooser = breed.SecondOrderPairIDChooser(final_mothers, final_fathers)
   >>> example_pop.evolve(
   ...  matingScheme=sim.HomoMating(
   ...      sim.PyParentsChooser(final_parent_chooser.snd_ord_id_pairs),
   ...      sim.OffspringGenerator(ops=[
   ...          sim.IdTagger(),
   ...          sim.PedigreeTagger(),
   ...          sim.Recombinator(rates=recom_map)],
   ...       numOffspring=1),
   ...      subPopSize=1000),
   ...   gen=1,
   ...  )
   1

.. _parameterization_of_additive_trait:

Additive Trait
##############

We will choose ``20`` loci to declare as quantitative trait loci with
exponentially distributed allele effects with mean equal to ``1``.

.. math::

   G \sim Exp(1)

We will use the same process in :ref:`additive_trait_parameterization`.

.. code-block:: python
   :caption: Choosing QTL and assigning effects

   >>> segregating_loci = sim.stat(example_pop, numOfSegSites=sim.ALL_AVAIL, vars=['segSites'])
   >>> qtl = sorted(random.sample(segregating_loci, 20))
   >>> trait = parameters.Trait()
   >>> ae_table = trait.construct_allele_effects_table(example_pop, qtl, random.expovariate, 1)
   >>> ae_array = trait.construct_ae_array(ae_table, qtl)
   >>> print(ae_array[qtl])


Opening the HDF5 File and Declaring Groups
##########################################

All of the data derived from the simulation will be stored in a single HDF5
file. Each type of data will have a separate :py:class:`h5py.Group`. HDF5
groups make it very easy to split data into categories.

.. code-block:: python
   :caption: Set up the HDF5 File

   >>> integrated_example_data = h5py.File('integrated_example_data.hdf5')
   >>> allele_group = integrated_example_data.create_group('allele')
   >>> genotype_group = integrated_example_data.create_group('genotype')
   >>> trait_group = integrated_example_data.create_group('trait')


.. _ten_generations_of_random_mating:

Ten Generations of Random Mating
################################

This example will simulate ten generations of random mating with a population
size of ``1000``.

.. _collect_and_store_data_by_generation:

Operator Forms for Storing Data from Each Generation
====================================================

Just as :py:mod:`simuPOP` has function forms of its operators. :py:mod:`saegus`
has operator forms of its functions. There are operators that collect each
type of data and store it in an HDF5 file.

Allele Frequencies, Genotype Frequencies, ``g`` and ``p``
---------------------------------------------------------

Each kind of data is stored by generation as specified in the data model.
The operators in sim.evolve each take an :py:class:`h5py.Group` and acquire
the generation from the :py:class:`Population`. These two items are enough
to specify a unique address for the data inside the HDF5 file. The data
are stored in the :py:class:`h5py.File` generation by generation.

.. code-block:: python
   :caption: Creating the allele data and frequency arrays

   >>> allele_data_table =

.. code-block:: python
   :caption: Storing ten generations of data

   >>> example_pop.evolve(
   >>>  preOps=[
   ...      sim.Stat(alleleFreq=sim.ALL_AVAIL), #calculate allele frequencies
   ...      sim.Stat(genoFreq=sim.ALL_AVAIL), # calculate genotype frequencies
   ...      operators.HDF5AlleleFrequencies(allele_group), #store allele frequencies
   ...      operators.HDF5GenotypeFrequencies(genotype_group), # store genotype frequencies
   ...      operators.GenoAdditiveArray(qtl, ae_array), # calculate g
   ...      operators.PhenoAdditive(), # calulcate p
   ...      operators.HDFTrait('g', trait_group), # store g
   ...      operators.HDFTrait('p', trait_group), # store p
   ...           ],
   ...  matingScheme=sim.RandomMating(ops=[
   ...      sim.IdTagger(),
   ...      sim.PedigreeTagger(),
   ...      sim.Recombinator(rates=recom_rates)],
   ...      subPopSize=1000),
   ...  gen=10,
   ... )
