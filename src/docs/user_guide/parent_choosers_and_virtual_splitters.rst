

###############
Parent Choosers
###############

Some of the most difficult parts of dealing with simuPOP derive from the
difficulty of properly implementing a non-random mating scheme. This is a log
of my attempts to use existing parent choosers and design new ones. I am still
trying to figure out the CombinedParentsChooser. I cannot quite seem
to grasp how to use separate virtual subpopulations. However, here are my
workarounds

Directional Selection
#####################

For a case of simple directional selection I can use the
:class:`PolyParentsChooser` by sorting the population according to the trait
under selection ``p`` . Then the :class:`ProportionSplitter` can be used to
determine the correct number of individuals to select.

.. code-block:: python
   :caption: Setting up data storage

   >>> data_file = h5py.File("822_selection_data_example.hdf5")
   >>> allele_group = data_file['allele']
   >>> data_file['allele'].attrs['columns'] = list(map(np.string_, ['locus', 'alpha', 'omega', 'minor', 'major']))
   >>> trait_group = data_file['trait']

.. code-block:: python
   :caption: Load the population and recombination map

   >>> example_pop = sim.loadPopulation('example_pop.pop')
   >>> example_pop.addInfoFields(['ind_id', 'mother_id', 'father_id', 'g', 'p', 'sample', 'selection'])
   >>> sim.tagID(example_pop)
   >>> sim.stat(example_pop, numOfSegSites=sim.ALL_AVAIL,
   ...      vars=['numOfSegSites','segSites', 'fixedSites'])
   >>> sim.stat(example_pop, alleleFreq=sim.ALL_AVAIL)
   >>> segregating_loci = example_pop.dvars().segSites
   >>> tf = parse.TusonFounders()
   >>> recom_map = tf.parse_recombination_rates("genetic_map.txt")

.. code-block:: python
   :caption: Allele data

   >>> allele_states = analyze.gather_allele_data(example_pop)
   >>> allele_frequencies = analyze.gather_allele_frequencies(example_pop,
   ...    allele_states)
   >>> gwas = analyze.GWAS(example_pop,
   ...    np.array(segregating_loci, dtype=np.int_),
   ...     allele_states[:, 3], 'example')

Calculating the Trait
=====================

To compute the trait values as usual.

.. code-block:: python
   :caption: Trait value

   >>> trait = parameters.Trait()
   >>> qtl = random.sample(segregating_loci, 30)
   >>> ae_table = trait.construct_allele_effects_table(example_pop, qtl, random.expovariate, 1)
   >>> ae_array = trait.construct_ae_array(ae_table, qtl)
   >>> operators.calculate_g(example_pop, ae_array)
   >>> operators.calculate_error_variance(example_pop, 0.7)
   >>> example_pop.dvars().epsilon
   >>> trait_group.attrs['error_variance'] = example_pop.dvars().epsilon # store the error variance
   >>> operators.calculate_p(example_pop)

Virtual Splitter on Proportion
==============================

.. code-block:: python
   :caption: Using a virtual splitter to facilitate selection

   >>> example_pop.sortIndividuals('p')
   >>> example_pop.setVirtualSplitter(
   ...     sim.ProportionSplitter([0.95, 0.05],
   ...         names=['not_selected', 'selected']))
   >>> example_pop.indInfo('p', [0, 1])
   (95.45393097104348,
    95.45455314831601,
    100.74438351659812,
    104.71413228108841,
    106.97581697850238)

Change the sex of all individuals in the 'selected' population to male. Change
the sex of all individual's in the 'non-selected' population to female.

.. code-block:: python

   >>> for ind in example_pop.individuals([0, 1]):
   ...    ind.setSex(1)
   >>> for ind in example_pop.individuals([0, 0]):
   ...    ind.setSex(2)


Using the PolyParentsChooser
############################

Now that the population has been divided according to male and female we can
easily use the :class:`PolyParentsChooser` to choose a random 'male' to mate
against :math:`n` random `females`.

.. note:: This is a far better implementation of the HalfSibBulkBalanceChooser

.. code-block:: python
   :caption: Using PolyParentsChooser

   >>> example_pop.evolve(
   ...   matingScheme=
   ...   sim.HomoMating(
   ...      sim.PolyParentsChooser(
   ...          polySex=sim.MALE, polyNum=10
   ...          ),
   ...      sim.OffspringGenerator(ops=[
   ...          sim.IdTagger(),
   ...          sim.PedigreeTagger(),
   ...          sim.Recombinator(rates=recom_map),
   ...      ]), subPopSize=100
   ...  ),
   ...  gen=1
   ... )

We can check the pedigree to make sure it worked properly.

.. code-block:: python
   :caption: Quick pedigree

   >>> pedigree = np.array([example_pop.indInfo('father_id'),
   ...        example_pop.indInfo('mother_id'),
   ...        example_pop.indInfo('ind_id')]).T
   >>> print(pedigree)
   [[  75.   80.  106.]
    [  75.    1.  110.]
    [  75.   73.  111.]
    [  75.   98.  112.]
    [  75.   92.  113.]
    [  75.   79.  114.]
    [  75.   11.  115.]
    [  75.   27.  116.]
    [  75.   63.  117.]
    [  75.   44.  118.]
    [  78.   77.  119.]
    [  78.   28.  120.]
    [  78.   45.  121.]
    [  78.   18.  122.]
    [  78.   18.  123.]
    [  78.   85.  124.]
    [  78.  100.  125.]
    [  78.   52.  126.]

So that's one way to do this kind of manipulation. However, the virtual
sub-population split is lost after each generation. So I can just re-apply it
after each period of breeding. The process is: compute ``p``,