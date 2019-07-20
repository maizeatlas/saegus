.. _structured_population_expansion:

####################################################
Expanding A Population According to Structure Matrix
####################################################

This example uses :py:mod:`simuPOP`'s capability to perform non-random
mating schemes We have inferred the population structure of ``example_pop``.
The goal of this example is to use the population structure to create a
simulated population with similar structure to the original population. The
genome of each individual in ``example_pop`` derives from a single
sub-population or a combination of sub-populations.

``population_structure_matrix.txt`` defines the proportion of the genome
inherited by an individual from each sub-population. For our example: each
individual is assigned a ``primary`` sub-population. The primary sub-population
is the sub-population from which the individual inherited the largest
proportion of their genome. The proportions of inheritance are interpreted as
probabilities for determining which sub-population the mate will derive from.

.. code-block:: python
   :caption: Module imports

   import simuOpt
   simuOpt.setOptions(alleleType='short', quiet=True)
   import simuPOP as sim
   import numpy as np, pandas as pd
   import collections as col
   from scipy import stats
   from saegus import parameters, breed, parse
   np.set_printoptions(suppress=True, precision=3)

We will continue to use the same population as the rest of our examples.

.. code-block:: python
   :caption: Load the population from ``example_pop.pop``

   example_pop = sim.loadPopulation('example_pop.pop')

We will add information fields to the population so that we can track the
pedigree. If we wanted to analyze the pedigree we could look at the mother
and father of each individual.

.. code-block:: python
   :caption: Adding information fields

   example_pop.addInfoFields(['ind_id', 'mother_id', 'father_id', 'primary'])
   sim.tagID(example_pop)
   example_pop.indInfo('ind_id')
   # (1.0, 2.0, 3.0, 4.0, 5.0, ... 105.0)

We will import a file that tells us the likely mating structure of each of the
105 individuals of our population.

.. code-block:: python
   :caption: Importing the population structure matrix

   structure_matrix = pd.read_csv('example_population_structure_matrix.txt', index_col=0)
   popst = parameters.PopulationStructure(example_pop)
   proportions = np.array(np.array(structure_matrix)[:, 1:7])
   print(proportions)
   # [[0.0, 0.9996, 0.0, 0.0004, 0.0, 0.0],
   # [0.011000000000000001, 0.0015, 0.0004, 0.1047, 0.0, 0.8824],
   # [0.0, 0.3832, 0.0, 0.6168, 0.0, 0.0],
   # ...
   # [0.0, 0.4602, 0.0, 0.5398, 0.0, 0.0]]

There is a very small amount of rounding error in the proportions for some
individuals. If the proportions do not sum to ``1`` then we cannot use
them to make a probability mass function. For example:

.. code-block:: python
   :caption: Example of rounding error

   proportions[33]
   # array([0.8856999999999998, 0.0016, 0.0009, 0.1065, 0.0042, 0.0011], dtype=object)
   sum(proportions[33])
   # 1.0000000000000002

So we will use a function to adjust the small difference from ``1`` by adding or
subtracting from the ``primary`` sub-population proportion.

.. code-block:: python
   :caption: Correcting the rounding error

   corrected_proportions = popst.correct_rounding_error(proportions)
   sum(corrected_proportions[33])
   # 0.9999999999999999

Apparently the result of ``0.9999999999999999`` is close enough for the
``scipy.stats`` module we are about to use. For peace of mind, we can use the
``name`` attribute of the ``stats.rv_discrete`` function to match the ``ind_id``
with the corresponding probabilities.

.. code-block:: python
   :caption: Creating the probability mass functions

   mating_pmfs = {}
   for i, ind in enumerate(example_pop.individuals()):
       mating_pmfs[ind.ind_id] = stats.rv_discrete(values=([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], 
           corrected_proportions[i]), name=str(ind.ind_id))
   
   example_pop.dvars().mating_probabilities = mating_pmfs

.. _validating_the_mating_probabilities:

Validating the Mating Probabilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before we proceed we should check the empirical distributions of the
probability mass functions. We will use an example individual who is quite
diverse in its lineage.

.. code-block:: python
   :caption: Comparing empirical results versus pmf

   corrected_proportions[5]
   # array([0.2195, 0.0198, 0.021, 0.2371, 0.1295, 0.3731], dtype=object)
   mating_pmfs[6].pk # corresponding mating pmf
   # array([0.2195, 0.0198, 0.021, 0.2371, 0.1295, 0.3731], dtype=object)
   mating_pmfs[6].name
   # 6.0

This individual is composed from all six sub-populations. We will draw
1000 times from the corresponding probability mass function and compare the
results.

.. code-block:: python
   :caption: Comparing empirical distribution

   draw_results = mating_pmfs[6].rvs(size=1000)
   draw_results
   # array([4, 3, 5, 3, 3, ... 4])
   draw_counts = col.Counter(draw_results)
   draw_frequencies = []
   for sp in range(6):
       draw_frequencies.append(draw_counts[sp]/1000)

Finally let's compare the ``1000`` draws with the probabilities.

.. code-block:: python
   :caption: Are they close?

   draw_frequencies
   # [0.219, 0.017, 0.021, 0.223, 0.148, 0.372]
   corrected_proportions[5]
   # array([0.2195, 0.0198, 0.021, 0.2371, 0.1295, 0.3731], dtype=object)

The draw frequencies are pretty close to the probability mass function. If we
increased the number of draws to 10,000 the differences would become even
smaller.

.. _assigning_primary_subpopulations:

Assigning Primary Subpopulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will continue by assigning each individual a primary sub-population. The
primary sub-population is the sub-population from which the majority of their
genome is derived.

.. code-block:: python
   :caption: Assignment of Primary Sub-Populations

   primary_subpops = {ind.ind_id: float(np.argmax(corrected_proportions[i]))
                      for i, ind in enumerate(example_pop.individuals())}

   for ind in example_pop.individuals():
       ind.primary = primary_subpops[ind.ind_id]

   example_pop.indInfo('primary')
   # (1.0, 5.0, 3.0, ..., 3.0)

Then we will use the virtual sub-population feature of ``simuPOP`` to group the
individuals without restricting mating between groups.

.. code-block:: python
   :caption: Split ``example_pop`` into virtual sub-populations

   primary_subpopulation_splitter = sim.InfoSplitter(field='primary', values=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
   example_pop.setVirtualSplitter(primary_subpopulation_splitter)

.. _parent_chooser_and_recombination_map:

Parent Chooser and Recombination Map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class containing the parent chooser function must be instantiated with the
expanded population size. The recombination map will be parsed with an older
function. We will explain in a later section more details about recombination
in :py:mod:`simuPOP`.

.. code-block:: python
   :caption: Instantiating parent chooser and parsing recombination map

   popst_parent_chooser = breed.ForcedPopulationStructureParentChooser(1000, mating_pmfs)
   tf = parse.TusonFounders()
   recom_rates = tf.parse_recombination_rates('genetic_map.txt')
   recom_rates
   # [0.0020926625899999962, 2.2615580000007186e-05, 0.00042822784999999361, ..., 0.0]

.. _expanding_the_population:

Expanding the Population
~~~~~~~~~~~~~~~~~~~~~~~~

Finally we can expand the population via mating according to the population
structure probability mass functions. Each mating event follows this process:

   1. Randomly draw the first parent
   2. Given the mating probability mass function of the first parent: draw the second parent from the probability mass function of the first parent
   3. Cross the two parents

This procedure is repeated 1,000 times because each mating event produces a
single offspring.

.. code-block:: python
   :caption: Expand the population to ``1000`` individuals

   example_pop.evolve(
       matingScheme=sim.HomoMating(
           sim.PyParentsChooser(popst_parent_chooser.forced_structure_parent_chooser),
           sim.OffspringGenerator(
               ops=[sim.IdTagger(), sim.PedigreeTagger(), sim.Recombinator(recom_rates)],
                   numOffspring=1),
           subPopSize=1000
       ),
       gen=1
   )

If we wanted to analyze the specific crosses we can create a pedigree using
the ``ind_id``, ``mother_id`` and ``father_id`` fields.

.. code-block:: python
   :caption: Create a pedigree

   pedigree = np.array((example_pop.indInfo('ind_id'),
                        example_pop.indInfo('mother_id'),
                        example_pop.indInfo('father_id'))).T
   print(pedigree)
   # [[  106.,    45.,    86.],
   # [  107.,    26.,    70.],
   # [  108.,    60.,    31.],
   # ...,
   # [ 1105.,    39.,    40.]]
