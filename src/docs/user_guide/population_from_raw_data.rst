.. _population-from-raw-data:

###################################
Creating a Population From Raw Data
###################################

In general ``saegus`` was designed in mind for importing experimentally
determined genotype data. This tutorial shows how to convert raw data into a
``Population``.

.. _parsing_raw_data:

Parsing Raw Data
################

Our ``simuPOP.Population`` requires the population size, the genomic structure
(i.e. number and length of chromosomes) and genotypes for each individual.
We will use ``pandas`` and ``numpy`` to make to read and manipulate the raw
data. The genotype matrix and genetic map give us all the required information.
The recombination rates in the genetic map will be used once we simulate
mating.

.. code-block:: python
   :caption: Module imports

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='short', quiet=True)
   >>> import simuPOP as sim
   >>> import pandas as pd
   >>> import numpy as np
   >>> import collections
   >>> np.set_printoptions(suppress=True, precision=3)

.. _genotype_data:

Genotype Data
=============

The genotype data is in a file called ``genotype_matrix.txt``. There are 44445
loci and 105 individuals. We are making use of the very convenient
:mod:`pandas` :func:`read_csv` function. However, the ``pandas.DataFrame``
is immediately converted into a ``numpy.array``.

.. code-block:: python
   :caption: Read the file with the genotype matrix using ``pandas``

   >>> genotypes = np.array(pd.read_csv('example_genotype_matrix.txt', sep='\t', index_col=0))
   >>> print(genotypes)
   [['2/1' '3/2' '2/3' ..., '1/2' '1/3' '3/1']
    ['2/2' '3/3' '2/2' ..., '1/2' '1/1' '3/3']
    ['2/1' '3/3' '2/2' ..., '1/2' '1/1' '3/3']
    ...,
    ['2/2' '3/2' '2/2' ..., '2/2' '1/3' '1/1']
    ['2/2' '2/2' '2/2' ..., '1/2' '1/1' '3/3']
    ['2/2' '3/3' '2/2' ..., '1/2' '1/3' '3/3']]

Each row represents an individual. Each column represents a locus. For
example this is a truncated representation of individual ``1``'s genotype.

.. code-block:: python
   :caption: Example genotype

   >>> print(genotypes[0])
   ['2/1' '3/2' '2/3' ..., '1/2' '1/3' '3/1']

We need to transform the genotype data into a format which is acceptable to
:mod:`simuPOP`.

.. code-block:: python
   :caption: Converting ``numpy.array`` into ``Python.list``

   >>> converted_genotypes = [
   ...  [int(genotypes[ind, :][i][0]) for i in range(genotypes.shape[1])] +
   ...   [int(genotypes[ind, :][i][-1]) for i in range(genotypes.shape[1])] for ind in range(105)
   ... ]


.. _genetic_map:

Genetic Map
===========

The genetic map is parsed the same way as the genotype matrix.
:file:`example_genetic_map.txt` has columns:

+ locus
+ chromosome
+ cM

.. code-block:: python
   :caption: Parsing the genetic map

   >>> genetic_map = np.array(pd.read_csv('example_genetic_map.txt', sep='\t'))
   >>> print(genetic_map)
   [[     1.         1.        -5.511]
    [     2.         1.        -5.302]
    [     3.         1.        -5.3  ]
    ...,
    [ 44443.        10.        89.659]
    [ 44444.        10.        89.682]
    [ 44445.        10.        89.77 ]]

The :mod:`collections` allows us to easily obtain the genomic structure from
the genetic map. We will count how many loci are on each chromosome by using a
:class:`Counter` from :mod:`collections`.

.. code-block:: python
   :caption: Counting loci per chromosome

   >>> chromosome_column = np.array(genetic_map[:, 1], dtype=np.int)
   >>> print(chromosome_column)
   [ 1  1  1 ..., 10 10 10]
   >>> loci_counts = collections.Counter(chromosome_column)
   Counter({1: 6939, 2: 5171, 3: 4974, 5: 4838,
      4: 4819, 8: 3849, 7: 3775, 6: 3570, 9: 3337, 10: 3173})
   >>> chromosome_lengths = [loci_counts[i] for i in range(1, 11)]
   >>> print(chromosome_lengths)
   [6939, 5171, 4974, 4819, 4838, 3570, 3775, 3849, 3337, 3173]


.. warning::

   ``Counter`` may not be ordered the same way the data was entered

.. _creating_and_saving_the_population:

Creating and Saving the Population
##################################

Finally create an "empty" ``Population`` object and set the genotypes. We can
save the :class:`Population` object in native :mod:`simuPOP` format so we
do not have to re-do this step every single time we want to work with the
same population.

.. code-block:: python
   :caption: Creating a :class:`Population` from parsed data

   >>> example_pop = sim.Population(size=105, ploidy=2, loci=chromosome_lengths)
   >>> for i, ind in enumerate(example_pop.individuals()):
   ...      ind.setGenotype(converted_genotypes[i])

Let's examine ``example_pop`` to get a feel for :mod:`simuPOP`. :mod:`simuPOP`
has a distinct *feel* compared to most other Python packages.
:mod:`simuPOP` has a Python interface but it is really a C++ program. If you
are like the author of this walkthrough and Python is your first language
:mod:`simuPOP` can be intimidating. However, every single moment of frustration
pays off in both expected and unexpected ways. Make sure to thank the author
`Bo Peng`_ for all of his hard work in creating :mod:`simuPOP`.

.. _`Bo Peng` : //github.com/BoPeng/simuPOP

.. code-block:: python
   :caption: Examining a :class:`Population`

   >>> example_pop
   <simuPOP.Population>
   >>> print(example_pop.popSize())
   105
   >>> print(example_pop.numChrom())
   10
   >>> print(example_pop.numLoci())
   (6939, 5171, 4974, 4819, 4838, 3570, 3775, 3849, 3337, 3173)

It seems like the :class:`Population` has the correct structure. Let's examine
an individual.

.. code-block:: python
   :caption: Genotype data can be easily subsetted

   >>> example_individual = example_pop.individual(0)
   >>> example_genotype = np.array(example_individual.genotype(ploidy=0, chroms=0))
   >>> print(example_genotype)
   [2 3 2 ..., 3 3 2]

The examples to come will deepen our understanding of :mod:`simuPOP`. Finally
let's save our population in native :mod:`simuPOP` format.

.. code-block:: python
   :caption: Saving population for re-use

   >>> example_pop.save('example_pop.pop')
