.. _collect_and_store_data:

######################
Collect and Store Data
######################

This is probably the most under-developed aspect of ``saegus``. I have written
up the data model (the way the data is stored) in a Google Doc. For the time
being we will design the functions to store and access the data around this
data model. I hope to move to a more permanent solution such as SciDB_ sooner
or later.

.. _SciDB: http://www.paradigm4.com/

.. _categories_of_data:

Categories of Data
##################

There is allele data, genotype data and possibly haplotype data. This example
creates an HDF5 file and adds each kind of data as it is generated. Then we
will show how to re-open the HDF5 file and retrieve the data in both ``Python``
and ``R``.

.. note:: The Google Doc about the data model goes into greater detail and gives examples.

This example is forcing me to fully develop the HDF5 storage functions.

.. _generating_the_data:

Generating the Data
===================

Since this is a simulation program we are free to make our own data up. Let's
import the necessary modules and load the example population like usual.

.. code-block:: python
   :caption: Module imports

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='short', quiet=True, numThreads=4)
   >>> import simuPOP as sim
   >>> import pandas as pd
   >>> import numpy as np
   >>> from saegus import analyze, parse
   >>> np.set_printoptions(precision=3, suppress=True)

Load the example population:

.. code-block:: python
   :caption: Load the example population, recombination map and initialize sex

   >>> example_pop = sim.loadPopulation('example_pop.pop')
   >>> sim.tagID(example_pop)
   >>> sim.initSex(example_pop)
   >>> tf = parse.TusonFounders()
   >>> recom_map = tf.parse_recombination_map('genetic_map.txt')

.. _allele_data:

Allele Data
===========

Now we have to determine the allele frequencies and declare the minor _vs_ major
alleles and alpha _vs_ omega alleles. A ``saegus`` function returns an array
with columns. Tied loci (loci with both alleles at 0.5 frequency) have the
alpha allele set at the minor allele and the omega allele set as the major
allele.

   locus | alpha | omega | minor | major
   =====================================

.. code-block:: python
   :caption: Allele data

   >>> allele_data_table = analyze.gather_allele_data(example_pop)
   >>> print(allele_data_table)
   [[     0.      1.      2.      1.      2.]
    [     1.      2.      3.      2.      3.]
    [     2.      2.      3.      3.      2.]
    ...,
    [ 44442.      1.      2.      2.      1.]
    [ 44443.      1.      3.      3.      1.]
    [ 44444.      1.      3.      1.      3.]]

The array for allele frequencies has the same ordering of columns:

   locus | alpha | omega | minor | major
   =====================================

.. code-block:: python
   :caption: Allele frequencies

   >>> allele_frequencies = analyze.gather_allele_frequencies(example_pop, allele_data_table)
   >>> print(allele_frequencies)
   [[     0.         0.319      0.681      0.319      0.681]
    [     1.         0.219      0.781      0.219      0.781]
    [     2.         0.938      0.062      0.062      0.938]
    ...,
    [ 44442.         0.533      0.467      0.467      0.533]
    [ 44443.         0.738      0.262      0.262      0.738]
    [ 44444.         0.267      0.733      0.267      0.733]]

.. _genotype_data:

Genotype Data
=============

Genotype data is stored in a different way than allele data. Genotype
frequencies are stored in a 3-dimensional array with axes:

   locus x alpha x omega

Where the frequency of genotype ``(1, 1)`` at locus ``0`` is ``(0, 1, 1)``. The
frequency data is stored in a ``numpy.ndarray``. We can collect the genotype
frequency array by using a ``saegus`` function.

.. code-block:: python
   :caption: Structure of genotype frequency data

   >>> genotype_frequencies = analyze.gather_genotype_frequencies(example_pop)
   >>> print(genotype_frequencies)
   [[[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.133  0.     0.     0.   ]
     [ 0.     0.371  0.495  0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.086  0.     0.   ]
     [ 0.     0.     0.267  0.648  0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.886  0.105  0.   ]
     [ 0.     0.     0.     0.01   0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    ...,
    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.305  0.457  0.     0.   ]
     [ 0.     0.     0.238  0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.562  0.     0.352  0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.086  0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.143  0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.248  0.     0.61   0.   ]
     [ 0.     0.     0.     0.     0.   ]]]

The syntax to access the frequency of genotype ``(1, 1)`` at locus ``0`` is

.. code-block:: python
   :caption: Accessing genotype frequencies

   >>> print(genotype_frequencies[0, 1, 1])
   0.133333333333

Unlike the allele data we do not have an organized array of genotypes by locus.
However, we can obtain all the genotypes as a set of coordinates by locus
using a very simple manipulation.

.. code-block:: python
   :caption: Genotypes as coordinates

   >>> genotypes_by_locus = np.array(np.ndarray.nonzero(genotype_frequencies)).T
   >>> print(genotypes_by_locus)
   [[    0     1     1]
    [    0     2     1]
    [    0     2     2]
    ...,
    [44444     1     1]
    [44444     3     1]
    [44444     3     3]]

.. note:: ``simuPOP`` considers ``(2, 1)`` and ``(1, 2)`` as distinct genotypes

This tells us that at locus ``0`` there are genotypes: ``(1, 1)``, ``(2, 1)``
and ``(2, 2)``. ``genotypes_by_locus`` is a 2-dimensional array. There are
a variable number of genotypes at each locus. At fixed sites there is only one
genotype. At segregating sites there may be up to ``4`` genotypes because
``simuPOP`` orders genotypes. Therefore, ``genotypes_by_locus`` has more
rows than the number of loci.

.. code-block:: python
   :caption: Variable number of genotypes by locus

   >>> print(genotypes_by_locus.shape)
   (122993, 3)

It is clear that the locus index will not match the ``genotypes_by_locus``
index. If we wanted to see the genotypes at a specific locus we can use the
``np.where`` function. For example if we wanted the genotypes present at locus
``5`` we would do:

.. code-block:: python
   :caption: Retrieve genotypes by locus

   >>> locus_five_genotypes = np.array(np.where(genotypes_by_locus[:, 0] == 5))
   >>> print(locus_five)
   [14, 15, 16]
   >>> print(genotypes_by_locus[locus_five_genotypes])
   [[5 1 1]
    [5 1 3]
    [5 3 3]]
   >>> print(genotypes_by_locus[locus_five_genotypes][:, 1:]) # without locus
   [[1 1]
    [1 3]
    [3 3]]

This tells us that at locus ``5`` there are genotypes ``(1, 1)``, ``(1, 3)``
and ``(3, 3)``. Let's check their frequencies.

.. code-block:: python
   :caption: Checking genotypic frequencies at locus ``5``

   >>> print(genotype_frequencies[5, 1, 1])
   0.904761904762
   >>> print(genotype_frequencies[5, 1, 3])
   0.0857142857143
   >>> print(genotype_frequencies[5, 3, 3])
   0.00952380952381

.. _storing_data_hdf5:

Storing Data in HDF5 Files
##########################

Our data take the form of arrays. Hierarchical Data Format 5 (``HDF5``) is a file
format optimized for 'lookup' operations. ``HDF5`` allow for
:math:`n`-dimensional arrays as well as metadata attached to HDF5 ``Groups``.
This part of this guide will demonstrate how to store allele data,
genotype data and the corresponding metadata.

.. _basics_of_hdf5:

Basics of Working with HDF5 and ``h5py``
========================================

HDF5 files can be navigated the same way as a directory. Every file has at
minimum a root directory: ``'/'``. ``numpy`` arrays can be directly stored
into HDF5 files as if you were working with a ``dict``.

.. code-block:: python
   :caption: Creating an HDF5 file

   >>> import h5py
   >>> example_data = h5py.File('example_data.hdf5')
   >>> allele_group = example_data.create_group('allele')
   >>> allele_group['states'] = allele_data # store data
   >>> allele_group['states']
   <HDF5 dataset "states": shape (44445, 5), type "<f8">
   >>> print(np.array(allele_group['states'])) # retrieve the data
   [[     0.      1.      2.      1.      2.]
    [     1.      2.      3.      2.      3.]
    [     2.      2.      3.      3.      2.]
    ...,
    [ 44442.      1.      2.      2.      1.]
    [ 44443.      1.      3.      3.      1.]
    [ 44444.      1.      3.      1.      3.]]


.. _groups_and_datasets:

Groups and Datasets
===================

HDF5 files are partitioned into any number of ``groups`` and ``datasets``.
Each ``group``` can have its own metadata stored as ``attributes``. A ``group``
can also contain any number of ``datasets``. But ``datasets`` cannot store
metadata as ``attributes``. We will store the column names associated with
the allele data array as an ``attribute``

.. code-block:: python
   :caption: HDF5 ``groups`` and attributes

   >>> allele_group = example_data.create_group('allele')
   >>> genotype_group = example_data.create_group('genotype')
   >>> allele_group.attrs['columns'] = list(map(np.string_, ['locus', 'alpha',
   ...                                           'omega', 'minor', 'major' ]))
   >>> print([name.decode('UTF-8') for name in allele_group.attrs['states']])
   ['locus', 'alpha', 'omega', 'minor', 'major']

We can store the allele frequency data and genotype frequency data in their
own groups.

.. code-block:: python
   :caption: Storing frequency data

   >>> allele_group['generation/founder'] = allele_frequencies
   >>> genotype_group = example_data.create_group('genotype')
   >>> genotype_group['generation/founder'] = genotype_frequencies # store
   >>> print(np.array(genotype_group['generation/founder'])) # retrieve
   [[[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.133  0.     0.     0.   ]
     [ 0.     0.371  0.495  0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.086  0.     0.   ]
     [ 0.     0.     0.267  0.648  0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.886  0.105  0.   ]
     [ 0.     0.     0.     0.01   0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    ...,
    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.305  0.457  0.     0.   ]
     [ 0.     0.     0.238  0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.562  0.     0.352  0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.     0.     0.086  0.   ]
     [ 0.     0.     0.     0.     0.   ]]

    [[ 0.     0.     0.     0.     0.   ]
     [ 0.     0.143  0.     0.     0.   ]
     [ 0.     0.     0.     0.     0.   ]
     [ 0.     0.248  0.     0.61   0.   ]
     [ 0.     0.     0.     0.     0.   ]]]

.. _creating_a_generation_of_data:

Data from Multiple Generations
==============================

We can easily store multiple generations of data