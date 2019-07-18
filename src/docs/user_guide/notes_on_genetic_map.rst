.. _notes_on_genetic_map:

######################################
Some Notes on Genetic Map in HDF5 File
######################################

A quick reference file for the properties of the genetic map file. The datasets
in ``genetic_map.hdf5`` are:

``genetic_map\``
   * chromosome_1
   * chromosome_2
   * chromosome_3
   * chromosome_4
   * chromosome_5
   * chromosome_6
   * chromosome_7
   * chromosome_8
   * chromosome_9
   * chromosome_10

Each ``chromosome`` dataset has 5 columns:

   absolute_locus, relative_locus, chromosome, cM, pairwise_distance

absolute_locus
==============

Absolute locus is a floating point number which refers to the index of the
locus of all 44, 445 loci. Note that the loci are indexed on the natural numbers
1, 2, 3 ... vs the Python index

relative_locus
==============

A floating point number which refers to the index of the locus on a particular
chromosome. Note that the loci are indexed on the natural numbers 1, 2, 3, ...
vs the Python index

chromsome
=========

A floating point number which refers to the chromosomes 1 through 10

cM
==

Position on genetic map given in centiMorgans

pairwise_distance
=================

A floating point number giving the absolute difference between a locus and
the locus immediately after. The difference is divided by 100 to yield a
probability of recombination for use in simuPOP.


Conversion from Text File to HDF5 File
######################################

The genetic map was given to me in text file format. I have recently converted
into an HDF5 file because I find myself referring to it frequently. Might
as well get totally on board with HDF5 and become comfortable with it.

Below is the code used to covert from text to HDF5

.. code-block:: python
   :caption: Creating the genetic map file

   >>> import numpy as np
   >>> import pandas as pd
   >>> import h5py
   >>> import simuOpt
   >>> simuOpt.setOptions(quiet=True, numThreads=4)
   >>> text_gen_map = pd.read_csv('example_genetic_map.txt', index_col=0, sep='\t')
   >>> text_gen_map.loc[text_gen_map['chr'] == 1]
   chr 	cM
   locus
   1 	1 	-5.511247
   2 	1 	-5.301981
   3 	1 	-5.299719
   4 	1 	-5.256896
   5 	1 	-5.225641
   6 	1 	-5.078748
   7 	1 	-5.057972
   8 	1 	-4.937512


.. code-block:: python
   :caption: Creating the HDF5 version

   >>> genetic_map = h5py.File("genetic_map.hdf5")
   >>> genetic_map.attrs['columns'] = list(map(np.string_, ['absolute_locus, relative_locus,'
   ...                                               'chromosome, cM, pairwise_distances']))
   >>> improved_genetic_map.attrs['columns']
   array([b'absolute_locus, relative_locus,chromosome, cM, pairwise_distances'],
      dtype='|S68')
   >>> for i in range(1, 11):
   ...      current_chromosome = np.array(text_gen_map.loc[text_gen_map['chr'] == str(i)])
   ...      improved_chromosome = np.zeros((current_chromosome.shape[0], 5))
   ...      improved_chromosome[:, 0] = current_chromosome[:, 0]
   ...      improved_chromosome[:, 1] = np.array([locus for locus in range(example_pop.numLoci(i-1))]) + 1
   ...      improved_chromosome[:, 2] = current_chromosome[:, 1]
   ...      improved_chromosome[:, 3] = current_chromosome[:, 2]
   ...      cM = current_chromosome[:, 2]
   ...      improved_chromosome[:, 4] = np.array([np.abs(cM[i-1] - cM[i]) for i in range(1, cM.shape[0])]
   ...                                + [0])/100
   ...      improved_genetic_map['chromosome_' + str(i)] = improved_chromosome
   >>> improved_genetic_map.close()
