.. _analysis_module:

==============
Analyze Module
==============




.. _gather_allele_data:

.. py:function:: gather_allele_data(pop)

   :param sim.Population pop: diploid ``simuPOP.population``
   :return: Array labeling alpha, omega, minor and major alleles

   .. warning:: Assumes alleles states are 1, 2, 3, 4

   Constructs a numpy.array with columns:
   locus   alpha   omega   minor   major

   Loci at 0.5 frequency have the minor allele set as the alpha allele and the
   major allele set as the omega allele. At present this function assumes that
   the population is diploid and each locus is bi-allelic. In the situation
   that a locus has a single allele an imaginary allele of ``0`` is assigned to
   the locus to make the rest of the data structures consistent.

   .. code-block:: python

      >>> allele_data = gather_allele_data(pop)
      >>> print(allele_data)
      [[     0.      1.      2.      1.      2.]
       [     1.      2.      3.      2.      3.]
       [     2.      2.      3.      3.      2.]
       ...,
       [ 44442.      1.      2.      2.      1.]
       [ 44443.      1.      3.      3.      1.]
       [ 44444.      1.      3.      1.      3.]]

.. _gather_allele_frequencies:

.. py:function:: gather_allele_frequencies(pop)

   :param sim.Population pop: diploid ``simuPOP.Population``
   :return: Array of allele frequencies of alpha, omega, minor and major

   .. warning:: Assumes allele states are 1, 2, 3, 4

   Constructs a ``np.array`` with columns:
   locus alpha_frequency   omega_frequency   minor_frequency   major_frequency

   Loci at 0.5 frequency have the minor allele set as the alpha allele and the
   major allele set as the omega allele. At present this function assumes that
   the population is diploid and each locus is bi-allelic. In the situation
   that a locus has a single allele an imaginary allele of ``0`` is assigned to
   the locus to make the rest of the data structures consistent.

.. code-block:: python

   >>> allele_data = gather_allele_data(pop)
   >>> allele_frequencies = gather_allele_frequencies(pop, allele_data)
   >>> print(allele_frequencies)
   [[     0.         0.319      0.681      0.319      0.681]
    [     1.         0.219      0.781      0.219      0.781]
    [     2.         0.938      0.062      0.062      0.938]
    ...,
    [ 44442.         0.533      0.467      0.467      0.533]
    [ 44443.         0.738      0.262      0.262      0.738]
    [ 44444.         0.267      0.733      0.267      0.733]]

.. _gather_genotype_frequencies:

.. py:function:: gather_genotype_frequencies(pop)

   :param sim.Population pop: diploid ``simuPOP.Population``
   :return: ``np.ndarray`` of genotype frequencies by locus

   Genotype data is stored in a different way than allele data. Genotype
   frequencies are stored in a 3-dimensional array with axes:

      locus x alpha x omega

   Where the frequency of genotype ``(1, 1)`` at locus ``0`` is ``(0, 1, 1)``. The
   frequency data is stored in a ``numpy.ndarray``. We can collect the genotype
   frequency array by using a ``saegus`` function.

.. code-block:: python

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



.. _single_generation:

SingleGeneration
================

.. py:class:: SingleGeneration






.. _gwas:

GWAS
====

.. py:class:: GWAS(pop, loci, run_id)


   .. py:method:: calculate_count_matrix(allele_subset, count_matrix_file_name=None)

      :parameter allele_subset:
      :parameter str count_matrix_file_name:

   .. py:method:: calc_kinship_matrix(allele_count_matrix, allele_frequencies, kinship_matrix_file_name)

      :parameter numpy.array allele_count_matrix: Minor/major allele copy number counts for each individual at each locus
      :parameter allele_frequencies: Minor/major allele frequencies for each locus. Used for Kinship (K) matrix count.
      :parameter kinship_matrix_file_name: Output file name to write TASSEL formatted K matrix with additional column for individual IDs

   .. py:method:: pop_struct_svd(count_matrix)

      :parameter count_matrix: numpy.array of minor allele counts of each individual

   .. py:method:: population_structure_formatter(eigen_data, pop_struct_file_name=None)

      :parameter dict eigen_data: Output from pop_struct_svd. Contains eigenvectors of PCA
      :parameter str pop_struct_file_name: File name to write first two components of PCA


   .. py:method:: hapmap_formatter(int_to_snp_conversions, hapmap_file_name)

      :parameter dict int_to_snp_conversions: Converts integer alleles to their corresponding string nucleotides
      :parameter str hapmap_file_name: Output file name to write tab-delimited columns

   .. py:method:: trait_formatter(trait_file_name=None)

      :parameter str trait_file_name: Output file name with tab-delimited columns and special TASSEL header.

   .. py:method:: replacement_trait_formatter(existing_trait_file_name, new_trait_file_name, new_trait_values)

      :parameter str existing_trait_file_name: Existing file of TASSEL formatted phenotype vector
      :parameter str new_trait_file_name: New file name written using replacement data
      :parameter new_trait_values: Values to replace existing phenotype values. Must be same number of values in existing_trait_file_name

.. _study:

Study
=====

.. py:class:: Study(run_id)

   .. py:method:: collect_samples(replicate_populations, sample_sizes)

      :parameter replicate_populations: simuPOP.Simulator with more than one population.
      :parameter sample_sizes: A list of integer valued sample sizes to take from each population. Multiple samples taken from each replicate.
      :return: Dictionary of lists of populations. Dictionary is keyed by ``population.dvars().rep``.

      .. code-block:: python
         :caption: Example of collect_samples

         >>> sample_sizes = [500, 600, 700, 800, 900, 1000,
         ...                    1100, 1200, 1300, 1400, 1500]
         >>> samples = Study.collect_samples(replicate_pops, sample_sizes)
         >>> samples
         {0: [<simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>],
         1: [<simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>,
         <simuPOP.Population>],

   .. py:method:: calculate_power_fpr(panel_map, sample_sizes, number_of_replicates, number_of_qtl)

      Determines the power by calculating number of detected loci divided by
      the number of loci with effects.

      :param panel_map: Dictionary of dictionaries of pandas.DataFrames. Keyed by panel_map[size][rep] = pd.DataFrame
      :param sample_sizes: List of integers corresponding to how many individuals are sampled from each replicate.
      :param number_of_replicates: Number of replicates in the run
      :param number_of_qtl: Loci declared as QTL and assigned an effect
      :return: pd.DataFrame summarizing power and false positive rate across replicates and sample sizes, lists of true positive loci detected in each run.


   .. py:method:: probability_of_detection(allele_effects_table, sample_sizes, number_of_replicates, true_positives_detected)

      Calculates the probability that a locus with an effect is detected.
      Probability of detection is defined as the number of times a locus is detected
      divided by the total number of realizations

      If the number of realizations is 200 and a locus is detected in all 200 realizations
      then its probability of detection is 1.0

      :param allele_effects_table: Allele effects table given by generate_allele_effects_table
      :param sample_sizes: List of number of individuals sampled from each replicate
      :param number_of_replicates: Number of replicates in the run
      :param true_positives_detected: Dictionary of lists of loci with effects that were detected.
      :return: Modified version of allele effects table which includes the probability of detection column.

      .. code-block:: python
         :caption: Example of the return value

         >>> prob_detection_table(aetable, sample_sizes, 20, true_positives_detected)
         <div>
         <table border="1" class="dataframe">
         <thead>
           <tr style="text-align: right;">
             <th></th>
             <th>locus</th>
             <th>alpha_allele</th>
             <th>alpha_effect</th>
             <th>beta_allele</th>
             <th>beta_effect</th>
             <th>difference</th>
             <th>detected</th>
           </tr>
         </thead>
         <tbody>
           <tr>
             <th>58</th>
             <td>96</td>
             <td>1</td>
             <td>3.079182</td>
             <td>3</td>
             <td>2.537866</td>
             <td>0.541317</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>274</th>
             <td>445</td>
             <td>0</td>
             <td>3.976630</td>
             <td>2</td>
             <td>5.201130</td>
             <td>1.224500</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>392</th>
             <td>619</td>
             <td>2</td>
             <td>2.087530</td>
             <td>3</td>
             <td>6.534154</td>
             <td>4.446624</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>431</th>
             <td>677</td>
             <td>2</td>
             <td>2.390493</td>
             <td>0</td>
             <td>4.353833</td>
             <td>1.963340</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>447</th>
             <td>703</td>
             <td>2</td>
             <td>4.543503</td>
             <td>0</td>
             <td>2.135412</td>
             <td>2.408091</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>620</th>
             <td>981</td>
             <td>0</td>
             <td>0.862903</td>
             <td>3</td>
             <td>4.536607</td>
             <td>3.673704</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>671</th>
             <td>1050</td>
             <td>3</td>
             <td>4.559900</td>
             <td>1</td>
             <td>0.713189</td>
             <td>3.846711</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>749</th>
             <td>1174</td>
             <td>2</td>
             <td>3.797462</td>
             <td>0</td>
             <td>1.208076</td>
             <td>2.589386</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>915</th>
             <td>1438</td>
             <td>2</td>
             <td>1.455625</td>
             <td>0</td>
             <td>2.069203</td>
             <td>0.613578</td>
             <td>0.0</td>
           </tr>
           <tr>
             <th>924</th>
             <td>1449</td>
             <td>0</td>
             <td>2.051093</td>
             <td>3</td>
             <td>0.869114</td>
             <td>1.181979</td>
             <td>0.0</td>
           </tr>
         </tbody>
         </table>
         </div>


.. _allele_data:

.. py:function:: allele_data(pop, alleles, loci)

   Determines the minor alleles, minor allele frequencies, major alleles and
   major allele frequencies.

   :parameter pop: Population intended for GWAS analysis
   :parameter list loci: Loci for which to calculate frequency
   :parameter dict alleles: Dictionary of alleles present at each locus

   This function is used to find the major/minor alleles of a Population
   ``pop`` given a list of ``alleles`` at each locus given in ``loci``.
   The output is intended to be used in other functions to determine the
   kinship matrix and population structure.

   Additionally this function will also resolve ties between the
   major and minor alleles which result when both alleles have exactly equal
   frequency i.e. 0.50.

.. code-block:: python

   pop = sim.loadPopulation('magic1478.pop')
   loci = list(range(pop.totNumLoci()))
   alleles = shelve.open('magic_1478_simulation_parameters')
   alleles

   {0: [1, 2],
    1: [1, 3],
    2: [3, 1],
    3: [0, 2],
    4: [2, 0],
    5: [0, 2],
    6: [0, 2],
    7: [3, 1],
    8: [0, 2],
    ...}

    af = analyze.allele_data(magic1478_2718, alleles, list(range(1478)))
    af

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: center;">
          <th></th>
          <th>minor_allele</th>
          <th>minor_frequency</th>
          <th>major_allele</th>
          <th>major_frequency</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>2</td>
          <td>0.00000</td>
          <td>1</td>
          <td>1.00000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>3</td>
          <td>0.13275</td>
          <td>1</td>
          <td>0.86725</td>
        </tr>
        <tr>
          <th>2</th>
          <td>1</td>
          <td>0.06575</td>
          <td>3</td>
          <td>0.93425</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2</td>
          <td>0.00000</td>
          <td>0</td>
          <td>1.00000</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0</td>
          <td>0.05675</td>
          <td>2</td>
          <td>0.94325</td>
        </tr>
        <tr>
          <th>5</th>
          <td>2</td>
          <td>0.24875</td>
          <td>0</td>
          <td>0.75125</td>
        </tr>
        <tr>
          <th>6</th>
          <td>2</td>
          <td>0.12300</td>
          <td>0</td>
          <td>0.87700</td>
        </tr>
        <tr>
          <th>7</th>
          <td>1</td>
          <td>0.00000</td>
          <td>3</td>
          <td>1.00000</td>
        </tr>
        <tr>
          <th>8</th>
          <td>2</td>
          <td>0.24000</td>
          <td>0</td>
          <td>0.76000</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
      </tbody>
    </table>
    <p>1478 rows × 4 columns</p>
    </div>

.. py:function:: rank_allele_effects(pop, loci, alleles, allele_effects)

   Collects information about alleles at quantitative trait loci into a
   dictionary. Determines favorable/unfavorable allele and corresponding
   frequency. Keys of quantitative_trait_alleles have similar hierarchy
   for both the alleles and their frequencies.

   :param pop:
   :param loci:
   :param alleles:
   :param allele_effects:

.. py:function:: allele_frq_table(pop, number_gens, allele_frq_data, recombination_rates, genetic_map)

   Tabulates useful information about each locus and allele frequency

   :param pop: Population with multiple sub-populations. Usually represents multiple generations of recurrent selection or drift.
   :param int number_gens: Number of generations of selection or drift
   :param dict allele_frq_data: Allele frequency data and the major/minor alleles at each locus.
   :param list recombination_rates: Recombination rates for each locus in order.
   :param genetic_map: Chromosome:cM position correspondence.


.. code-block:: python

   allele_data = analyze.Frq(pop, triplet_qtloci[0], alleles_by_locus, qt_allele_effects[0])
   allele_frequencies = allele_data.allele_frequencies(pop, range(pop.totNumLoci())
   allele_frequency_table = selection_qtd.allele_frq_table(pop, 10, allele_frq_data, recombination_rates,
                                                         genetic_map)


.. py:function:: generate_allele_effects_table(qtl, founder_alleles, allele_effects)

 Creates a simple pd.DataFrame for allele effects. Hard-coded
 for bi-allelic case.

    :parameter list qtl: List of loci declared as QTL
    :parameter np.array alleles: Array of alleles at each locus
    :parameter dict allele_effects: Mapping of effects for alleles at each QTLocus

.. code-block:: python
   :caption: Example of an allele effects table

   >>> alleles
   array([[1, 2],
        [1, 3],
        [3, 1],
        ...,
        [1, 0],
        [3, 0],
        [3, 1]], dtype=int64)

   >>> qtl
   [44, 103, 168, 340, 488, 639, 737, 819, 981, 1065]

   >>> allele_effects
   {44: {0: 5.629446187924926, 2: 1.8962727055819322},
   103: {0: 1.3097813991257303, 2: 6.14070564290979},
   168: {2: 6.718096248082958, 3: 4.697238579652859},
   340: {1: 1.521689147484636, 2: 2.2131077852927032},
   488: {1: 2.512286137462885, 3: 2.486777318327935},
   639: {0: 1.1268072986309254, 3: 1.3391282487711016},
   737: {0: 1.4879865577936147, 1: 1.607534785598338},
   819: {1: 2.2153417608326986, 3: 0.20077940947200731},
   981: {0: 3.9513501430851568, 3: 1.78843909724396},
   1065: {0: 0.998194377898828, 2: 1.5139052352904945}}

    >>> aeframe

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>locus</th>
          <th>alpha_allele</th>
          <th>alpha_effect</th>
          <th>beta_allele</th>
          <th>beta_effect</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>44</td>
          <td>0</td>
          <td>5.629446</td>
          <td>2</td>
          <td>1.896273</td>
        </tr>
        <tr>
          <th>1</th>
          <td>103</td>
          <td>0</td>
          <td>1.309781</td>
          <td>2</td>
          <td>6.140706</td>
        </tr>
        <tr>
          <th>2</th>
          <td>168</td>
          <td>2</td>
          <td>6.718096</td>
          <td>3</td>
          <td>4.697239</td>
        </tr>
        <tr>
          <th>3</th>
          <td>340</td>
          <td>2</td>
          <td>2.213108</td>
          <td>1</td>
          <td>1.521689</td>
        </tr>
        <tr>
          <th>4</th>
          <td>488</td>
          <td>3</td>
          <td>2.486777</td>
          <td>1</td>
          <td>2.512286</td>
        </tr>
        <tr>
          <th>5</th>
          <td>639</td>
          <td>0</td>
          <td>1.126807</td>
          <td>3</td>
          <td>1.339128</td>
        </tr>
        <tr>
          <th>6</th>
          <td>737</td>
          <td>1</td>
          <td>1.607535</td>
          <td>0</td>
          <td>1.487987</td>
        </tr>
        <tr>
          <th>7</th>
          <td>819</td>
          <td>1</td>
          <td>2.215342</td>
          <td>3</td>
          <td>0.200779</td>
        </tr>
        <tr>
          <th>8</th>
          <td>981</td>
          <td>0</td>
          <td>3.951350</td>
          <td>3</td>
          <td>1.788439</td>
        </tr>
        <tr>
          <th>9</th>
          <td>1065</td>
          <td>2</td>
          <td>1.513905</td>
          <td>0</td>
          <td>0.998194</td>
        </tr>
      </tbody>
    </table>
    </div>


.. _multi_generation:

MultiGeneration
===============

.. py:class:: MultiGeneration(run_id)


   .. _multi_generation_collect_allele_frequency_data:

   .. py:method:: collect_allele_frequency_data(meta_population_library, minor_alleles)

      :parameter dict meta_population_library: Dictionary of lists of simuPOP.Populations
      :parameter minor_alleles: A tuple, list or array of the minor alleles at each locus

      Generates an array of the minor allele frequencies of each replicate at each
      generation. This is the *old* way of doing things. But it is still useful because
      it is designed to be written to a text file.

      Columns are: replicate, generation, locus1, locus2, ..., locusN

      .. code-block:: py
         :caption: Collecting allele frequency data for a writable text file

         >>> mafs = collect_allele_frequency_data(meta_populations, minor_alleles)
         >>> print(mafs)
         [[  0.   ,   0.   ,   0.325, ...,   0.435,   0.27 ,   0.255],
          ...,
          [  4.   ,  10.   ,   0.165, ...,   0.465,   0.035,   0.035]]

   .. _multi_generation_store_allele_frequency_data:

   .. py:method:: store_allele_frequency_data(meta_population_library, hdf_file_name)

      :parameter meta_population_library: Dict of lists of simuPOP.Populations
      :parameter str hdf_file_name: File name to write output

       Collects minor allele frequency data of a multiple generation
       population library. Stores the allele frequency data in an
       HDF5 file.

       af/replicate_id/generation_id

      .. code-block:: py
         :caption: Storing and accessing alelle frequency data in an HDF5 file

         >>> minor_af_data = h5py.File("example_af_data.hdf5")
         >>> minor_af_data
         <HDF5 file "example_af_data.hdf5" (mode r+)>
         >>> list(minor_af_data.keys())
         ['af']
         >>> minor_af_data['af']['0'] # replicate 0
         <HDF5 group "/af/0" (6 members)>

      If we wanted to make an array out of all the generations within a replicate
      we can use a generator expression, list comprehension or a loop to make a
      list of lists. For example if we wanted to put the generational data into
      a :py:class:`np.array`.

      .. warning::

         HDF5 files do not store data in the same order it was inserted.
         If we want to have the generations in order we need to do an
         extra step.

      .. code-block:: py
         :caption: Extract allele frequencies into a numpy array

         >>> generations = tuple(map(str, range(0, 11, 2)))
         >>> generations
         ('0', '2', '4', '6', '8', '10')
         >>> minor_allele_frequencies = np.asarray((tuple(np.asarray(minor_af_data['af']['0']) for gen in generations)))
         >>> minor_allele_frequencies # the rows are generations columns are loci
         array([[ 0.325,  0.18 ,  0.05 , ...,  0.435,  0.27 ,  0.255],
          [ 0.275,  0.255,  0.07 , ...,  0.36 ,  0.095,  0.08 ],
          [ 0.315,  0.175,  0.105, ...,  0.34 ,  0.125,  0.09 ],
          [ 0.32 ,  0.13 ,  0.115, ...,  0.275,  0.02 ,  0.015],
          [ 0.34 ,  0.185,  0.215, ...,  0.35 ,  0.025,  0.   ],
          [ 0.375,  0.075,  0.26 , ...,  0.315,  0.   ,  0.   ]])

   .. _collect_heterozygote_frequency_data:

   .. py:method:: collect_heterozygote_frequency_data(meta_population_library)

      :parameter meta_population_library: Dictionary of lists of simuPOP.Populations

      Collects heterozygote frequency data from the
      populations in ``meta_population_library``. The data is collected
      into a :class:`np.array` which is suitable for writing to a text file. The
      columns of the array are:

      + replicate
      + generation
      + locus1
      + locus2
      + so on and so forth

      .. code-block:: py
         :caption: Collecting heterozygote data from samples

         >>> hetf = collect_heterozygote_frequency_data(meta_population_library)
         >>> print(hetf)
         [[  0.     0.     0.45 ...,   0.39   0.26   0.31]
         [  0.     2.     0.35 ...,   0.46   0.19   0.16]
         [  0.     4.     0.51 ...,   0.44   0.21   0.14]
         ...,
         [  4.     6.     0.26 ...,   0.5    0.09   0.09]
         [  4.     8.     0.39 ...,   0.46   0.14   0.14]
         [  4.    10.     0.31 ...,   0.51   0.07   0.07]]

   .. _store_heterozygote_frequency_data:

   .. py:method:: store_heterozygote_frequency_data(meta_population_library, hdf_file_name)

      :parameter meta_population_library: Dict of lists of simuPOP.Populations
      :parameter str hdf_file_name: Output file name

      Stores heterozygote frequency data in and HDF5 file. The data are stored
      keyed as

         hetf/replicate/generation


      :parameter meta_population_library: Dict of lists of simuPOP.Populations
      :parameter str hdf_file_name: File name to write output

       Collects minor allele frequency data of a multiple generation
       population library. Stores the allele frequency data in an
       HDF5 file.

       hetf/replicate_id/generation_id

      .. code-block:: py
         :caption: Storing and accessing heterozygote frequency data in an HDF5 file

         >>> store_heterozygote_frequency_data(meta_population_library, "example_hetf_data.hdf5")
         >>> hetf_data = h5py.File("example_hetf_data.hdf5")
         >>> hetf_data
         <HDF5 file "example_hetf_data.hdf5" (mode r+)>
         >>> list(hetf_data.keys())
         ['hetf']
         >>> hetf_data['hetf']['0'] # replicate 0
         <HDF5 group "/hetf/0" (6 members)>

      If we wanted to make an array out of all the generations within a replicate
      we can use a generator expression, list comprehension or a loop to make a
      list of lists. For example if we wanted to put the generational data into
      a :py:class:`np.array`.

      .. warning::

         HDF5 files do not store data in the same order it was inserted.
         If we want to have the generations in order we need to do an
         extra step.

      .. code-block:: py
         :caption: Extract heterozygote frequencies into a numpy array

         >>> hetf_data = h5py.File("example_hetf_data.hdf5")
         >>> generations = tuple(map(str, range(0, 11, 2)))
         >>> generations
         ('0', '2', '4', '6', '8', '10')
         >>> het_frequencies = np.asarray((tuple(np.asarray(hetf_data['af']['0']) for gen in generations)))
         >>> het_frequencies # the rows are generations columns are loci
         array([[ 0.325,  0.18 ,  0.05 , ...,  0.435,  0.27 ,  0.255],
          [ 0.275,  0.255,  0.07 , ...,  0.36 ,  0.095,  0.08 ],
          [ 0.315,  0.175,  0.105, ...,  0.34 ,  0.125,  0.09 ],
          [ 0.32 ,  0.13 ,  0.115, ...,  0.275,  0.02 ,  0.015],
          [ 0.34 ,  0.185,  0.215, ...,  0.35 ,  0.025,  0.   ],
          [ 0.375,  0.075,  0.26 , ...,  0.315,  0.   ,  0.   ]])
         >>> hetf_data.close()

.. _definition_collect_genotype_phenotype_data:

   .. py:method:: collect_genotype_phenotype_data(meta_population_library)

      :parameter meta_population_library: Dict of lists of simuPOP.Populations

      Collects the genotype and phenotype data of a multiple replicate
      multiple sample population dictionary. The resulting data is
      a single array. Each row has ind_id, replicate, generation, g and p.

      .. note::

         Assumes that the population has infoFields ``g`` and ``p`` defined.

      .. code-block:: py
         :caption: Example of input and output

         >>> meta_population_library
         {0: [<simuPOP.Population>, ..., <simuPOP.Population>],
         ...,
         1: [<simuPOP.Population>, ..., <simuPOP.Population>]}
         >>> geno_pheno_data = collect_genotype_phenotype_data(meta_population_library)
         >>> print(geno_pheno_data)
         [[   117.         0.         0.        90.311     62.455]
          [   122.         0.         0.        90.889    101.073]
          [   126.         0.         0.        90.194     77.146]
          ...,
          [ 80084.         4.        10.       124.4      148.832]
          [ 80096.         4.        10.       129.004    100.359]
          [ 80100.         4.        10.       123.914    133.201]]

   .. _definition_store_genotype_phenotype_data:

   .. py:method:: store_genotype_phenotype_data(meta_population_library, hdf5_file_name)

      :parameter meta_population_library: Dict of lists of simuPOP.Populations
      :parameter str hdf5_file_name: Output file name

      Collects the genotype and phenotype data of a multiple replicate
      multiple sample population dictionary. Stores the results in
      an HDF5 file.

      Keyed as

         geno_pheno/replicate_id/generation_id

      .. code-block:: py
         :caption: Storing and accessing geno pheno data in an HDF5 file

         >>> store_genotype_phenotype_data(meta_population_library, "example_geno_pheno_data.hdf5")
         >>> gp_data = h5py.File("example_geno_pheno_data.hdf5")
         >>> gp_data
         <HDF5 file "example_geno_pheno_data.hdf5" (mode r+)>
         >>> list(gp_data.keys())
         ['geno_pheno']
         >>> gp_data['hetf']['0'] # replicate 0
         <HDF5 group "/geno_pheno/0" (6 members)>

      If we wanted to make an array out of all the generations within a replicate
      we can use a generator expression, list comprehension or a loop to make a
      list of lists. For example if we wanted to put the generational data into
      a :py:class:`np.array`. The resulting array has dimensions

         generations x sample_size x data_columns

      .. warning::

         HDF5 files do not store data in the same order it was inserted.
         If we want to have the generations in order we need to do an
         extra step.

      .. code-block:: py
         :caption: Extract genotype/phenotype data into a numpy array

         >>> gp_data = h5py.File("example_geno_pheno_data.hdf5")
         >>> generations = tuple(map(str, range(0, 11, 2)))
         >>> generations
         ('0', '2', '4', '6', '8', '10')
         >>> gp_zero = np.asarray((tuple(np.asarray(gp_data['geno_pheno']['0'])
         ...                          for gen in generations)))
         >>> print(gp_zero)
         [[[   117.         0.         0.        90.311     62.455]
           ...,
           [  1102.         0.         0.        83.207     98.937]]

          [[ 12631.         0.         2.       116.315    102.098]
           ...,
           [ 14084.         0.         2.        96.314     96.24 ]]

          [[ 27620.         0.         4.       117.47     133.751]
           ...,
           [ 29098.         0.         4.       114.059    109.896]]

          [[ 42609.         0.         6.       122.617    117.903]
           ...,
           [ 44077.         0.         6.       120.406    120.769]]

          [[ 57615.         0.         8.       123.669    163.46 ]
           ...,
           [ 59084.         0.         8.       124.701    123.834]]

          [[ 72622.         0.        10.       122.074    135.145]
           ...,
           [ 74059.         0.        10.       122.845    118.8  ]]]
         >>> gp_data.close()

      We can use the ``with`` key word so we don't have to worry about closing the
      file after we are done with it.

      .. code-block:: py
         :caption: Accessing data using the context manger: ``with``

         >>> with h5py.File('example_geno_pheno_data.hdf5') as exgp_file:
         ...   gp_zero = np.asarray(tuple(exgp_file['geno_pheno']['0'][gen] for gen in generations))

.. _definition_store_genotype_frequency_data:

   .. py:method:: store_genotype_frequency_data(meta_population_library, minor_alleles, hdf_file_name)

      :parameter meta_population_library: Dict of lists of simuPOP.Populations
      :parameter minor_alleles: A list of the minor alleles at each locus.
      :parameter str hdf_file_name: Output file name

      Collects the frequency of the minor allele homozygote data
      of a multiple replicate multiple sample population dictionary. The minor
      allele genotypes are created using the ``minor_alleles`` parameter.
      Stores the results in an HDF5 file.

      Keyed by

        homf/replicate_id/generation_id

   .. code-block:: py
      :caption: Example of storing genotype frequency data


.. _definition_generate_allele_effects_table:

.. py:function:: generate_allele_effects_table(population_allele_frequencies, allele_array, allele_effect_array):

   :parameter dict population_allele_frequencies: Allele frequencies keyed by locus
   :parameter np.array allele_array: Array where rows are loci and columns are alleles
   :parameter np.array allele_effect_array: Array where rows are loci and the columns are effects.

   Creates a pandas DataFrame with the columns:
   + alpha allele
   + alpha allele effect
   + alpha allele frequency
   + beta allele
   + beta allele effect
   + beta allele frequency

   .. warning::

      Assumes di-allelic case


   .. code-block:: py
      :caption: Examples of input parameters

      >>> population_allele_frequencies
      {0: defdict({1: 0.9807692307692307, 2: 0.019230769230769232}),
      1: defdict({1: 0.8461538461538461, 3: 0.15384615384615385}),
      2: defdict({1: 0.07692307692307693, 3: 0.9230769230769231}),
      3: defdict({0: 0.9230769230769231, 2: 0.07692307692307693}),
      4: defdict({0: 0.019230769230769232, 2: 0.9807692307692307}),
      5: defdict({0: 0.9230769230769231, 2: 0.07692307692307693}),
      6: defdict({0: 0.75, 2: 0.25}),
      ...,
      }
      >>> print(allele_array)
      [[1 2]
       [1 3]
       [3 1]
       ...,
       [1 0]
       [3 0]
       [3 1]]
      >>> qtl = sorted(tuple(random(sample(range(1478), 10)))
      >>> print(allele_effect_array[qtl])
      [[ 1.892  0.179  0.     0.     0.     0.   ]
       [ 0.92   1.     0.     0.     0.     0.   ]
       [ 0.079  0.     0.     1.653  0.     0.   ]
       [ 0.118  1.263  0.     0.     0.     0.   ]
       [ 3.731  0.     2.626  0.     0.     0.   ]
       [ 0.     0.673  0.     0.417  0.     0.   ]
       [ 0.418  0.     0.     1.94   0.     0.   ]
       [ 0.     0.6    0.     0.175  0.     0.   ]
      ...,
      ]

   .. code-block:: py
      :caption: Example usage

      >>> generate_allele_effects_table(population_allele_frequencies, allele_array, allele_effect_array)

   .. raw:: html

      <table border="1" class="dataframe">
        <thead>
          <tr style="text-align: right;">
            <th></th>
            <th>alpha</th>
            <th>alpha_effect</th>
            <th>alpha_frequency</th>
            <th>beta</th>
            <th>beta_effect</th>
            <th>beta_frequency</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th>0</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.980769</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.019231</td>
          </tr>
          <tr>
            <th>1</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.846154</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.153846</td>
          </tr>
          <tr>
            <th>2</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.923077</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>3</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.923077</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>4</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.980769</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.019231</td>
          </tr>
          <tr>
            <th>5</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.923077</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>6</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.750000</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.250000</td>
          </tr>
          <tr>
            <th>7</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>8</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.846154</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.153846</td>
          </tr>
          <tr>
            <th>9</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>10</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.730769</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.269231</td>
          </tr>
          <tr>
            <th>11</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.923077</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>12</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.788462</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.211538</td>
          </tr>
          <tr>
            <th>13</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>14</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>15</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.538462</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.461538</td>
          </tr>
          <tr>
            <th>16</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.923077</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>17</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>18</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>19</th>
            <td>5</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>4</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>20</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.538462</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.461538</td>
          </tr>
          <tr>
            <th>21</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.769231</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.230769</td>
          </tr>
          <tr>
            <th>22</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.980769</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.019231</td>
          </tr>
          <tr>
            <th>23</th>
            <td>5</td>
            <td>0.000000</td>
            <td>0.519231</td>
            <td>4</td>
            <td>0.000000</td>
            <td>0.480769</td>
          </tr>
          <tr>
            <th>24</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>25</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.692308</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.307692</td>
          </tr>
          <tr>
            <th>26</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.980769</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.019231</td>
          </tr>
          <tr>
            <th>27</th>
            <td>0</td>
            <td>1.891549</td>
            <td>0.980769</td>
            <td>1</td>
            <td>0.179440</td>
            <td>0.019231</td>
          </tr>
          <tr>
            <th>28</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.884615</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.115385</td>
          </tr>
          <tr>
            <th>29</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.653846</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.346154</td>
          </tr>
          <tr>
            <th>...</th>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
          </tr>
          <tr>
            <th>1448</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.923077</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>1449</th>
            <td>0</td>
            <td>0.415928</td>
            <td>0.653846</td>
            <td>3</td>
            <td>0.921988</td>
            <td>0.346154</td>
          </tr>
          <tr>
            <th>1450</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.730769</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.269231</td>
          </tr>
          <tr>
            <th>1451</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.730769</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.269231</td>
          </tr>
          <tr>
            <th>1452</th>
            <td>4</td>
            <td>0.000000</td>
            <td>0.980769</td>
            <td>5</td>
            <td>0.000000</td>
            <td>0.019231</td>
          </tr>
          <tr>
            <th>1453</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>1454</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.923077</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>1455</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.538462</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.461538</td>
          </tr>
          <tr>
            <th>1456</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.942308</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.057692</td>
          </tr>
          <tr>
            <th>1457</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.942308</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.057692</td>
          </tr>
          <tr>
            <th>1458</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.769231</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.230769</td>
          </tr>
          <tr>
            <th>1459</th>
            <td>4</td>
            <td>0.000000</td>
            <td>0.942308</td>
            <td>5</td>
            <td>0.000000</td>
            <td>0.057692</td>
          </tr>
          <tr>
            <th>1460</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.807692</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.192308</td>
          </tr>
          <tr>
            <th>1461</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.980769</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.019231</td>
          </tr>
          <tr>
            <th>1462</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>1463</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.538462</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.461538</td>
          </tr>
          <tr>
            <th>1464</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>1465</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>1466</th>
            <td>1</td>
            <td>1.176202</td>
            <td>0.923077</td>
            <td>3</td>
            <td>0.260720</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>1467</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.884615</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.115385</td>
          </tr>
          <tr>
            <th>1468</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.634615</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.365385</td>
          </tr>
          <tr>
            <th>1469</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>1470</th>
            <td>2</td>
            <td>0.000000</td>
            <td>0.692308</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.307692</td>
          </tr>
          <tr>
            <th>1471</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.923077</td>
            <td>3</td>
            <td>0.000000</td>
            <td>0.076923</td>
          </tr>
          <tr>
            <th>1472</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>1473</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.865385</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.134615</td>
          </tr>
          <tr>
            <th>1474</th>
            <td>0</td>
            <td>0.000000</td>
            <td>0.807692</td>
            <td>2</td>
            <td>0.000000</td>
            <td>0.192308</td>
          </tr>
          <tr>
            <th>1475</th>
            <td>1</td>
            <td>0.000000</td>
            <td>0.884615</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.115385</td>
          </tr>
          <tr>
            <th>1476</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>0</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
          <tr>
            <th>1477</th>
            <td>3</td>
            <td>0.000000</td>
            <td>0.961538</td>
            <td>1</td>
            <td>0.000000</td>
            <td>0.038462</td>
          </tr>
        </tbody>
      </table>


.. _definition_minor_allele_frequencies_table:

.. py:method:: minor_allele_frequencies_table(population_allele_frequencies, minor_alleles)

   :parameter dict population_allele_frequencies: Allele frequencies by locus
   :parameter minor_alleles: Array or list of minor alleles

   Returns a pandas DataFrame of the minor alleles and their frequencies.
   Expects a set of allele frequencies from simuPOP's Stat class.

   .. code-block:: py
      :caption: Example usage

      >>> sim.stat(pop, alleleFreq=sim.ALL_AVAIL)
      >>> population_allele_frequencies = pop.dvars().alleleFreq
      >>> population_allele_frequencies
      {0: defdict({1: 0.9807692307692307, 2: 0.019230769230769232}),
       1: defdict({1: 0.8461538461538461, 3: 0.15384615384615385}),
       2: defdict({1: 0.07692307692307693, 3: 0.9230769230769231}),
       3: defdict({0: 0.9230769230769231, 2: 0.07692307692307693}),
      ...,
      }
      >>> mafrqs = minor_allele_frequencies_table(population_allele_frequencies, minor_alleles)
      >>> print(mafrq)
            minor_allele  minor_frequency
      0                2         0.019231
      1                3         0.153846
      2                1         0.076923
      3                2         0.076923
      4                0         0.019231
      5                2         0.076923
      ...

      [1478 columns x 2 rows]

