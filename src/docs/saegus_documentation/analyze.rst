==============
:mod:`analyze`
==============



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

Usage
#####

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

**Example**

.. code-block:: python

   alleles
   array([[1, 2],
        [1, 3],
        [3, 1],
        ...,
        [1, 0],
        [3, 0],
        [3, 1]], dtype=int64)

   qtl
   [44, 103, 168, 340, 488, 639, 737, 819, 981, 1065]

   allele_effects
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

    aeframe

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


