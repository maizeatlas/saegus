===============================================
Example of a Complete Run: Simulator to Results
===============================================


.. code:: python

    import simuOpt
    simuOpt.setOptions(alleleType='short', optimized=True, numThreads=4, quiet=True)
    import simuPOP as sim
    import pandas as pd
    from saegus import breed, operators, simulate, analyze, parse, parameters
    import shelve
    import numpy as np
    import random
    import collections as col
    np.set_printoptions(suppress=True, precision=3)


.. _initial_parameters:

Initial Parameters
==================

.. code-block:: python
   :caption: Setting parameters

   >>> run_id='demonstration'
   >>> number_of_replicates = 5
   >>> operating_population_size = 2000
   >>> sample_sizes=[250, 500, 750, 1000]
   >>> number_of_qtl = 10
   >>> founders = [[2, 26], [3, 25], [4, 24], [5, 23]]
   >>> offspring_per_pair = 500
   >>> recombination_rates = [0.01]*1478

.. _magic_internals:

An Inside Look at the MAGIC Mating Scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python
   :caption: Internals of MAGIC breeding scheme

   >>> prefounders = sim.loadPopulation('prefounders1478.pop')
   >>> multi_prefounders = sim.Simulator(prefounders, number_of_replicates, stealPops=False)
   >>> magic = breed.MAGIC(multi_prefounders, founders, recombination_rates)
   >>> sim.tagID(prefounders, reset=27) # first 26 ids are reserved for the prefounders
   >>> magic.generate_f_one(founders, offspring_per_pair)
   >>> mrc = breed.MultiRandomCross(multi_prefounders, 4, offspring_per_pair)
   >>> mother_choices, father_choices = mrc.determine_random_cross()
   >>> multi_snd_ord_chooser = breed.MultiSecondOrderPairIDChooser(mother_choices, father_choices)
   >>> multi_prefounders.evolve(
   ...        matingScheme=sim.HomoMating(
   ...            sim.PyParentsChooser(multi_snd_ord_chooser.snd_ord_id_pairs),
   ...            sim.OffspringGenerator(ops=[
   ...                sim.IdTagger(),
   ...                sim.PedigreeTagger(),
   ...                sim.Recombinator(rates=recombination_rates)
   ...            ],
   ...                numOffspring=1),
   ...            subPopSize=[operating_population_size],
   ...        ),
   ...        gen=1,
   ...    )

   >>> final_mrc = breed.MultiRandomCross(multi_prefounders, 2, 1000)
   >>> final_mothers, final_fathers = final_mrc.determine_random_cross()
   >>> final_multi_snd_ord_chooser = breed.MultiSecondOrderPairIDChooser(final_mothers, final_fathers)
   >>> multi_prefounders.evolve(
   ...      matingScheme=sim.HomoMating(
   ...            sim.PyParentsChooser(final_multi_snd_ord_chooser.snd_ord_id_pairs),
   ...            sim.OffspringGenerator(ops=[
   ...                sim.IdTagger(),
   ...                sim.PedigreeTagger(),
   ...                sim.Recombinator(rates=0.01)
   ...            ],
   ...                numOffspring=1),
   ...            subPopSize=[operating_population_size],
   ...        ),
   ...        gen=1,
   ...    )
   >>> mater = breed.MAGIC(multi_prefounders, founders, recombination_rates)
   >>> mater.random_mating(3, operating_population_size)
   Initiating random mating for 3 generations.
    
.. _sample_collection:

Sample Collection
=================

After our simulation finishes we collect all of the samples at once and store
them in a dictionary keyed by replicate ID. You can see below that each
entry of ``sample_library`` is a list of populations.

.. code-block:: python

   >>> demonstration = analyze.Study(run_id)
   >>> sample_library = demonstration.collect_samples(multi_prefounders, sample_sizes)
   >>> sample_library
    {0: [<simuPOP.Population>,
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
      <simuPOP.Population>],
     2: [<simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>],
     3: [<simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>],
     4: [<simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>,
      <simuPOP.Population>]}



    >>> alleles = np.array(pd.read_hdf('parameters\\alleles_at_1478_loci.hdf'))
    >>> alleles
    array([[1, 2],
           [1, 3],
           [3, 1],
           ..., 
           [1, 0],
           [3, 0],
           [3, 1]], dtype=int64)


.. _storing_allele_frequencies:

Allele Frequency
================

Allele frequency data is stored in an hdf5 file. Allele frequencies are collected
from all replicates and all samples from that replicate.

.. code-block:: python
   :caption: Storing allele frequencies

   >>> analyze.store_allele_frequencies(sample_library, alleles, run_id)

Organization of Allele Frequency Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example of a single allele frequency table for ``replicate 0`` at
``sample size 500`` for the ``run_id demonstration``. The table is indexed by
``absolute locus index``.


.. code-block:: python
   :caption: The structure of an allele frequency table

   >>> reloaded_allele_frequencies = pd.read_hdf('demonstration_allele_frequency_storage.h5')

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
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
          <td>2.0</td>
          <td>0.000</td>
          <td>1.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>3.0</td>
          <td>0.141</td>
          <td>1.0</td>
          <td>0.859</td>
        </tr>
        <tr>
          <th>2</th>
          <td>1.0</td>
          <td>0.125</td>
          <td>3.0</td>
          <td>0.875</td>
        </tr>
        <tr>
          <th>3</th>
          <td>2.0</td>
          <td>0.105</td>
          <td>0.0</td>
          <td>0.895</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.0</td>
          <td>0.043</td>
          <td>2.0</td>
          <td>0.957</td>
        </tr>
        <tr>
          <th>5</th>
          <td>2.0</td>
          <td>0.219</td>
          <td>0.0</td>
          <td>0.781</td>
        </tr>
        <tr>
          <th>6</th>
          <td>2.0</td>
          <td>0.272</td>
          <td>0.0</td>
          <td>0.728</td>
        </tr>
        <tr>
          <th>7</th>
          <td>1.0</td>
          <td>0.000</td>
          <td>3.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>8</th>
          <td>2.0</td>
          <td>0.089</td>
          <td>0.0</td>
          <td>0.911</td>
        </tr>
        <tr>
          <th>9</th>
          <td>3.0</td>
          <td>0.000</td>
          <td>1.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>10</th>
          <td>3.0</td>
          <td>0.413</td>
          <td>1.0</td>
          <td>0.587</td>
        </tr>
        <tr>
          <th>11</th>
          <td>3.0</td>
          <td>0.000</td>
          <td>1.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>12</th>
          <td>1.0</td>
          <td>0.366</td>
          <td>3.0</td>
          <td>0.634</td>
        </tr>
        <tr>
          <th>13</th>
          <td>0.0</td>
          <td>0.090</td>
          <td>2.0</td>
          <td>0.910</td>
        </tr>
        <tr>
          <th>14</th>
          <td>0.0</td>
          <td>0.128</td>
          <td>3.0</td>
          <td>0.872</td>
        </tr>
        <tr>
          <th>15</th>
          <td>1.0</td>
          <td>0.401</td>
          <td>3.0</td>
          <td>0.599</td>
        </tr>
        <tr>
          <th>16</th>
          <td>3.0</td>
          <td>0.130</td>
          <td>2.0</td>
          <td>0.870</td>
        </tr>
        <tr>
          <th>17</th>
          <td>3.0</td>
          <td>0.000</td>
          <td>2.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>18</th>
          <td>0.0</td>
          <td>0.000</td>
          <td>2.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>19</th>
          <td>4.0</td>
          <td>0.000</td>
          <td>5.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>20</th>
          <td>0.0</td>
          <td>0.388</td>
          <td>3.0</td>
          <td>0.612</td>
        </tr>
        <tr>
          <th>21</th>
          <td>1.0</td>
          <td>0.123</td>
          <td>2.0</td>
          <td>0.877</td>
        </tr>
        <tr>
          <th>22</th>
          <td>3.0</td>
          <td>0.000</td>
          <td>2.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>23</th>
          <td>5.0</td>
          <td>0.437</td>
          <td>4.0</td>
          <td>0.563</td>
        </tr>
        <tr>
          <th>24</th>
          <td>3.0</td>
          <td>0.061</td>
          <td>1.0</td>
          <td>0.939</td>
        </tr>
        <tr>
          <th>25</th>
          <td>0.0</td>
          <td>0.244</td>
          <td>2.0</td>
          <td>0.756</td>
        </tr>
        <tr>
          <th>26</th>
          <td>0.0</td>
          <td>0.075</td>
          <td>3.0</td>
          <td>0.925</td>
        </tr>
        <tr>
          <th>27</th>
          <td>1.0</td>
          <td>0.058</td>
          <td>0.0</td>
          <td>0.942</td>
        </tr>
        <tr>
          <th>28</th>
          <td>1.0</td>
          <td>0.000</td>
          <td>2.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>29</th>
          <td>0.0</td>
          <td>0.368</td>
          <td>2.0</td>
          <td>0.632</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>1448</th>
          <td>0.0</td>
          <td>0.244</td>
          <td>3.0</td>
          <td>0.756</td>
        </tr>
        <tr>
          <th>1449</th>
          <td>3.0</td>
          <td>0.395</td>
          <td>0.0</td>
          <td>0.605</td>
        </tr>
        <tr>
          <th>1450</th>
          <td>2.0</td>
          <td>0.219</td>
          <td>0.0</td>
          <td>0.781</td>
        </tr>
        <tr>
          <th>1451</th>
          <td>3.0</td>
          <td>0.000</td>
          <td>1.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1452</th>
          <td>5.0</td>
          <td>0.000</td>
          <td>4.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1453</th>
          <td>3.0</td>
          <td>0.000</td>
          <td>1.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1454</th>
          <td>1.0</td>
          <td>0.107</td>
          <td>3.0</td>
          <td>0.893</td>
        </tr>
        <tr>
          <th>1455</th>
          <td>2.0</td>
          <td>0.467</td>
          <td>3.0</td>
          <td>0.533</td>
        </tr>
        <tr>
          <th>1456</th>
          <td>2.0</td>
          <td>0.141</td>
          <td>0.0</td>
          <td>0.859</td>
        </tr>
        <tr>
          <th>1457</th>
          <td>0.0</td>
          <td>0.036</td>
          <td>2.0</td>
          <td>0.964</td>
        </tr>
        <tr>
          <th>1458</th>
          <td>1.0</td>
          <td>0.470</td>
          <td>0.0</td>
          <td>0.530</td>
        </tr>
        <tr>
          <th>1459</th>
          <td>5.0</td>
          <td>0.000</td>
          <td>4.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1460</th>
          <td>0.0</td>
          <td>0.146</td>
          <td>3.0</td>
          <td>0.854</td>
        </tr>
        <tr>
          <th>1461</th>
          <td>1.0</td>
          <td>0.054</td>
          <td>2.0</td>
          <td>0.946</td>
        </tr>
        <tr>
          <th>1462</th>
          <td>3.0</td>
          <td>0.102</td>
          <td>1.0</td>
          <td>0.898</td>
        </tr>
        <tr>
          <th>1463</th>
          <td>3.0</td>
          <td>0.263</td>
          <td>1.0</td>
          <td>0.737</td>
        </tr>
        <tr>
          <th>1464</th>
          <td>1.0</td>
          <td>0.000</td>
          <td>3.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1465</th>
          <td>0.0</td>
          <td>0.147</td>
          <td>2.0</td>
          <td>0.853</td>
        </tr>
        <tr>
          <th>1466</th>
          <td>3.0</td>
          <td>0.000</td>
          <td>1.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1467</th>
          <td>0.0</td>
          <td>0.000</td>
          <td>2.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1468</th>
          <td>0.0</td>
          <td>0.326</td>
          <td>2.0</td>
          <td>0.674</td>
        </tr>
        <tr>
          <th>1469</th>
          <td>2.0</td>
          <td>0.000</td>
          <td>1.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1470</th>
          <td>0.0</td>
          <td>0.219</td>
          <td>2.0</td>
          <td>0.781</td>
        </tr>
        <tr>
          <th>1471</th>
          <td>3.0</td>
          <td>0.000</td>
          <td>0.0</td>
          <td>1.000</td>
        </tr>
        <tr>
          <th>1472</th>
          <td>0.0</td>
          <td>0.103</td>
          <td>1.0</td>
          <td>0.897</td>
        </tr>
        <tr>
          <th>1473</th>
          <td>1.0</td>
          <td>0.049</td>
          <td>0.0</td>
          <td>0.951</td>
        </tr>
        <tr>
          <th>1474</th>
          <td>2.0</td>
          <td>0.237</td>
          <td>0.0</td>
          <td>0.763</td>
        </tr>
        <tr>
          <th>1475</th>
          <td>0.0</td>
          <td>0.236</td>
          <td>1.0</td>
          <td>0.764</td>
        </tr>
        <tr>
          <th>1476</th>
          <td>0.0</td>
          <td>0.137</td>
          <td>3.0</td>
          <td>0.863</td>
        </tr>
        <tr>
          <th>1477</th>
          <td>1.0</td>
          <td>0.000</td>
          <td>3.0</td>
          <td>1.000</td>
        </tr>
      </tbody>
    </table>
    <p>1478 rows × 4 columns</p>
    </div>


Allele Effects and Segregating Loci
===================================

If at all possible we would like to have a common set of segregating loci across
all replicates and all samples. All the samples have been collected into the
variable ``sample_library``. So we can collect all segregating loci of all samples
and examine if there are any differences. If there is more than one value then
there is more than one set of segregating loci.

.. code-block:: python
   :capation: Determine segregating loci in all samples

   >>> sets_of_segregating_loci = demonstration.seg_loci_among_samples(sample_library)
   >>> sets_of_segregating_loci.values()
   dict_values([30])

   >>> concordant_segregating_loci = list(sets_of_segregating_loci.keys())[0]

Determine QTL From Concordant Set of Segregating Loci
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the time being I am assigning QTL only to loci which are segregating. Assigning
QTL to segregating loci is not realistic biologically; however, it simplifies
the analysis. Later versions will assign effects to prefounder populations.

.. code-block:: python
   :caption: Choosing QTL

   >>> qtl = sorted(random.sample(concordant_segregating_loci, number_of_qtl))
   >>> qtl
   [246, 432, 527, 783, 965, 998, 1035, 1056, 1245, 1444]

Assigning Additive Allele Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example of assigning allele effects via random draws from an exponential
distribution with lambda = 1. The py:func:`random.expovariate` is an argument
in the function. We can change the allele effect distribution by replacing
py:func:`random.expovariate` with whatever distribution function we like and its
parameters (multiple parameter functions work as well).

.. code-block::
   :caption: Assigning additive allele effects

   >>> add_trait = parameters.Trait()
   >>> allele_effects = add_trait.assign_allele_effects(alleles, qtl, random.expovariate, 1, multiplicity=3)
   >>> allele_effects
   {246: {1: 3.13370150370361, 3: 2.3333776978977627},
   432: {0: 3.307659276528477, 3: 2.3923475464249715},
   527: {0: 1.8558917885028081, 2: 1.5406900580075562},
   783: {4: 2.317132355134784, 5: 1.3295667375269518},
   965: {0: 2.31035019629015, 2: 6.22957905138777},
   998: {1: 1.1739532295469035, 3: 1.2072378820811571},
   1035: {1: 4.493406487495378, 3: 1.1529343427499426},
   1056: {4: 1.8568520871689185, 5: 5.06545115412201},
   1245: {2: 3.458945179018148, 3: 1.5286068388242993},
   1444: {4: 4.227937082576118, 5: 3.3236868346837367}}

Store Metadata and Useful Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is probably a good idea to store the allele effects under the
run id ``demonstration``.

.. code-block:: python
   :caption: Storing allele effects

   >>> allele_effects_store = shelve.open('allele_effects_storage')
   >>> allele_effects_store['demonstration'] = allele_effects
   >>> allele_effects_store.close()

Because we drop the non-segregating loci and relabel the segregating loci before
we run TASSEL on our data we need to be able to convert TASSEL loci to their
original ``saegus`` loci.

.. code-block:: python
   :caption: Storing forward and backward maps

   >>> loci_conversions = shelve.open('demonstration_loci_conversions')
   >>> saegus_to_tassel_loci = {}
   >>> tassel_to_saegus_loci = {}
   >>> for idx, locus in enumerate(concordant_segregating_loci):
   ...  saegus_to_tassel_loci[locus] = idx
   ...  tassel_to_saegus_loci[idx] = locus
   >>> loci_conversions['saegus_to_tassel'] = saegus_to_tassel_loci
   >>> loci_conversions['tassel_to_saegus'] = tassel_to_saegus_loci
   >>> loci_conversions.close()

Formatting Data in Human Readable Tables
========================================

A ``python`` :py:class:`dict` is not a human friendly way of presenting data.


analyze.store_allele_effect_frequency_tables(sample_library, alleles,
                                             qtl,
                                             exponential_allele_effects,
                                            run_id, 'exponential')


.. code-block:: python
   :caption: Allele effect table.

   >>> analyze.store_allele_effect_frequency_tables(sample_library, alleles, qtl, allele_effects, run_id, 'exponential')

.. raw:: html

   <div>
   <table border="1" class="dataframe">
   <thead>
     <tr style="text-align: right;">
       <th></th><th>locus</th>
       <th>tassel_locus</th>
       <th>alpha_allele</th>
       <th>alpha_effect</th>
       <th>beta_allele</th>
       <th>beta_effect</th>
       <th>difference</th>
     </tr>
   </thead>
   <tbody>
     <tr>
       <th>0</th>
       <td>246</td>
       <td>141</td>
       <td>1</td>
       <td>3.133702</td>
       <td>3</td>
       <td>2.333378</td>
       <td>0.800324</td>
     </tr>
     <tr>
       <th>1</th>
       <td>432</td>
       <td>263</td>
       <td>0</td>
       <td>3.307659</td>
       <td>3</td>
       <td>2.392348</td>
       <td>0.915312</td>
     </tr>
     <tr>
       <th>2</th>
       <td>527</td>
       <td>333</td>
       <td>2</td>
       <td>1.540690</td>
       <td>0</td>
       <td>1.855892</td>
       <td>0.315202</td>
     </tr>
     <tr>
       <th>3</th>
       <td>783</td>
       <td>498</td>
       <td>4</td>
       <td>2.317132</td>
       <td>5</td>
       <td>1.329567</td>
       <td>0.987566</td>
     </tr>
     <tr>
       <th>4</th>
       <td>965</td>
       <td>611</td>
       <td>2</td>
       <td>6.229579</td>
       <td>0</td>
       <td>2.310350</td>
       <td>3.919229</td>
     </tr>
     <tr>
       <th>5</th>
       <td>998</td>
       <td>632</td>
       <td>1</td>
       <td>1.173953</td>
       <td>3</td>
       <td>1.207238</td>
       <td>0.033285</td>
     </tr>
     <tr>
       <th>6</th>
       <td>1035</td>
       <td>662</td>
       <td>3</td>
       <td>1.152934</td>
       <td>1</td>
       <td>4.493406</td>
       <td>3.340472</td>
     </tr>
     <tr>
       <th>7</th>
       <td>1056</td>
       <td>675</td>
       <td>4</td>
       <td>1.856852</td>
       <td>5</td>
       <td>5.065451</td>
       <td>3.208599</td>
     </tr>
     <tr>
       <th>8</th>
       <td>1245</td>
       <td>794</td>
       <td>2</td>
       <td>3.458945</td>
       <td>3</td>
       <td>1.528607</td>
       <td>1.930338</td>
     </tr>
     <tr>
       <th>9</th>
       <td>1444</td>
       <td>919</td>
       <td>5</td>
       <td>3.323687</td>
       <td>4</td>
       <td>4.227937</td>
       <td>0.904250</td>
     </tr>
   </tbody>
   </table>
   </div>

.. note::

   If you need the chromosome and the relative locus on that chromosome we can use
   a simuPOP function called :py:func:`Population.chromLocPair`()`

Computing TASSEL Input and Writing Input Files
==============================================

All of the calculations or manipulations required TASSEL input are handled by
the single function :py:func:`write_multiple_sample_analyzer()`.


.. code-block:: python
   :caption: TASSEL input

   >>> analyze.write_multiple_sample_analyzer(sample_library, sample_sizes, qtl, alleles,
   ... allele_effects, 0.7,  concordant_segregating_loci,
   ... run_id='demonstration',
   ... allele_frequency_hdf='demonstration_storage.h5')

The output by default is set to go into ``C:\\tassel\\input``. Work in progress
to be able to dynamically and easily specify the routing of input and output.

.. note::

   Intermission of running TASSEL. Running TASSEL is documented elsewhere.

Analyzing TASSEL Output
=======================

After TASSEL finishes running we have the results deposited into
``C:\\tassel\\output``. We need to calculate the ``q-values`` for the TASSEL
output. The ``q-value`` calculation is performed in R.

.. note::

   Eventually I will convert the ``q-value`` calculation into Python.

.. code-block:: rconsole

   > library(qvalue)
   > library(ggplot2)
   > library(gap)
   > setwd("C:/tassel/output")

   > sample_sizes = seq(from = 250, to = 1000, by = 250)
   > replicates = seq(from = 0, to = 49, by = 250)
   > base_file_prefix = "demonstration_"
   > base_file_suffix = "_out_2.txt"
   > base_qvalue_file_prefix = "demonstration_"
   > qvalue_file_suffix = "_qvalues.txt"

   > for(rep in replicates){
      for(sample in sample_sizes){

         input_file_name = paste(base_file_prefix, rep, sep="", '_', sample, base_file_suffix)
         output_file_name = paste(base_file_prefix, rep, sep="", '_', sample, q_value_suffix)

         results_header = scan(input_file_name, what="character", nlines=1, sep="\t")
         gwas_results = read.table(input_file_name, header=F, row.names=NULL, skip=2)
         colnames(gwas_results) = results_header

         pvalues = gwas_results$p

         qobj = qvalue(p = pvalues)
         qvalues = data.frame(qobj$qvalues)
         write.table(qvalues, output_file_name, sep="\t")
      }

   }

After we have our ``q-values`` we collect the data into super tables which
is a large table which contains all the of the data we have computed so far.




.. code-block:: python
   :caption: Collection of power and false positive rate data into a dictionary

   >>> demonstration = analyze.Study(run_id)
   >>> power_fpr_raw_data = demonstration.collect_power_analysis_data(sample_sizes, number_of_replicates, expanded)
   >>> power_fpr_raw_data

.. parsed-literal::
   :caption: Example of TASSEL output + q-values + allele effect difference

   {500:
   {0:       Chr  df         F        p        q  difference
   0     1.0   2  0.198370  0.82013  0.99982     0.00000
   1     1.0   2  0.507910  0.60207  0.99982     0.00000
   2     1.0   2  1.032440  0.35690  0.99982     0.00000
   3     1.0   1  0.141970  0.70649  0.99982     0.00000
   4     1.0   2  0.609030  0.54429  0.99982     0.00000
   5     1.0   2  0.677510  0.50835  0.99982     0.00000
   6     1.0   2  1.128290  0.32442  0.99982     0.00000
   7     1.0   2  0.933870  0.39372  0.99982     0.00000
   8     1.0   2  1.131010  0.32354  0.99982     0.00000
   9     1.0   2  0.026590  0.97376  0.99982     0.00000
   10    1.0   2  0.113300  0.89291  0.99982     0.00000
   11    1.0   2  0.770060  0.46354  0.99982     0.00000
   12    1.0   2  0.258250  0.77251  0.99982     0.00000
   13    1.0   2  0.482280  0.61767  0.99982     0.00000
   14    1.0   2  0.883250  0.41409  0.99982     0.00000
   15    1.0   2  0.733070  0.48095  0.99982     0.00000
   16    1.0   2  1.523410  0.21899  0.99982     0.00000
   17    1.0   2  0.150900  0.85998  0.99982     0.00000
   18    1.0   2  0.174470  0.83995  0.99982     0.00000
   19    1.0   2  0.134980  0.87376  0.99982     0.00000
   20    1.0   2  0.361180  0.69703  0.99982     0.00000
   21    1.0   2  0.082190  0.92111  0.99982     0.00000
   22    1.0   2  0.018710  0.98146  0.99982     0.00000
   23    1.0   2  2.435070  0.08864  0.99982     0.00000
   24    1.0   2  1.108880  0.33075  0.99982     0.00000
   25    1.0   2  0.517210  0.59650  0.99982     0.00000
   26    1.0   2  1.188360  0.30559  0.99982     0.00000
   27    1.0   2  0.039310  0.96145  0.99982     0.00000
   28    1.0   2  0.190020  0.82700  0.99982     0.00000
   29    1.0   2  1.908650  0.14937  0.99982     0.00000
   ...    ...  ...       ...      ...      ...         ...
   913  10.0   2  0.774170  0.46165  0.99982     0.00000
   914  10.0   2  0.923280  0.39790  0.99982     0.00000
   915  10.0   2  0.793190  0.45297  0.99982     0.00000
   916  10.0   2  2.949290  0.05330  0.99982     0.00000
   917  10.0   2  0.138340  0.87084  0.99982     0.00000
   918  10.0   2  2.770770  0.06359  0.99982     0.00000
   919  10.0   2  1.052820  0.34973  0.99982     0.90425
   920  10.0   2  0.778770  0.45953  0.99982     0.00000
   921  10.0   2  0.762870  0.46687  0.99982     0.00000
   922  10.0   2  0.196290  0.82184  0.99982     0.00000
   923  10.0   2  1.057310  0.34817  0.99982     0.00000
   924  10.0   2  1.394370  0.24896  0.99982     0.00000
   925  10.0   2  0.186340  0.83005  0.99982     0.00000
   926  10.0   2  1.116590  0.32822  0.99982     0.00000
   927  10.0   2  0.216190  0.80566  0.99982     0.00000
   928  10.0   2  0.285520  0.75175  0.99982     0.00000
   929  10.0   1  0.000109  0.99166  0.99982     0.00000
   930  10.0   2  0.319140  0.72692  0.99982     0.00000
   931  10.0   2  0.751450  0.47222  0.99982     0.00000
   932  10.0   2  2.080780  0.12592  0.99982     0.00000
   933  10.0   2  1.988760  0.13796  0.99982     0.00000
   934  10.0   2  0.472400  0.62379  0.99982     0.00000
   935  10.0   2  0.610200  0.54365  0.99982     0.00000
   936  10.0   2  1.887420  0.15255  0.99982     0.00000
   937  10.0   2  1.459130  0.23344  0.99982     0.00000
   938  10.0   2  1.366180  0.25604  0.99982     0.00000
   939  10.0   1  3.595660  0.05851  0.99982     0.00000
   940  10.0   2  0.971010  0.37942  0.99982     0.00000
   941  10.0   2  2.079250  0.12611  0.99982     0.00000
   942  10.0   2  1.127040  0.32482  0.99982     0.00000

   [943 rows x 6 columns],
   1:       Chr  df        F        p         q  difference
   0     1.0   2  1.73693  0.17713  0.990741     0.00000
   1     1.0   2  2.09598  0.12404  0.990741     0.00000
   2     1.0   2  0.67511  0.50957  0.990741     0.00000
   3     1.0   2  0.04883  0.95235  0.991774     0.00000
   4     1.0   2  0.32188  0.72494  0.990741     0.00000
   5     1.0   2  0.33153  0.71799  0.990741     0.00000
   6     1.0   2  0.57059  0.56556  0.990741     0.00000
   7     1.0   2  0.34987  0.70496  0.990741     0.00000
   8     1.0   2  2.76706  0.06382  0.990741     0.00000
   9     1.0   2  0.30501  0.73725  0.990741     0.00000
   10    1.0   2  0.00240  0.99760  0.999900     0.00000
   11    1.0   2  0.35998  0.69787  0.990741     0.00000
   12    1.0   2  0.22766  0.79648  0.990741     0.00000
   13    1.0   2  1.00135  0.36813  0.990741     0.00000
   14    1.0   2  0.31618  0.72908  0.990741     0.00000
   15    1.0   2  1.17775  0.30883  0.990741     0.00000
   16    1.0   2  0.83942  0.43257  0.990741     0.00000
   17    1.0   2  1.57616  0.20780  0.990741     0.00000
   18    1.0   2  2.60308  0.07506  0.990741     0.00000
   19    1.0   2  0.50234  0.60542  0.990741     0.00000
   20    1.0   2  0.46813  0.62645  0.990741     0.00000
   21    1.0   2  0.29473  0.74486  0.990741     0.00000
   22    1.0   2  0.31265  0.73165  0.990741     0.00000
   23    1.0   2  1.52864  0.21785  0.990741     0.00000
   24    1.0   2  0.15071  0.86014  0.990741     0.00000
   25    1.0   2  1.50446  0.22315  0.990741     0.00000
   26    1.0   2  1.69603  0.18447  0.990741     0.00000
   27    1.0   2  0.69291  0.50060  0.990741     0.00000
   28    1.0   2  0.09875  0.90598  0.990741     0.00000
   29    1.0   2  0.93787  0.39216  0.990741     0.00000
   ...    ...  ...      ...      ...       ...         ...
   913  10.0   2  3.22115  0.04075  0.990741     0.00000
   914  10.0   2  0.88494  0.41339  0.990741     0.00000
   915  10.0   2  0.36133  0.69693  0.990741     0.00000
   916  10.0   2  0.99232  0.37145  0.990741     0.00000
   917  10.0   2  0.19015  0.82690  0.990741     0.00000
   918  10.0   2  1.84142  0.15968  0.990741     0.00000
   919  10.0   2  5.64215  0.00378  0.356164     0.90425
   920  10.0   2  0.04951  0.95170  0.991774     0.00000
   921  10.0   2  0.02037  0.97984  0.996577     0.00000
   922  10.0   2  0.60646  0.54568  0.990741     0.00000
   923  10.0   2  1.17492  0.30970  0.990741     0.00000
   924  10.0   2  0.15283  0.85832  0.990741     0.00000
   925  10.0   2  0.09760  0.90703  0.990741     0.00000
   926  10.0   2  2.00713  0.13547  0.990741     0.00000
   927  10.0   2  1.06895  0.34416  0.990741     0.00000
   928  10.0   2  0.10252  0.90258  0.990741     0.00000
   929  10.0   2  0.09430  0.91002  0.990741     0.00000
   930  10.0   2  0.65323  0.52081  0.990741     0.00000
   931  10.0   2  0.29747  0.74283  0.990741     0.00000
   932  10.0   2  1.06281  0.34627  0.990741     0.00000
   933  10.0   2  0.17255  0.84157  0.990741     0.00000
   934  10.0   2  0.06837  0.93392  0.990741     0.00000
   935  10.0   2  0.40211  0.66912  0.990741     0.00000
   936  10.0   2  0.61450  0.54132  0.990741     0.00000
   937  10.0   2  0.42815  0.65195  0.990741     0.00000
   938  10.0   2  0.14843  0.86210  0.990741     0.00000
   939  10.0   2  1.96086  0.14183  0.990741     0.00000
   940  10.0   2  0.81085  0.44507  0.990741     0.00000
   941  10.0   2  0.63080  0.53259  0.990741     0.00000
   942  10.0   2  0.28102  0.75514  0.990741     0.00000


Finally we can calculate the power and false positive rate.

.. code-block:: python
   :caption: Calculating power from the results.

   >>> results, true_positives, false_positives = demonstration.calculate_power_fpr(power_fpr_raw_data,
   ...                                              sample_sizes, number_of_replicates,
   ...                                              number_of_replicates,
   ...                                              number_of_qtl)
   results

.. raw:: html

   <div>
   <table border="1" class="dataframe">
   <thead>
     <tr style="text-align: right;">
       <th></th>
       <th>power_500</th>
       <th>fpr_500</th>
       <th>power_600</th>
       <th>fpr_600</th>
       <th>power_700</th>
       <th>fpr_700</th>
       <th>power_800</th>
       <th>fpr_800</th>
       <th>power_900</th>
       <th>fpr_900</th>
       <th>power_1000</th>
       <th>fpr_1000</th>
     </tr>
   </thead>
   <tbody>
     <tr>
       <th>0</th>
       <td>0.4</td>
       <td>0</td>
       <td>0.4</td>
       <td>0</td>
       <td>0.4</td>
       <td>0</td>
       <td>0.4</td>
       <td>0</td>
       <td>0.4</td>
       <td>0</td>
       <td>0.4</td>
       <td>0.00107181</td>
     </tr>
     <tr>
       <th>1</th>
       <td>0.3</td>
       <td>0.00107181</td>
       <td>0.4</td>
       <td>0.00107181</td>
       <td>0.3</td>
       <td>0</td>
       <td>0.4</td>
       <td>0.00107181</td>
       <td>0.4</td>
       <td>0.00214362</td>
       <td>0.4</td>
       <td>0.00107181</td>
     </tr>
     <tr>
       <th>2</th>
       <td>0.5</td>
       <td>0</td>
       <td>0.4</td>
       <td>0.00214362</td>
       <td>0.4</td>
       <td>0.00321543</td>
       <td>0.4</td>
       <td>0.00214362</td>
       <td>0.4</td>
       <td>0.00107181</td>
       <td>0.4</td>
       <td>0.00107181</td>
     </tr>
     <tr>
       <th>3</th>
       <td>0.5</td>
       <td>0.00107181</td>
       <td>0.3</td>
       <td>0.00107181</td>
       <td>0.4</td>
       <td>0</td>
       <td>0.4</td>
       <td>0.00214362</td>
       <td>0.4</td>
       <td>0.00321543</td>
       <td>0.4</td>
       <td>0.00107181</td>
     </tr>
     <tr>
       <th>4</th>
       <td>0.4</td>
       <td>0</td>
       <td>0.4</td>
       <td>0</td>
       <td>0.4</td>
       <td>0.00214362</td>
       <td>0.5</td>
       <td>0.00107181</td>
       <td>0.4</td>
       <td>0.00107181</td>
       <td>0.4</td>
       <td>0.00214362</td>
     </tr>
   </tbody>
   </table>
   </div>

And the true positive results.

.. code-block:: python

   >>> true_positives
   {(500, 0):      Chr   df         F             p             q  difference
   611  7.0  2.0  27.67028  4.059600e-12  3.828203e-09    3.919229
   662  7.0  2.0  23.63099  1.577700e-10  7.438855e-08    3.340472
   675  7.0  2.0  19.24244  8.955000e-09  2.814855e-06    3.208599
   794  9.0  2.0   9.94932  5.803000e-05  1.368057e-02    1.930338,
   (500, 1):      Chr   df         F             p             q  difference
   611  7.0  2.0  27.58324  4.390300e-12  4.140053e-09    3.919229
   662  7.0  2.0  13.73347  1.567600e-06  4.927489e-04    3.340472
   675  7.0  2.0  21.73976  8.921000e-10  4.206252e-07    3.208599,
   (500, 2):       Chr   df         F             p         q  difference
   611   7.0  1.0  36.84348  2.545700e-09  0.000002    3.919229
   662   7.0  2.0  16.64202  1.011800e-07  0.000047    3.340472
   675   7.0  2.0  10.48535  3.468100e-05  0.006468    3.208599
   794   9.0  2.0  13.96430  1.259800e-06  0.000392    1.930338
   919  10.0  2.0  12.83702  3.670500e-06  0.000856    0.904250,
   (500, 3):      Chr   df         F             p             q  difference
   263  3.0  2.0   8.35635  2.696500e-04  4.304795e-02    0.915312
   611  7.0  2.0  17.13919  6.352700e-08  2.995298e-05    3.919229
   662  7.0  2.0  26.81431  8.777500e-12  8.277183e-09    3.340472
   675  7.0  2.0   9.89850  6.093600e-05  1.436566e-02    3.208599
   794  9.0  2.0  10.07199  5.157700e-05  1.436566e-02    1.930338,
   (500, 4):      Chr   df         F             p             q  difference
   611  7.0  2.0  22.34613  5.112200e-10  4.458647e-07    3.919229
   662  7.0  2.0  17.00027  7.234300e-08  2.103151e-05    3.340472
   675  7.0  2.0  17.11715  6.485000e-08  2.103151e-05    3.208599
   794  9.0  2.0   8.57896  2.174300e-04  4.740834e-02    1.930338,
   (600, 0):      Chr   df         F             p             q  difference
   611  7.0  2.0  29.43589  6.452900e-13  6.085085e-10    3.919229
   662  7.0  2.0  17.15029  5.734600e-08  2.703864e-05    3.340472
   675  7.0  2.0  10.56975  3.085200e-05  7.273359e-03    3.208599
   794  9.0  2.0  12.54095  4.626000e-06  1.454106e-03    1.930338,
   (600, 1):      Chr   df         F             p             q  difference
   611  7.0  2.0  23.25339  1.889800e-10  1.782081e-07    3.919229
   662  7.0  2.0  16.10681  1.540600e-07  4.842619e-05    3.340472
   675  7.0  2.0  21.64993  8.393300e-10  3.957441e-07    3.208599
   794  9.0  2.0  10.31065  3.962700e-05  7.473652e-03    1.930338,

