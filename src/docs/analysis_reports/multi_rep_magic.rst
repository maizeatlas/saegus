.. _multi_rep_magic:

============================================
Creating Multiple Replicate MAGIC Population
============================================

I have defined functions to handle making multiple replicates of a MAGIC
population from a common group of prefounders. This is a record of the exact
steps I used to make a multiple replicate population from 8 prefounders.

Imports
=======

.. code-block:: python

   import simuOpt
   simuOpt.setOptions(alleleType='short', numThreads=4, quiet=True)
   import simuPOP as sim
   import pandas as pd
   from saegus import breed, operators, simulate, analyze, parse, parameters
   import shelve
   import numpy as np
   import random
   np.set_printoptions(suppress=True, precision=3)

Determine Founders
==================

I used prefounders 1-8 from the NAM prefounders. We are using
the 1478 integer-valued positions of the genetic map. I created a
population of 2000 individuals in the first step and kept the population
size constant throughout all subsequent steps. So 500 offspring per pair
of individuals and 4 pairs results in 2000 offspring.

.. code-block:: python

   >>> prefounders = sim.loadPopulation('prefounders1478.pop')
   >>> multi_prefounders = sim.Simulator(prefounders, 50, stealPops=False)
   >>> founders = [[1, 2], [3, 4], [5, 6], [7, 8]]
   >>> os_per_pair = 500

The regular :class:`MAGIC` works with simuPOP.Simulator objects so no
modification was needed for this step. simuPOP.tagID(prefounders, reset=27)
resets the system-wide index which determines which number to
start ``ind_id`` for each generation.

.. code-block:: python

   >>> magic = breed.MAGIC(multi_prefounders, [0.01]*1478)
   >>> sim.tagID(prefounders, reset=27)
   >>> magic.generate_f_one(founders, os_per_pair)
   Generation: 0
   Generation: 0
   .....
   Generation: 0

.. _second_step_magic:

Second Step of MAGIC
====================

The first round of mating is unique in that it is very simple. For
the subsequent rounds of mating I determine which individuals cross
for each replicate before mating occurs. The pre-determined crosses are
created by randomly sampling individuals from each sub-population of hybrids.
That is:

   Hybrid sub-population 1 is made from founder pair [1, 2]
   Hybrid sub-population 2 is made from founder pair [3, 4]
   Hybrid sub-population 1 is made from founder pair [5, 6]
   Hybrid sub-population 1 is made from founder pair [7, 8]

simuPOP maintains this order in generating the offspring so the first 500 individuals
are offspring of [1, 2]. The next 500 individuals are offspring of [3, 4] so on
and so forth.

Hence I split the first hybrid offspring population in 4 sub-populations which
each have size 500. From each sub-population we sample twice the number of individuals
in that sub-population so 1000 individuals. This process is completed with
a new class :py:class:`MultiSecondOrderPairIDChooser`

.. code-block:: python

   >>> mrc = breed.MultiRandomCross(multi_prefounders, 4, 500)
   >>> mother_choices, father_choices = mrc.determine_random_cross()
   >>> multi_snd_ord_chooser = breed.MultiSecondOrderPairIDChooser(mother_choices, father_choices)
   >>> multi_prefounders.evolve(
   ...     matingScheme=sim.HomoMating(
   ...         sim.PyParentsChooser(multi_snd_ord_chooser.snd_ord_id_pairs),
   ...         sim.OffspringGenerator(ops=[
   ...             sim.IdTagger(),
   ...             sim.PedigreeTagger(),
   ...             sim.Recombinator(rates=0.01)
   ...          ],
   ...              numOffspring=1),
   ...          subPopSize=[2000],
   ...      ),
   ...      gen=1,
   ...  )
   (1,
   1,
   1,
   1,
   ........
   1)

.. _final_step_magic:

Final Step of MAGIC
===================

The third and final step is nearly the same as the `second_step_magic`_
In the final step we split the population into 2 sub-populations with 1000 individuals
each. Our new function will randomly sample twice the size of the sub-population size
so 2000 individuals are sampled from each sub-population.

.. code-block:: python

   >>> final_mrc = breed.MultiRandomCross(multi_prefounders, 2, 1000)
   >>> final_mothers, final_fathers = final_mrc.determine_random_cross()
   >>> final_multi_snd_ord_chooser = breed.MultiSecondOrderPairIDChooser(final_mothers, final_fathers)
   >>> multi_prefounders.evolve(
   ...     matingScheme=sim.HomoMating(
   ...         sim.PyParentsChooser(final_multi_snd_ord_chooser.snd_ord_id_pairs),
   ...         sim.OffspringGenerator(ops=[
   ...             sim.IdTagger(),
   ...             sim.PedigreeTagger(),
   ...             sim.Recombinator(rates=0.01)
   ...         ],
   ...             numOffspring=1),
   ...         subPopSize=[2000],
   ...     ),
   ...     gen=1,
   ... )
   (1,
   1,
   1,
   1,
   ........
   1)

Random Mating
=============

After the MAGIC protocol is complete we are left with 50 replicates of a
maize MAGIC population. We have decided to do three generations of random
mating and then analyze the resulting population.

.. code-block:: python

   >>> rmating = breed.MAGIC(multi_prefounders, recombination_rates)
   >>> rmating.random_mating
   Generation: 3
   Generation: 3
   Generation: 3
   .....
   Generation: 4
   Generation: 4
   Generation: 4
   .....
   Generation: 5
   Generation: 5
   Generation: 5

Sample Analysis
===============

In this particular case we were interested in analyzing three different sample
sizes for each replicate. Each replicate we sampled 100, 250 and 500 individuals
and ran each sample through TASSEL. So 50 replicates x 3 sample sizes equals
150 samples run through TASSEL. Kinship matrices etc are calculated for
each sample. The saegus side computations and output are performed by a single
function which calls many other functions.

.. code-block:: python

   >>> sample_sizes = [100, 250, 500]
   >>> analyze.multiple_sample_analyzer(multi_prefounders, sample_sizes, qtl, alleles, allele_effects, 0.7, segregating_loci)

Unfortunately when I run this through TASSEL I am unable to
carry the actual saegus loci so I have to preserve the saegus absolute indexes of
the segregating loci and vice versa. I am unsure why this error occurs.

.. code-block:: python

   >>> saegus_to_tassel_loci = {}
   >>> tassel_to_saegus_loci = {}
   >>> for idx, locus in enumerate(segregating_loci):
   >>>    saegus_to_tassel_loci[locus] = idx
   >>>    tassel_to_saegus_loci[idx] = locus
   >>> infinite_loci_conversions = shelve.open('infinite_loci_conversions')
   >>> infinite_loci_conversions['saegus_to_tassel'] = saegus_to_tassel_loci
   >>> infinite_loci_conversions['tassel_to_saegus'] = tassel_to_saegus_loci
   >>> infinite_loci_conversions.close()

RunID: Infinite Output
======================

I ran TASSEL mixed linear modeling on replicates 0 through 49. I drew independent
samples from each replicate of sizes: 100, 250 and 500. I am hoping to see
increasing statistical power with the increase in sample size for each replicate.
The next set of steps are performed in R using RStudio. I am using
``R Version: 3.2.0`` and the qvalue package to get the qvalues from the ``p`` column
from TASSEL output. I need to evaluate the results using Python.




Formatting Results and Qvalues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The results from TASSEL are edited to contain the following columns:

   + Chr
   + df
   + F
   + p
   + q
   + difference

An example output would be:

.. _results_dataframe_example:

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Chr</th>
          <th>df</th>
          <th>F</th>
          <th>p</th>
          <th>q</th>
          <th>difference</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>1</td>
          <td>2</td>
          <td>1.32812</td>
          <td>0.26592</td>
          <td>0.980281</td>
          <td>0</td>
        </tr>
        <tr>
          <th>1</th>
          <td>1</td>
          <td>2</td>
          <td>2.67238</td>
          <td>0.07008</td>
          <td>0.980281</td>
          <td>0</td>
        </tr>
        <tr>
          <th>2</th>
          <td>1</td>
          <td>2</td>
          <td>1.46384</td>
          <td>0.23235</td>
          <td>0.980281</td>
          <td>0</td>
        </tr>
        <tr>
          <th>3</th>
          <td>1</td>
          <td>2</td>
          <td>1.31119</td>
          <td>0.27043</td>
          <td>0.980281</td>
          <td>0</td>
        </tr>
        <tr>
          <th>4</th>
          <td>1</td>
          <td>2</td>
          <td>1.10830</td>
          <td>0.33094</td>
          <td>0.980281</td>
          <td>0</td>
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
          <th>865</th>
          <td>10</td>
          <td>2</td>
          <td>0.38228</td>
          <td>0.68250</td>
          <td>0.993854</td>
          <td>0</td>
        </tr>
      </tbody>
    </table>
    <p>866 rows × 6 columns</p>
    </div>



Analyzing TASSEL Results: Power and False Positive Rate
=======================================================

In order to determine if the simulation is working as expected I calculuated
`power` and `false positive rate`. Power is defined as the probability to
detect a locus which truly has an effect. In other words it is :math:`1 - B`
where :math:`B` is the `false negative rate`. For `power` I simply counted
the how many QTL were detected and divided by the number of QTL: 10 in this case.
For the false positive rate I counted the number of loci which were declared significant
but had no effect. Shown below is the criterion in code.

.. code-block:: python

   >>> for tds, panel in zip(testing_data_sets, sample_panels):
   ...    for rep in range(50):
   ...        tds[rep, 0] = len(panel[rep][(panel[rep].ix[:, 'q'] < 0.05)
   ...                                     & (panel[rep].ix[:, 'difference'] > 0.0)]) / 10
   ...        tds[rep, 1] = len(panel[rep][(panel[rep].ix[:, 'q'] < 0.05)
   ...                                     & (panel[rep].ix[:, 'difference'] == 0.0)]) / 855


I analyzed each set of results from at each sample size. The results show an
increase in `power` as sample size increases. Albeit the power is low we have
captured a general trend.

Power and FPR Data Format
~~~~~~~~~~~~~~~~~~~~~~~~~

At present the data is in 3 `hdf` files which contain 50 data frames each.
At present things are extremely unorganized and messy; however, I will clean
everything up very soon. The main thing is that I keep everything documented.


Results of Run Infinite
=======================

I corrected the false positive rate calculation and re-ran the tests. The power
increased with increasing sample size; however, the results were very low. The results
are represented in a table below.

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>size_100_power</th>
          <th>size_100_fpr</th>
          <th>size_250_power</th>
          <th>size_250_fpr</th>
          <th>size_500_power</th>
          <th>size_500_fpr</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.2</td>
          <td>0.000000</td>
          <td>0.3</td>
          <td>0.002339</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.003509</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.001170</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>10</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.2</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>11</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>12</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>13</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>14</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>15</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.001170</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>16</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>17</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>18</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.3</td>
          <td>0.002339</td>
        </tr>
        <tr>
          <th>19</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>20</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.001170</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>21</th>
          <td>0.1</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.3</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>22</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>23</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.2</td>
          <td>0.001170</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>24</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>25</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.3</td>
          <td>0.002339</td>
        </tr>
        <tr>
          <th>26</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.001170</td>
          <td>0.2</td>
          <td>0.002339</td>
        </tr>
        <tr>
          <th>27</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>28</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>29</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>30</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.3</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>31</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.001170</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>32</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>33</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>34</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.002339</td>
        </tr>
        <tr>
          <th>35</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>36</th>
          <td>0.1</td>
          <td>0.0</td>
          <td>0.2</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>37</th>
          <td>0.1</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>38</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>39</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>40</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.001170</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>41</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>42</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.002339</td>
          <td>0.2</td>
          <td>0.002339</td>
        </tr>
        <tr>
          <th>43</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>44</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>45</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>46</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.2</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>47</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>48</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.000000</td>
          <td>0.0</td>
          <td>0.001170</td>
        </tr>
        <tr>
          <th>49</th>
          <td>0.0</td>
          <td>0.0</td>
          <td>0.1</td>
          <td>0.000000</td>
          <td>0.1</td>
          <td>0.000000</td>
        </tr>
      </tbody>
    </table>
    </div>


End of Run: Infinite
====================

The overall results of this run were encouraging but warranted a more in depth
analysis of the results. The next run is "Heaven Denies". I will repeat the same
analysis with larger sample sizes.

