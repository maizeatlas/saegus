.. _tuson_with_selection:

########################################
Tuson Revisited: Simulation of Selection
########################################

Many months ago we simulated genetic drift in the Tuson population. For comparison
we will now simulate truncation selection which mirrors the experimental protocol.
In other simulations of selection we have randomly assigned QTL and allele effects
were determined by random draws from some distribution (i.e. exponential). In
this case we will be using the experimentally estimated allele effects of every
single locus.

Initial Setup and Data Collection
=================================

We will use G:sub:`0` as founder population. In G:sub:`0` some sites are fixed
which are segregating in later generations. In other words the G:sub:`0` sample
failed to capture all of the segregating loci. We assume that this implies the
sites which are fixed in the sample have low frequency alleles.

Conversion of Fixed Sites to Segregating Sites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We determine which sites are fixed in the sample and then randomly assign a
single individual an allele which is **not** present at the fixed locus. This
converts a fixed site into a segregating site of low frequency.

.. code-block:: py
   :caption: Convert fixed sites to low frequency segregating sites

   >>> sim.stat(tuson, numOfSegSites=sim.ALL_AVAIL, vars=['numOfSegSites', 'segSites', 'numOfFixedSites', 'fixedSites'])
   >>> parameters.randomly_convert_fixed_sites(tuson, tuson.dvars().fixedSites)
   >>> tuson.dvars().numOfFixedSites
   0
   >>> tuson.dvars().numOfSegSites
   44445

Interesting Data to Store
~~~~~~~~~~~~~~~~~~~~~~~~~

Since each conversion of fixed sites to low frequency segregating sites is a
different random outcome I store each realization before carrying the simulation
forward. We store:

+ Alleles at each locus
+ Allele frequencies
+ Heterozygote frequencies

Alleles at Each Locus
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: py
   :caption: Alleles
   >>> alleles
   array([[1, 2],
       [2, 3],
       [2, 3],
       ...,
       [1, 2],
       [1, 3],
       [1, 3]], dtype=int8)

.. code-block:: py
   :caption: Allele and heterozygote frequencies

   >>> eaf

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
          <th>heterozygote_frequency</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>1</td>
          <td>0.319048</td>
          <td>2</td>
          <td>0.680952</td>
          <td>0.371429</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2</td>
          <td>0.219048</td>
          <td>3</td>
          <td>0.780952</td>
          <td>0.266667</td>
        </tr>
        <tr>
          <th>2</th>
          <td>3</td>
          <td>0.061905</td>
          <td>2</td>
          <td>0.938095</td>
          <td>0.104762</td>
        </tr>
        <tr>
          <th>3</th>
          <td>1</td>
          <td>0.061905</td>
          <td>3</td>
          <td>0.938095</td>
          <td>0.104762</td>
        </tr>
        <tr>
          <th>4</th>
          <td>3</td>
          <td>0.309524</td>
          <td>1</td>
          <td>0.690476</td>
          <td>0.619048</td>
        </tr>
        <tr>
          <th>5</th>
          <td>3</td>
          <td>0.052381</td>
          <td>1</td>
          <td>0.947619</td>
          <td>0.085714</td>
        </tr>
        <tr>
          <th>6</th>
          <td>1</td>
          <td>0.204762</td>
          <td>3</td>
          <td>0.795238</td>
          <td>0.314286</td>
        </tr>
        <tr>
          <th>7</th>
          <td>1</td>
          <td>0.128571</td>
          <td>3</td>
          <td>0.871429</td>
          <td>0.200000</td>
        </tr>
        <tr>
          <th>8</th>
          <td>1</td>
          <td>0.133333</td>
          <td>3</td>
          <td>0.866667</td>
          <td>0.209524</td>
        </tr>
        <tr>
          <th>9</th>
          <td>3</td>
          <td>0.180952</td>
          <td>2</td>
          <td>0.819048</td>
          <td>0.266667</td>
        </tr>
        <tr>
          <th>10</th>
          <td>3</td>
          <td>0.461905</td>
          <td>1</td>
          <td>0.538095</td>
          <td>0.923810</td>
        </tr>
        <tr>
          <th>11</th>
          <td>1</td>
          <td>0.461905</td>
          <td>2</td>
          <td>0.538095</td>
          <td>0.923810</td>
        </tr>
        <tr>
          <th>12</th>
          <td>1</td>
          <td>0.090476</td>
          <td>3</td>
          <td>0.909524</td>
          <td>0.161905</td>
        </tr>
        <tr>
          <th>13</th>
          <td>3</td>
          <td>0.114286</td>
          <td>1</td>
          <td>0.885714</td>
          <td>0.114286</td>
        </tr>
        <tr>
          <th>14</th>
          <td>2</td>
          <td>0.004762</td>
          <td>1</td>
          <td>0.995238</td>
          <td>0.009524</td>
        </tr>
        <tr>
          <th>15</th>
          <td>1</td>
          <td>0.004762</td>
          <td>3</td>
          <td>0.995238</td>
          <td>0.009524</td>
        </tr>
        <tr>
          <th>16</th>
          <td>3</td>
          <td>0.004762</td>
          <td>1</td>
          <td>0.995238</td>
          <td>0.009524</td>
        </tr>
        <tr>
          <th>17</th>
          <td>2</td>
          <td>0.038095</td>
          <td>1</td>
          <td>0.961905</td>
          <td>0.076190</td>
        </tr>
        <tr>
          <th>18</th>
          <td>2</td>
          <td>0.128571</td>
          <td>3</td>
          <td>0.871429</td>
          <td>0.238095</td>
        </tr>
        <tr>
          <th>19</th>
          <td>1</td>
          <td>0.004762</td>
          <td>2</td>
          <td>0.995238</td>
          <td>0.009524</td>
        </tr>
        <tr>
          <th>20</th>
          <td>1</td>
          <td>0.295238</td>
          <td>3</td>
          <td>0.704762</td>
          <td>0.323810</td>
        </tr>
        <tr>
          <th>21</th>
          <td>3</td>
          <td>0.423810</td>
          <td>1</td>
          <td>0.576190</td>
          <td>0.428571</td>
        </tr>
        <tr>
          <th>22</th>
          <td>1</td>
          <td>0.214286</td>
          <td>3</td>
          <td>0.785714</td>
          <td>0.257143</td>
        </tr>
        <tr>
          <th>23</th>
          <td>1</td>
          <td>0.042857</td>
          <td>3</td>
          <td>0.957143</td>
          <td>0.085714</td>
        </tr>
        <tr>
          <th>24</th>
          <td>0</td>
          <td>0.004762</td>
          <td>1</td>
          <td>0.995238</td>
          <td>0.009524</td>
        </tr>
        <tr>
          <th>25</th>
          <td>1</td>
          <td>0.309524</td>
          <td>3</td>
          <td>0.690476</td>
          <td>0.142857</td>
        </tr>
        <tr>
          <th>26</th>
          <td>1</td>
          <td>0.223810</td>
          <td>3</td>
          <td>0.776190</td>
          <td>0.447619</td>
        </tr>
        <tr>
          <th>27</th>
          <td>2</td>
          <td>0.004762</td>
          <td>3</td>
          <td>0.995238</td>
          <td>0.009524</td>
        </tr>
        <tr>
          <th>28</th>
          <td>1</td>
          <td>0.104762</td>
          <td>3</td>
          <td>0.895238</td>
          <td>0.171429</td>
        </tr>
        <tr>
          <th>29</th>
          <td>1</td>
          <td>0.004762</td>
          <td>3</td>
          <td>0.995238</td>
          <td>0.009524</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>44415</th>
          <td>1</td>
          <td>0.066667</td>
          <td>3</td>
          <td>0.933333</td>
          <td>0.114286</td>
        </tr>
        <tr>
          <th>44416</th>
          <td>3</td>
          <td>0.319048</td>
          <td>1</td>
          <td>0.680952</td>
          <td>0.295238</td>
        </tr>
        <tr>
          <th>44417</th>
          <td>3</td>
          <td>0.333333</td>
          <td>1</td>
          <td>0.666667</td>
          <td>0.285714</td>
        </tr>
        <tr>
          <th>44418</th>
          <td>3</td>
          <td>0.328571</td>
          <td>1</td>
          <td>0.671429</td>
          <td>0.276190</td>
        </tr>
        <tr>
          <th>44419</th>
          <td>1</td>
          <td>0.147619</td>
          <td>3</td>
          <td>0.852381</td>
          <td>0.200000</td>
        </tr>
        <tr>
          <th>44420</th>
          <td>1</td>
          <td>0.419048</td>
          <td>3</td>
          <td>0.580952</td>
          <td>0.380952</td>
        </tr>
        <tr>
          <th>44421</th>
          <td>1</td>
          <td>0.071429</td>
          <td>3</td>
          <td>0.928571</td>
          <td>0.123810</td>
        </tr>
        <tr>
          <th>44422</th>
          <td>1</td>
          <td>0.419048</td>
          <td>3</td>
          <td>0.580952</td>
          <td>0.419048</td>
        </tr>
        <tr>
          <th>44423</th>
          <td>1</td>
          <td>0.166667</td>
          <td>3</td>
          <td>0.833333</td>
          <td>0.238095</td>
        </tr>
        <tr>
          <th>44424</th>
          <td>3</td>
          <td>0.080952</td>
          <td>1</td>
          <td>0.919048</td>
          <td>0.161905</td>
        </tr>
        <tr>
          <th>44425</th>
          <td>3</td>
          <td>0.295238</td>
          <td>1</td>
          <td>0.704762</td>
          <td>0.304762</td>
        </tr>
        <tr>
          <th>44426</th>
          <td>3</td>
          <td>0.180952</td>
          <td>1</td>
          <td>0.819048</td>
          <td>0.285714</td>
        </tr>
        <tr>
          <th>44427</th>
          <td>1</td>
          <td>0.028571</td>
          <td>2</td>
          <td>0.971429</td>
          <td>0.019048</td>
        </tr>
        <tr>
          <th>44428</th>
          <td>3</td>
          <td>0.171429</td>
          <td>1</td>
          <td>0.828571</td>
          <td>0.209524</td>
        </tr>
        <tr>
          <th>44429</th>
          <td>1</td>
          <td>0.080952</td>
          <td>2</td>
          <td>0.919048</td>
          <td>0.104762</td>
        </tr>
        <tr>
          <th>44430</th>
          <td>3</td>
          <td>0.457143</td>
          <td>1</td>
          <td>0.542857</td>
          <td>0.457143</td>
        </tr>
        <tr>
          <th>44431</th>
          <td>1</td>
          <td>0.385714</td>
          <td>3</td>
          <td>0.614286</td>
          <td>0.771429</td>
        </tr>
        <tr>
          <th>44432</th>
          <td>1</td>
          <td>0.176190</td>
          <td>3</td>
          <td>0.823810</td>
          <td>0.295238</td>
        </tr>
        <tr>
          <th>44433</th>
          <td>1</td>
          <td>0.376190</td>
          <td>3</td>
          <td>0.623810</td>
          <td>0.409524</td>
        </tr>
        <tr>
          <th>44434</th>
          <td>1</td>
          <td>0.457143</td>
          <td>3</td>
          <td>0.542857</td>
          <td>0.457143</td>
        </tr>
        <tr>
          <th>44435</th>
          <td>3</td>
          <td>0.438095</td>
          <td>1</td>
          <td>0.561905</td>
          <td>0.476190</td>
        </tr>
        <tr>
          <th>44436</th>
          <td>1</td>
          <td>0.380952</td>
          <td>3</td>
          <td>0.619048</td>
          <td>0.419048</td>
        </tr>
        <tr>
          <th>44437</th>
          <td>1</td>
          <td>0.385714</td>
          <td>3</td>
          <td>0.614286</td>
          <td>0.771429</td>
        </tr>
        <tr>
          <th>44438</th>
          <td>1</td>
          <td>0.457143</td>
          <td>3</td>
          <td>0.542857</td>
          <td>0.457143</td>
        </tr>
        <tr>
          <th>44439</th>
          <td>3</td>
          <td>0.457143</td>
          <td>1</td>
          <td>0.542857</td>
          <td>0.457143</td>
        </tr>
        <tr>
          <th>44440</th>
          <td>1</td>
          <td>0.457143</td>
          <td>3</td>
          <td>0.542857</td>
          <td>0.457143</td>
        </tr>
        <tr>
          <th>44441</th>
          <td>2</td>
          <td>0.457143</td>
          <td>1</td>
          <td>0.542857</td>
          <td>0.457143</td>
        </tr>
        <tr>
          <th>44442</th>
          <td>2</td>
          <td>0.466667</td>
          <td>1</td>
          <td>0.533333</td>
          <td>0.457143</td>
        </tr>
        <tr>
          <th>44443</th>
          <td>3</td>
          <td>0.261905</td>
          <td>1</td>
          <td>0.738095</td>
          <td>0.352381</td>
        </tr>
        <tr>
          <th>44444</th>
          <td>1</td>
          <td>0.266667</td>
          <td>3</td>
          <td>0.733333</td>
          <td>0.247619</td>
        </tr>
      </tbody>
    </table>
    <p>44445 rows × 5 columns</p>
    </div>

Recombination Rates
===================

We use a simple function to obtain the recombination rates as a list of floats.
Recombination rates are derived from the genetic map we were provided with.

.. code-block:: py
   :caption:

   >>> recom_rates = parameters.parse_recombination_rates('raw_genetic_map'.txt)
   >>> recom_rates
   [0.0020926625899999962,
   2.2615580000007186e-05,
   0.00042822784999999361,
   0.00031254837999999729,
   0.0014689310100000075,
   0.00020776456000000111,
   0.0012046017399999975,
   0.0004001773199999992,
   0.0023329853400000022,
   0.00084844494999999573,
   0.00020627060000000697,
   0.0034117589199999989,
   ...,
   ]

Expanding G\:sub:`0` According to Population Structure
======================================================

Expansion of the Tuson founders is performed according to the same method used
in the previous drift simulation.

.. code-block:: py
   :caption: Functions to perform population expansion

   >>> popst = parameters.PopulationStructure(tuson, 'population_structure_matrix.xlsx', 0.01, 1.0)
   >>> struct_mating_probs = popst.generate_population_structure()
   >>> formatted_mating_pmfs = popst.format_mating_pmfs(struct_mating_probs)
   >>> popst.assign_primary_subpopulation(tuson, struct_mating_probs)
   >>> tuson.dvars().mating_pmfs = formatted_mating_pmfs
   >>> popst_expansion = breed.ForcedPopulationStructureParentChooser(expanded_pop_size, formatted_mating_pmfs)
   >>> primary_subpop_splitter = sim.InfoSplitter(field='primary', values=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
   >>> tuson.setVirtualSplitter(primary_subpop_splitter)
   >>> tuson.numVirtualSubPop()
   6

We will use the ``popst_expansion`` to obtain the actual parent_chooser function
required for the simuPOP mating scheme.

Selection
=========

I have not performed a selection simulation in quite some time. There are
many functions and operators which have changed in :py:mod:`saegus` so I
will have to modify the functions which perform recurrent selection. The first
and most obvious modification is that outdated or pointless functions have to
be removed i.e. :py:function:`operators.StoreStatistics`. The prior
version of selection operators created a simuPOP.Population object to store
the meta-population samples. I decided to change this to a dictionary of lists
of populations. I still retain all the functionality but the objects are in a
native Python data structure. As usual the entire selection process is performed
with a single function.

.. code-block:: py
   :caption: Recurrent selection on several replicates

   >>> simulate.replicate_selection(multi_son, meta_pop_sample_library, qtl, allele_effects)

Now I am waiting on Randy for the allele effects and heritability.









Note on :class:`simuPOP.Pedigree`
=================================

I did not understand the purpose of the :class:`Pedigree` until today. A
:class:`Pedigree` object can save the genotypes of the population to a file in
a nice and neat way. The file can be loaded and converted into a population
object. Provides an easy way for non-simuPOP users to have access to the exact
populations.

