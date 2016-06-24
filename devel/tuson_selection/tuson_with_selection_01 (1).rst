
.. code:: python

    import simuOpt
    simuOpt.setOptions(quiet=True, optimized=True, numThreads=4)
    import simuPOP as sim
    import os, numpy as np, pandas as pd, collections as col
    from saegus import analyze, simulate, parameters, breed
    from scipy import stats
    import random

.. code:: python

    tuson = sim.loadPopulation('tuson.pop')

.. code:: python

    artemis = analyze.Study('artemis')

.. code:: python

    sim.stat(tuson, numOfSegSites=sim.ALL_AVAIL, vars=['numOfSegSites', 'segSites', 'numOfFixedSites', 'fixedSites'])
    parameters.randomly_convert_fixed_sites(tuson, tuson.dvars().fixedSites)
    sim.stat(tuson, numOfSegSites=sim.ALL_AVAIL, vars=['numOfSegSites', 'segSites', 'numOfFixedSites', 'fixedSites'])

.. code:: python

    tuson.dvars().numOfFixedSites




.. parsed-literal::

    0




.. code:: python

    sim.stat(tuson, alleleFreq=sim.ALL_AVAIL)

.. code:: python

    sim.stat(tuson, homoFreq=sim.ALL_AVAIL)

.. code:: python

    sim.stat(tuson, heteroFreq=sim.ALL_AVAIL)


.. code:: python

    alleles  = np.array([list(tuson.dvars().alleleFreq[locus].keys()) for locus in range(tuson.totNumLoci())], dtype=np.int8)

.. code:: python

    alleles




.. parsed-literal::

    array([[1, 2],
           [2, 3],
           [2, 3],
           ..., 
           [1, 2],
           [1, 3],
           [1, 3]], dtype=int8)



Alleles for Fixed Sites chosen at random.
-----------------------------------------

.. code:: python

    np.savetxt('alleles_of_tuson_founders.txt', alleles, fmt='%d', delimiter='\t')

.. code:: python

    af = analyze.allele_data(tuson, alleles, range(tuson.totNumLoci()))

.. code:: python

    def expanded_allele_data(pop, allele_data_structure):
        sim.stat(pop, heteroFreq=sim.ALL_AVAIL)
        hetero_frqs = np.array(list(pop.dvars().heteroFreq.values()))
        hetero_column = pd.DataFrame(hetero_frqs, columns=['heterozygote_frequency'])
        return allele_data_structure.join(hetero_column)

.. code:: python

    eaf = expanded_allele_data(tuson, af)

.. code:: python

    formed_mia = np.array(af['minor_allele'], dtype=np.int8)
    formed_maj = np.array(af['major_allele'], dtype=np.int8)

.. code:: python

    eaf['minor_allele'] = formed_mia
    eaf['major_allele'] = formed_maj

.. code:: python

    eaf.to_csv('expanded_tuson_founder_allele_frqs.txt', sep='\t')

.. code:: python

    tuson.addInfoFields(['generation', 'g', 'p'])

.. code:: python

    tuson.save('working_tuson.pop')

.. code:: python

    tuson.asPedigree()

.. code:: python

    type(tuson)




.. parsed-literal::

    simuPOP.simuPOP_op.Pedigree



.. code:: python

    tuson.save("tuson_pedigree.txt", infoFields=['g', 'p'], loci=sim.ALL_AVAIL)

.. code:: python

    tuson.popSize()




.. parsed-literal::

    105



.. code:: python

    alleles




.. parsed-literal::

    array([[1, 2],
           [2, 3],
           [2, 3],
           ..., 
           [1, 2],
           [1, 3],
           [1, 3]], dtype=int8)



.. code:: python

    tuson.asPopulation()

.. code:: python

    type(tuson)




.. parsed-literal::

    simuPOP.simuPOP_op.Pedigree



.. code:: python

    eaf




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



The Tuson Genetic Map
---------------------

.. code:: python

    def parse_recombination_rates(genetic_map_filename):
        """
        Returns a list of crossover probabilities from a genetic map measured in centimorgans.
        """
        genetic_map = pd.read_csv(genetic_map_filename, sep='\t', index_col=None)
        genetic_map.drop(['locus', 'agpv2', 'namZmPRDA', 'namZmPRDS'], axis=1, inplace=True)
        genetic_map = np.array(genetic_map)
        recombination_rates = col.OrderedDict()
        for i in range(1, len(genetic_map), 1):
            if genetic_map[i-1][0] == genetic_map[i][0]:
                recombination_rates[i] = np.divide(np.abs(genetic_map[i][1] - genetic_map[i-1][1]), 100)
            elif genetic_map[i-1][0] != genetic_map[i][0]:
                recombination_rates[i] = 0.0
        recombination_rates[len(genetic_map)] = 0.0
        return list(recombination_rates.values())
    

.. code:: python

    recom_rates = parse_recombination_rates('raw_genetic_map.txt')

Using the parameters.PopulationStructure class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    popst = parameters.PopulationStructure(tuson, 'population_structure_matrix.xlsx', 0.01, 1.0)

.. code:: python

    struct_mating_probs = popst.generate_population_structure()

.. code:: python

    def format_mating_pmfs(population_structure_dict):
        mating_pmfs = {}
        for ind, probabilities in population_structure_dict.items():
            for i, prob in enumerate(probabilities):
                values = []
                probabilites = []
                for i, prob in enumerate(struct_mating_probs[ind]):
                    values.append(i)
                    probabilites.append(prob)
                pmf_values = (values, probabilites)
                mating_pmfs[ind] = stats.rv_discrete(values=pmf_values)
        return mating_pmfs

.. code:: python

    formed_mating_pmfs = format_mating_pmfs(struct_mating_probs)

.. code:: python

    def assign_primary_subpopulation(pop, struct_mating_probabilities):
        primary_subpop = {}
        for ind_id, inheritance_proportions in struct_mating_probabilities.items():
            primary_subpop[ind_id] = float(np.argmax(inheritance_proportions))
        for ind in pop.individuals():
            ind.primary = primary_subpop[ind.ind_id]

.. code:: python

    assign_primary_subpopulation(tuson, struct_mating_probs)

.. code:: python

    tuson.dvars().mating_pmfs = formed_mating_pmfs

.. code:: python

    pop_struct_expansion = breed.ForcedPopulationStructureParentChooser(10000, formed_mating_pmfs)

.. code:: python

    primary_subpopulation_splitter = sim.InfoSplitter(field='primary',
                                                      values=[0.0, 1.0, 2.0, 3.0,
                                                              4.0, 5.0])
    tuson.setVirtualSplitter(primary_subpopulation_splitter)
    

.. code:: python

    sim.tagID(tuson, reset=False)

.. code:: python

    multi_son = sim.Simulator(tuson, rep=5)



.. code:: python

    multi_son.evolve(
        matingScheme=sim.HomoMating(
            sim.PyParentsChooser(pop_struct_expansion.forced_structure_parent_chooser),
            sim.OffspringGenerator(ops=[sim.IdTagger(), sim.ParentsTagger(), sim.PedigreeTagger(),
                                       sim.Recombinator(recom_rates)], numOffspring=1),
                subPopSize=1000),
        gen=1
    )

.. code:: python

    multi_son.evolve(
        matingScheme=sim.RandomMating(ops=[sim.IdTagger(), sim.ParentsTagger(), sim.PedigreeTagger(),
                                       sim.Recombinator(recom_rates)], numOffspring=1,
                subPopSize=1000),
        gen=1,
    )

.. code:: python

    for pop in multi_son.populations():
        print(pop.popSize())


.. code:: python

    trun = simulate.Truncation(4, 1, 1000, 0.05, 0.50, 5, 0.7, 50, 5)

.. code:: python

    print(trun)

.. code:: python

    fi

.. code:: python

    def recurrent_truncation_selection(self, pop, meta_pop, qtl, aes,
                                       recombination_rates):
        """
        Sets up and runs recurrent selection for a number of generations for a
        single replicate population. Samples individuals at specified
        intervals to make a ``meta_pop``.
    
        :param pop: Population which undergoes selection.
        :param meta_pop: Population into which sampled individuals are
        deposited
        :param qtl: List of loci to which allele effects have been assigned
        :param aes: Dictionary of allele effects
        """
    
        pop.dvars().gen = 0
        meta_pop.dvars().gen = 0
    
        sizes = [individuals_per_breeding_subpop] \
                * number_of_breeding_subpops + \
                [number_of_nonbreeding_individuals]
        offspring_pops = [offspring_per_breeding_subpop] \
                         * number_of_breeding_subpops + [0]
    
        assert len(sizes) == len(offspring_pops), "Number of parental " \
                                                  "subpopulations must equal " \
                                                  "the number of offspring " \
                                                  "subpopulations"
    
        sampling_generations = [i for i in range(2, generations_of_selection, 2)]
    
        pc = breed.HalfSibBulkBalanceChooser(individuals_per_breeding_subpop, offspring_per_female)
    
        pop.evolve(
            initOps=[
                sim.InitInfo(0, infoFields=['generation']),
                operators.GenoAdditive(qtl, aes),
                operators.CalculateErrorVariance(heritability),
                operators.PhenotypeCalculator(proportion_of_individuals_saved),
                operators.MetaPopulation(meta_pop,
                                         self.meta_pop_sample_sizes),
                sim.PyEval(r'"Initial: Sampled %d individuals from generation '
                           r'%d Replicate: %d.\n" % (ss, gen_sampled_from, '
                           r'rep)'),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                sim.MergeSubPops(),
                operators.Sorter('p'),
            ],
            
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
                operators.GenoAdditive(qtl, aes, begin=1),
                sim.InfoExec('generation=gen'),
                operators.PhenotypeCalculator(proportion_of_individuals_saved, begin=1),
                operators.MetaPopulation(meta_pop, meta_pop_sample_sizes, at=sampling_generations),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                sim.MergeSubPops(),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=sizes, randomize=False),
            ],
            
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(pc.recursive_pairwise_parent_chooser),
                sim.OffspringGenerator(
                    ops=[sim.IdTagger(), sim.PedigreeTagger(),
                         sim.Recombinator(rates=recombination_rates)],
                    numOffspring=1),
                subPopSize=offspring_pops,
                subPops=list(range(1, self.number_of_breeding_subpops, 1))
            ),
            postOps=[
                sim.MergeSubPops(),
                operators.DiscardRandomOffspring(number_of_offspring_discarded),
            ],
            finalOps=[
                sim.InfoExec('generation=gen'),
                operators.GenoAdditive(qtl, aes),
                operators.PhenotypeCalculator(proportion_of_individuals_saved),
                operators.MetaPopulation(meta_pop, self.meta_pop_sample_sizes),
                sim.PyEval(
                    r'"Final: Sampled %d individuals from generation %d\n" '
                    r'% (ss, gen_sampled_from)'),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                operators.Sorter('p'),
                sim.MergeSubPops(),
                operators.Sorter('p'),
            ],
            gen=self.generations_of_selection)
    

.. code:: python

            sampling_generations = [i for i in range(2,
                                                     self.generations_of_selection,
                                                     2)]

.. code:: python

    class Trun(object):
        
        def __init__(self, generations_selection, start_gen, step):
            self.generations_selection = generations_selection
            self.sampling_generations = [i for i in range(start_gen, generations_selection, step)]

.. code:: python

    tr = Trun(10, 2, 2)

.. code:: python

    tr.sampling_generations

.. code:: python

    tuson
