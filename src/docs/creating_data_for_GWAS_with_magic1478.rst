
.. _magic1478-input-for-TASSEL:

==========================================
Generating Input for TASSEL with MAGIC1478
==========================================

This document records the exact steps I took to create a set
of data which serves as input to TASSEL for mixed-linear-modeling.
In this case we are simply performing random mating to create the population
which we will analyze. Other times we will use a recurrently selected population for analysis.


.. _allele-effects-magic1478:

Allele Effects for MAGIC1478
============================

Our usual population uses all 7386 loci and assigns appropriate ``recombination_rates`` to
create triplets of non-recombining loci. In an attempt to make MAGIC1478 and MAGIC7386 comparable
I have assigned allele effects in MAGIC1478 as three independent draws for each allele at each locus.
For example:

.. code-block:: python

   qtl = [2, 10, 20]

The alleles at loci ``qtl`` are:

.. code-block:: python

   alleles[2], alleles[10], alleles[20]

   ([3, 1], [1, 3], [3, 0])

So then we assign an effect to each ``allele`` at each ``locus`` as such:

.. code-block:: python

   allele_effects = {}
   for locus in qtl:
      allele_effects[locus] = {}
      for allele in alleles[locus]:
         allele_effects[locus][allele] = sum([random.expovariate(1) for i in range(3)])

.. parsed-literal::

   allele_effects

   {
      2: {1: 3.57, 3: 1.874},
      10: {1: 3.29, 3: 3.44},
      20: {0: 2.05, 3: 3.03},
   }




.. _random-mating-magic1478:

Random Mating MAGIC1478 for GWAS
================================





.. code:: python

    def assign_allele_effects(alleles, qtl, multiplicity, distribution_function,
                             *distribution_function_parameters):
        allele_effects = {}
        for locus in qtl:
            allele_effects[locus] = {}
            for allele in alleles[locus]:
                allele_effects[locus][allele] = sum([distribution_function(*distribution_function_parameters) 
                                          for i in range(multiplicity)])
        return allele_effects

.. code:: python

    def assign_additive_g(pop, qtl, allele_effects):
        """
        Calculates genotypic contribution ``g`` by summing the effect of each
        allele at each QTL triplet.
        """
        for ind in pop.individuals():
            genotypic_contribution = \
                sum([
                        allele_effects[locus][ind.genotype(ploidy=0)[locus]] +\
                        allele_effects[locus][ind.genotype(ploidy=1)[locus]]
                     for locus
                     in qtl])
            ind.g = genotypic_contribution

.. code:: python

    magic1478 = sim.loadPopulation('magic_1478.pop')

.. code:: python

    sim.tagID(magic1478, reset=False)

.. code:: python

    genetic_map = shelve.open('magic_1478_genetic_map')
    history = shelve.open('magic_1478_history')
    simulation = shelve.open('magic_1478_simulation_parameters')
    trait = shelve.open('magic_1478_trait_parameters')

.. code:: python

    locus_1478_names = list(range(1478))
    pos_1478_column = list(range(1478))

.. code:: python

    import importlib as imp
    imp.reload(breed)

.. code:: python

    sim.tagID()

.. code:: python

    breed_magic_1478 = breed.MAGIC(magic1478, simulation['recombination_rates'])
    breed_magic_1478.interim_random_mating(3, 2000)

.. parsed-literal::

    Initiating interim random mating for 3 generations.
    Generation: 3
    Generation: 4
    Generation: 5
    


Determining G and P in the 1478 Population
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    ae = assign_allele_effects(simulation['alleles'], trait['qtl'], 3, random.expovariate, 1)

.. code:: python

    ae




.. parsed-literal::

    {2: {1: 1.6039383268614498, 3: 2.795016834003455},
     10: {1: 3.3259920171422936, 3: 3.1695014054478565},
     20: {0: 2.4204478909872953, 3: 4.269861858273051}}



.. code:: python

    assign_additive_g(magic1478, qtl, ae)

.. code:: python

    def calculate_error_variance(pop, heritability):
        """
        Calculates the parameter ``epsilon`` to be used as the variance
        of the error distribution. The error distribution generates noise
        found in real experiments.
        """
        variance_of_g = np.var(pop.indInfo('g'))
        epsilon = variance_of_g*(1/heritability - 1)
        pop.dvars().epsilon = epsilon
    
    def phenotypic_effect_calculator(pop):
        """
        Simulate measurement error by adding random error to genotypic
        contribution.
        """
        for ind in pop.individuals():
            ind.p = ind.g + random.normalvariate(0, pop.dvars().epsilon)

.. code:: python

    heritability = 0.7
    calculate_error_variance(magic1478, heritability)
    print(magic1478.dvars().epsilon)
    phenotypic_effect_calculator(magic1478)


.. parsed-literal::

    0.849336297482
    

.. code:: python

    trait['heritability'] = heritability
    trait['epsilon'] = magic1478.dvars().epsilon
    trait['qtl'] = qtl
    trait['allele_effects'] = ae
    trait['g'] = list(magic1478.indInfo('g'))
    trait['p'] = list(magic1478.indInfo('p'))

.. code:: python

    trait.close()


.. code:: python

    print(np.var(pop.indInfo('p')), np.mean(pop.indInfo('p')))

.. code:: python

    np.var(pop.indInfo('p'))

.. code:: python

    import collections as col

.. code:: python

    Design = col.namedtuple("Design", ["genetic", "history", "simulation", "trait"])

.. code:: python

    trait[]

.. code:: python

    trait['allele_effects'] = allele_effects
    trait['qtl'] = qtl

.. code:: python

    trait['heritability'] = heritability
    trait['epsilon'] = pop.dvars().epsilon
    trait['g'] = list(pop.indInfo('g'))
    trait['p'] = list(pop.indInfo('p'))

.. code:: python

    trait.close()

.. code:: python

    trait = shelve.open('magic_1478_trait_parameters')

.. code:: python

    history.close()

.. code:: python

    def calc_error_variance(pop, heritability, *args, **kwargs):
        operators.CalculateErrorVariance(heritability, *args, **kwargs).apply(pop)

.. code:: python

    for magic_rep in multi_std_pop.populations():
        calc_error_variance(magic_rep, 0.7)


.. code:: python

    multi_g = {0: list(multi_std_pop.population(0).indInfo('g')), 
               1: list(multi_std_pop.population(1).indInfo('g'))}

.. code:: python

    multi_p = {0: list(multi_std_pop.population(0).indInfo('p')), 
               1: list(multi_std_pop.population(1).indInfo('p'))}

.. code:: python

    multi_g[1] == multi_g[0]

.. code:: python

    for magic_rep in multi_std_pop.populations():
        pheno_calc(magic_rep, 0.05)

.. code:: python

    trait_parameter_storeage = shelve.open("magic_random_trait_parameters")
    trait_parameter_storeage['triplet_qtl'] = triplet_qtl
    trait_parameter_storeage['allele_effects'] = allele_effects
    trait_parameter_storeage['epsilon'] = epsilon_reps
    trait_parameter_storeage['g'] = multi_g
    trait_parameter_storeage['p'] = multi_p
    trait_parameter_storeage.close()

.. code:: python

    import importlib as imp

.. code:: python

    for i in range(2):
        sim.stat(multi_std_pop.population(i), alleleFreq=sim.ALL_AVAIL, vars=[''])

.. code:: python

    imp.reload(analyze)

.. code:: python

    pop.dvars().qtl = qtl
    pop.dvars().allele_effects = allele_effects

.. code:: python

    alleles = simrams['alleles']

.. code:: python

    alleles

.. code:: python

    sim

.. code:: python

    pop

.. code:: python

    simupop.stat(pop, alleleFreq=simupop.ALL_AVAIL)
    

.. code:: python

    frq = analyze.Frq(magic1478)

.. code:: python

    af = frq.allele_frequencies(magic1478, alleles, list(range(1478)))

.. code:: python

    minor_alleles = np.array([af['minor', 'alleles'][locus] for locus in range(1478)])

.. code:: python

    ties = [locus for locus in range(magic1478.totNumLoci())
            if af['minor', 'alleles'][locus] == af['major', 'alleles'][locus]]
    
    for st in ties:
        af['major', 'alleles'][st] = list(magic1478.dvars().alleleFreq[st])[0]
        af['minor', 'alleles'][st] = list(magic1478.dvars().alleleFreq[st])[1]
    major_minor_allele_conflicts = sum(np.equal(list(af['minor', 'alleles'].values()),
                                                list(af['major', 'alleles'].values())))

.. code:: python

    pca = analyze.PCA(magic1478, range(magic1478.totNumLoci()), frq)
    minor_ac = pca.calculate_count_matrix(magic1478, af['minor', 'alleles'],
                                  'minor_allele_count.txt')
    eigendata = pca.svd(magic1478, minor_ac)

.. code:: python

    simulation['snp_to_integer'] = snp_to_integer
    simulation['integer_to_snp'] = integer_to_snp

.. code:: python

    individual_names = {}
    for ind in magic1478.individuals():
        individual_names[ind.ind_id] = str(ind.ind_id)


.. code:: python

    gwas = analyze.GWAS(magic1478, individual_names, locus_1478_names, pos_1478_column)
    hmap = gwas.hapmap_formatter(integer_to_snp, 'simulated_hapmap.txt')
    phenos = gwas.trait_formatter('phenotype_vector.txt')
    kinship_matrix = gwas.calc_kinship_matrix(minor_ac, af, 'kinship_matrix.txt')
    pop_struct_matrix = gwas.population_structure_formatter(eigendata, 'structure_matrix.txt')
    #pd.DataFrame(multipop.population(i).dvars().statistics).to_csv(prefix + 'means_and_vars.txt', sep='\t')
    analyze.generate_tassel_gwas_configs('', 
                                         '', 
                                         '', 
                                         '', 
                                         'sim_mlm_gwas_pipeline.xml')

.. code:: python

    simulation.close()

.. code:: python

    hmap = gwas.hapmap_formatter(integer_to_snp, 'simulated_hapmap.txt')



.. code:: python

    def pre_GWAS_grinder(multi_pop, founder_alleles, info_prefix, triplet_qtl, allele_effects):
        for i, pop_rep in enumerate(multi_pop.populations()):
            pop_rep_id = str(pop_rep.dvars().rep)
            prefix = info_prefix + str(pop_rep_id) + '_'
            qtl = triplet_qtl[i][1:-1:3]
            pop_rep.dvars().qtl = qtl
            pop_rep.dvars().triplet_qtl = triplet_qtl[i]
            pop_rep.dvars().allele_effects = allele_effects[i]
            frq = analyze.Frq(pop_rep, )
            af = frq.allele_frequencies(pop_rep, alleles, range(pop_rep.totNumLoci()))
            
            ties = [locus for locus in range(pop_rep.totNumLoci())
                    if af['minor', 'alleles'][locus] == af['major', 'alleles'][locus]]
            
            for st in ties:
                af['major', 'alleles'][st] = list(pop_rep.dvars().alleleFreq[st])[0]
                af['minor', 'alleles'][st] = list(pop_rep.dvars().alleleFreq[st])[1]
            major_minor_allele_conflicts = sum(np.equal(list(af['minor', 'alleles'].values()),
                                                        list(af['major', 'alleles'].values())))
    
            assert major_minor_allele_conflicts == 0, "There is a tie in at least one locus."
            
            pca = analyze.PCA(pop_rep, range(pop_rep.totNumLoci()), frq)
        
    
            minor_ac = pca.calculate_count_matrix(pop_rep, af['minor', 'alleles'],
                                              prefix + 'minor_allele_count.txt')
    
            eigendata = pca.svd(pop_rep, minor_ac)
    
    
            individual_names = {ind.ind_id: info_prefix + pop_rep_id +'_G' +
                            str(int(ind.generation)) +
                            '_I'+str(int(ind.ind_id))
                            for ind in pop_rep.individuals()}
            
            
    
    
            gwas = analyze.GWAS(pop_rep, individual_names, locus_names, pos_column)
            hmap = gwas.hapmap_formatter(integer_to_snp, prefix + 'simulated_hapmap.txt')
            phenos = gwas.trait_formatter(prefix + 'phenotype_vector.txt')
            kinship_matrix = gwas.calc_kinship_matrix(minor_ac, af, prefix + 'kinship_matrix.txt')
            pop_struct_matrix = gwas.population_structure_formatter(eigendata, prefix + 'structure_matrix.txt')
            #pd.DataFrame(multipop.population(i).dvars().statistics).to_csv(prefix + 'means_and_vars.txt', sep='\t')
            analyze.generate_tassel_gwas_configs('', 
                                                 '', 
                                                 '', 
                                                 prefix, 
                                                 'sim_mlm_gwas_pipeline.xml')


.. code:: python

    rd_sample.indInfo('ind_id')


.. code:: python

    for st in ties:
        af['major', 'alleles'][st] = list(pop_rep.dvars().alleleFreq[st])[0]
        af['minor', 'alleles'][st] = list(pop_rep.dvars().alleleFreq[st])[1]
    major_minor_allele_conflicts = sum(np.equal(list(af['minor', 'alleles'].values()),
                                                list(af['major', 'alleles'].values())))
    
    assert major_minor_allele_conflicts == 0, "There is a tie in at least one locus."
    


.. code:: python

    pre_GWAS_grinder(multi_std_pop, alleles, 'magic_rdm_mating_', triplet_qtl, allele_effects)












.. code:: python

    for i, magic_rep in enumerate(multi_std_pop.populations()):
        
        magic_rep_id = str(magic_rep.dvars().rep)
        prefix = 'magic_RM_L10_H07_R' + str(meta_rep_id) + '_'
        
        qtl = triplet_qtl[i][1:-1:3]
        
        meta_rep.dvars().qtl = qtl
        meta_rep.dvars().triplet_qtl = triplet_qtl[i]
        meta_rep.dvars().allele_effects = allele_effects[i]
        
        
        frq = analyze.Frq(meta_rep, triplet_qtl[i], alleles, allele_effects[i])
        af = frq.allele_frequencies(meta_rep, range(meta_rep.totNumLoci()))
        
        
        #qtalleles = frq.rank_allele_effects(meta_rep, triplet_qtl[i], alleles, allele_effects[i])
        ties = [locus for locus in range(meta_rep.totNumLoci())
                if af['minor', 'alleles'][locus] == af['major', 'alleles'][locus]]
    
        for st in ties:
            af['major', 'alleles'][st] = list(meta_rep.dvars().alleleFreq[st])[0]
            af['minor', 'alleles'][st] = list(meta_rep.dvars().alleleFreq[st])[1]
        major_minor_allele_conflicts = sum(np.equal(list(af['minor', 'alleles'].values()),
                                                    list(af['major', 'alleles'].values())))
    
        assert major_minor_allele_conflicts == 0, "There is a tie in at least one locus."
    
        #af_table = frq.allele_frq_table(meta_rep, meta_rep.numSubPop(), af, 
        #                                                       recombination_rates, genetic_map)
        #qtaf_table = analyze.qt_allele_table(meta_rep, qtalleles, allele_effects[i], 10)
        
        #af_table.to_csv(various_simulation_info_prefix + prefix + 'allele_frequency_table.txt', sep=',', index=0)
        #qtaf_table.to_csv(various_simulation_info_prefix + prefix + 'qt_allele_info.txt', sep=',', index=0)
        
        #del af_table, qtaf_table
    
    
    
        pca = analyze.PCA(meta_rep, range(meta_rep.totNumLoci()), frq)
        
    
        minor_ac = pca.calculate_count_matrix(meta_rep, af['minor', 'alleles'],
                                              prefix + 'minor_allele_count.txt')
    
        eigendata = pca.svd(meta_rep, minor_ac)
    
    
        individual_names = {ind.ind_id: 'rs_L10_H07_R'+ meta_rep_id +'_G' +
                            str(int(ind.generation)) +
                            '_I'+str(int(ind.ind_id))
                            for ind in meta_rep.individuals()}
    
        ind_names_for_gwas[meta_rep_id] = individual_names
    
        #meta_rep.save(populations_prefix + prefix + 'metapopulation.pop')
        
    #    names_filename = prefix + 'individual_names.yaml'
    #    with open(names_filename, 'w') as name_stream:
    #        yaml.dump(individual_names, name_stream)
    
    
        
        
    
        gwas = analyze.GWAS(meta_rep, individual_names, locus_names, pos_column)
        hmap = gwas.hapmap_formatter(integer_to_snp, prefix + 'simulated_hapmap.txt')
        phenos = gwas.trait_formatter(prefix + 'phenotype_vector.txt')
        kinship_matrix = gwas.calc_kinship_matrix(minor_ac, af, prefix + 'kinship_matrix.txt')
        pop_struct_matrix = gwas.population_structure_formatter(eigendata, prefix + 'structure_matrix.txt')
        pd.DataFrame(multipop.population(i).dvars().statistics).to_csv(prefix + 'means_and_vars.txt', sep='\t')
        analyze.generate_tassel_gwas_configs('', 
                                             '', 
                                             '', 
                                             prefix, 
                                             'sim_mlm_gwas_pipeline.xml')








.. code:: python

    minor_af_vector = np.zeros(7386)
    minor_af_vector[:] = [meta_rep.dvars(0).alleleFreq[locus][af['minor', 'alleles', 0][locus]] 
                          for locus in range(meta_rep.totNumLoci())]
    
    minor_alleles = np.zeros((7386), dtype=np.int8)
    major_alleles = np.zeros((7386), dtype=np.int8)
    minor_alleles[:] = [af['minor', 'alleles', 0][locus] 
                          for locus in range(meta_rep.totNumLoci())]
    major_alleles[:] = [af['major', 'alleles', 0][locus]
                           for locus in range(meta_rep.totNumLoci())]
    
    minor_ae = np.zeros(7386)
    major_ae = np.zeros(7386)
    for locus in qtl:
        minor_ae[locus] = allele_effects[0][locus][minor_alleles[locus]]
        major_ae[locus] = allele_effects[0][locus][major_alleles[locus]]
    
    
    
    avg_locus_effects = minor_af_vector*minor_ae + (1-minor_af_vector)*major_ae


