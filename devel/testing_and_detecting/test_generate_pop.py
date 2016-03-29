import pytest
import simuOpt
simuOpt.setOptions(alleleType='short', quiet=True)
import simuPOP as sim
import pandas as pd
from saegus import breed, operators, simulate, analyze, parse, parameters
import shelve
import numpy as np
import random
np.set_printoptions(suppress=True, precision=3)


def test_generate_operating_population():

    genetic_map = pd.read_csv('nam_prefounders_genetic_map.txt', index_col=None,
                             sep='\t')

    pf_map = shelve.open('pf_map')
    misc_gmap = shelve.open('misc_gmap')
    uniparams = shelve.open('uniparams')

    locus_names = uniparams['locus_names']
    pos_column = uniparams['pos_column']
    allele_names = uniparams['allele_names']
    snp_to_integer = uniparams['snp_to_integer']
    integer_to_snp = uniparams['integer_to_snp']

    alleles = misc_gmap['alleles']
    chr_cM_positions = misc_gmap['chr_cM_positions']
    cM_positions = misc_gmap['cM_positions']
    integral_valued_loci = misc_gmap['integral_valued_loci']
    relative_integral_valued_loci = misc_gmap['relative_integral_valued_loci']
    recombination_rates = misc_gmap['recombination_rates']

    nam = sim.loadPopulation(uniparams['prefounder_file_name'])
    sim.tagID(nam, reset=True)
    nam.setSubPopName('maize_nam_prefounders', 0)

    selection_statistics = {
        'aggregate': {},
        'selected': {},
        'non-selected': {}
    }

    ind_names_for_gwas = {i: {} for i in range(uniparams[
        'number_of_replicates'])}
    uniparams['meta_pop_sample_sizes'] = {i: 100 for i in
                                          range(0, uniparams['generations_of_selection'] + 1, 2)
                                          }

    s = simulate.Truncation(uniparams['generations_of_selection'],
                           uniparams['generations_of_random_mating'],
                           uniparams['operating_population_size'],
                            uniparams[
                                'proportion_of_individuals_saved'],
                           uniparams['overshoot_as_proportion'],
                       uniparams['individuals_per_breeding_subpop'],
                           uniparams['heritability'],
                           uniparams['meta_pop_sample_sizes'],
                           uniparams['number_of_replicates'])

    ind_names_for_gwas = {i: {} for i in range(uniparams[
        'number_of_replicates'])}

    founders = uniparams['founders']
    replicated_nam = sim.Simulator(nam, rep=2, stealPops=False)
    pop = replicated_nam.extract(0)

    assert pop.popSize() == 26, "Population is too large."

    s.generate_f_one(pop, recombination_rates, founders, 100)

    assert pop.popSize() == 400, "Population should have size: {} after the F_1 mating " \
                                               "procedure." \
                                               "".format(len(founders) * 100)

    #pop.splitSubPop(0, [100] * 4)
    #subpop_list = list(range(pop.numSubPop()))

    intmd_os_struct = s.restructure_offspring(pop, 100, 4)
    snd_order = breed.SecondOrderPairIDChooser(intmd_os_struct, 1)

    pop.evolve(
        preOps=[sim.MergeSubPops()],
        matingScheme=sim.HomoMating(
            sim.PyParentsChooser(snd_order.snd_ord_id_pairs),
            sim.OffspringGenerator(ops=[
                sim.IdTagger(),
                sim.ParentsTagger(),
                sim.PedigreeTagger(),
                sim.Recombinator(rates=recombination_rates)
            ],
                numOffspring=1),
            subPopSize=[200],
        ),
        gen=1,
    )

    assert pop.popSize() == 1, "Population does not have correct size after second round of mating."

    second_intmd_os_struct = s.restructure_offspring(pop, 100, 2)
    third_order = breed.SecondOrderPairIDChooser(second_intmd_os_struct, 1)


    pop.evolve(
        preOps=[sim.MergeSubPops()],
        matingScheme=sim.HomoMating(
            sim.PyParentsChooser(third_order.snd_ord_id_pairs),
            sim.OffspringGenerator(ops=[
                sim.IdTagger(),
                sim.ParentsTagger(),
                sim.PedigreeTagger(),
                sim.Recombinator(rates=recombination_rates)
            ],
                numOffspring=1),
            subPopSize=[100],
        ),
        gen=1,
    )

    assert pop.popSize() == 100, "Second merge of breeding sub-populations. Offspring population does not have " \
                                 "correct size"




