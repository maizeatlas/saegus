#! usr/bin/python
__author__ = 'John J Dougherty III'
import simuPOP as sim
import random
import collections as col
from . import breed, operators


class Truncation(object):
    """
    ``Truncation`` is a class which encapsulates all parameters and
    functions to perform recurrent truncation selection on a quantitative
    trait. The trait is assumed to have a simple additive basis so all
    phenotypes are calculated by adding contribution of ``n`` loci plus error.
    """

    def __init__(self, generations_of_selection=1,
                 generations_of_random_mating=1,
                 operating_population_size=2000,
                 proportion_of_individuals_saved=0.05,
                 overshoot_as_proportion=0.50,
                 individuals_per_breeding_subpop=5,
                 heritability=0.7,
                 meta_pop_sample_sizes=100,
                 number_of_replicates=1,
                 prefounders_file_name='',
                 ae_file_name='',
                 ):
        self.generations_of_selection = generations_of_selection
        self.generations_of_random_mating = generations_of_random_mating
        self.operating_population_size = operating_population_size
        self.proportion_of_individuals_saved = proportion_of_individuals_saved
        self.overshoot_as_proportion = overshoot_as_proportion
        self.individuals_per_breeding_subpop = individuals_per_breeding_subpop
        self.heritability = heritability
        self.meta_pop_sample_sizes = meta_pop_sample_sizes
        self.number_of_replicates = number_of_replicates

        self.prefounders_file_name = prefounders_file_name
        self.ae_file_name = ae_file_name

        self.breeding_parameters = col.OrderedDict()
        self.number_of_breeding_individuals = int(
            proportion_of_individuals_saved * operating_population_size)
        self.breeding_parameters['number_of_breeding_individuals'] = \
            self.number_of_breeding_individuals
        self.number_of_breeding_subpops = \
            int(
                self.number_of_breeding_individuals / self.individuals_per_breeding_subpop)
        self.breeding_parameters['number_of_breeding_subpops'] = \
            self.number_of_breeding_subpops
        self.total_number_of_offspring_per_generation = \
            int(operating_population_size * (1 + self.overshoot_as_proportion))
        self.breeding_parameters['total_number_of_offspring_per_generation'] = \
            self.total_number_of_offspring_per_generation
        self.offspring_per_breeding_subpop = \
            int(self.total_number_of_offspring_per_generation
                / self.number_of_breeding_subpops)
        self.breeding_parameters['offspring_per_breeding_subpop'] = \
            self.offspring_per_breeding_subpop
        self.offspring_per_female = \
            int(self.offspring_per_breeding_subpop /
                self.individuals_per_breeding_subpop)
        self.breeding_parameters['offspring_per_female'] = \
            self.offspring_per_female
        self.number_of_nonbreeding_individuals = \
            int(self.operating_population_size -
                self.number_of_breeding_individuals)
        self.breeding_parameters['number_of_nonbreeding_individuals'] = \
            self.number_of_nonbreeding_individuals
        self.number_of_offspring_discarded = int(
            self.overshoot_as_proportion * self.operating_population_size)
        self.breeding_parameters['number_of_offspring_discarded'] = \
            self.number_of_offspring_discarded

    @staticmethod
    def pairwise_merge_protocol(pop: sim.Population):
        """
        Convenience function to merge half-populations into breeding subpopulations. Does not test for an even number
        of subpopulations. Testing for an even number of subpopulations is carried out by a different method.
        :param pop: sim.Population split into multiple subpopulations
        :return: Population with half the number of subpopulations

        """
        k = pop.numSubPop() - 1
        j = k - 1
        while j > -1:
            pop.mergeSubPops([j, k], toSubPop=j)
            k -= 2
            j -= 2

    def pop_halver(self, pop):
        """
        Discards half of each sub-population at random and merges the results.
        """
        num_sub_pops = pop.numSubPop()
        size_of_each_sub_pop = pop.subPopSize(subPop=0)
        number_discarded = int(size_of_each_sub_pop / 2)
        for i in range(num_sub_pops):
            removed_ids = random.sample(pop.indInfo('ind_id', subPop=i),
                                        number_discarded)
            pop.removeIndividuals(IDs=removed_ids)
            # self.pairwise_merge_protocol(pop)

    @staticmethod
    def odd_to_even(pop: sim.Population):
        """
        Tests number of subpopulations in 'pop'. If odd then subpop 0 is paired with subpop -1.
        :param pop: sim.Population with more than one subpopulation
        :return: sim.Population with even number of subpopulations.
        """
        if pop.numSubPop() % 2 != 0 and pop.numSubPop() != 1:
            dummy_pop = pop.clone()
            sub_pop_indexes = list(range(pop.numSubPop()))
            dummy_pop.removeSubPops(sub_pop_indexes[1:])
            pop.addIndFrom(dummy_pop)

    def generate_f_one(self, pop, recombination_rates, parental_id_pairs):
        """
        Crosses pairs of founders as they are listed in founder indices
        using breed.PairwiseIDChooser
        """

        founder_chooser = breed.PairwiseIDChooser(parental_id_pairs)
        if len(parental_id_pairs) % 2 != 0:
            parental_id_pairs.append(random.choice(parental_id_pairs))
        os_size = len(parental_id_pairs)

        print("Creating the F_one population from selected founders.")
        # while pop.popSize() > 1:
        pop.evolve(
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
            ],
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(founder_chooser.by_id_pairs),
                sim.OffspringGenerator(ops=[
                    sim.IdTagger(), sim.ParentsTagger(), sim.PedigreeTagger(),
                    sim.Recombinator(rates=recombination_rates)],
                    numOffspring=1),
                subPopSize=os_size,
            ),
            gen=1,
        )

    def recombinatorial_convergence(self, pop, recombination_rates):
        """
        Implements the MAGIC breeding scheme of breeding single individuals
        in pairs determined by the offspring of the initial population. The
        initial population is given by generate_f_one.
        :param pop:
        :type pop:
        :param recombination_rates:
        :type recombination_rates:
        :return:
        :rtype:
        """
        while pop.popSize() > 1:
            new_parents = list(pop.indInfo('ind_id'))
            new_parent_id_pairs = [(pid, pid + 1) for pid in new_parents[::2]]

            if len(new_parent_id_pairs) % 2 != 0 and len(
                    new_parent_id_pairs) != 1:
                new_parent_id_pairs.append(random.choice(new_parent_id_pairs))
            new_os_size = len(new_parent_id_pairs)

            new_founder_chooser = breed.PairwiseIDChooser(new_parent_id_pairs)

            pop.evolve(
                preOps=[
                    sim.PyEval(r'"Generation: %d\t" % gen'),
                    sim.Stat(popSize=True, numOfMales=True),
                    sim.PyEval(r'"popSize: %d\n" % popSize'),
                ],
                matingScheme=sim.HomoMating(
                    sim.PyParentsChooser(new_founder_chooser.by_id_pairs),
                    sim.OffspringGenerator(ops=[
                        sim.IdTagger(), sim.ParentsTagger(),
                        sim.PedigreeTagger(),
                        sim.Recombinator(rates=recombination_rates)],
                        numOffspring=1),
                    subPopSize=new_os_size,
                ),
                gen=1,
            )

    def expand_by_selfing(self, pop, recombination_rates):
        """
        Specific for plant populations capable of selfing.
        Creates an F2 subpopulations generation by selfing the individuals of
        'pop'. Works on a population with one or more subpopulations.
        """
        # self.odd_to_even(pop)
        num_sub_pops = pop.numSubPop()
        progeny_per_individual = int(self.operating_population_size / 2)
        print("Creating the F_two population.")
        return pop.evolve(
            preOps=[
                sim.MergeSubPops(),
                sim.PyEval(r'"Generation: %d\n" % gen'),
                sim.SplitSubPops(sizes=[1] * num_sub_pops, randomize=False),
            ],
            matingScheme=sim.SelfMating(subPopSize=[
                                                       progeny_per_individual] * num_sub_pops,
                                        numOffspring=progeny_per_individual,
                                        ops=[
                                            sim.Recombinator(
                                                rates=recombination_rates),
                                            sim.IdTagger(),
                                            sim.PedigreeTagger()],
                                        ),
            gen=1,
        )

    def mate_and_merge(self, pop, recombination_rates):
        """
        mate_and_merge was designed to handle scenarios where there is more
        than a single pair of prefounders used from the NAM population and
        where the number of prefounders used is not a power of two.

        :param pop:
        :type pop:
        :param recombination_rates:
        :type recombination_rates:
        :return:
        :rtype:
        """
        starting_gen = pop.vars()['gen']
        print(
            "Initiating recombinatorial convergence at generation: %d" % pop.dvars().gen)
        while pop.popSize() > 1:
            # self.pop_halver(pop)
            # self.odd_to_even(pop)
            # self.pairwise_merge_protocol(pop)
            pop.evolve(
                preOps=[
                    sim.MergeSubPops(),
                    sim.PyEval(r'"Generation: %d\n" % gen'),
                ],
                matingScheme=sim.RandomMating(ops=[
                    sim.IdTagger(), sim.PedigreeTagger(),
                    sim.Recombinator(rates=recombination_rates)]
                ),
                gen=1,
            )

    def interim_random_mating(self, pop, recombination_rates):
        """
        Randomly mates 'pop' for 'gens_of_random_mating' generations to further recombine founder genomes and dissolve
        population structure.
        :param pop: Founder population after mate_and_merge procedure
        :return: Population ready to be subjected to selection
        """
        print("Initiating interim random mating for {} generations.".format(
            self.generations_of_random_mating))
        pop.evolve(
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
            ],
            matingScheme=sim.RandomMating(
                subPopSize=self.operating_population_size,
                ops=[sim.IdTagger(), sim.PedigreeTagger(),
                     sim.Recombinator(
                         rates=recombination_rates)]),
            gen=self.generations_of_random_mating,
        )

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

        sizes = [self.individuals_per_breeding_subpop] \
                * self.number_of_breeding_subpops + \
                [self.number_of_nonbreeding_individuals]
        offspring_pops = [self.offspring_per_breeding_subpop] \
                         * self.number_of_breeding_subpops + [0]

        assert len(sizes) == len(offspring_pops), "Number of parental " \
                                                  "subpopulations must equal " \
                                                  "the number of offspring " \
                                                  "subpopulations"

        sampling_generations = [i for i in range(2,
                                                 self.generations_of_selection,
                                                 2)]

        pc = breed.HalfSibBulkBalanceChooser(
            self.individuals_per_breeding_subpop, self.offspring_per_female)

        pop.evolve(
            initOps=[
                sim.InitInfo(0, infoFields=['generation']),
                operators.GenoAdditive(qtl, aes),
                operators.CalculateErrorVariance(self.heritability),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved),
                operators.MetaPopulation(meta_pop,
                                         self.meta_pop_sample_sizes),
                sim.PyEval(r'"Initial: Sampled %d individuals from generation '
                           r'%d Replicate: %d.\n" % (ss, gen_sampled_from, '
                           r'rep)'),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp']),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp']),
                operators.StoreStatistics(),
                sim.MergeSubPops(),
                operators.Sorter('p'),
            ],
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
                operators.GenoAdditive(qtl, aes, begin=1),
                sim.InfoExec('generation=gen'),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved, begin=1),
                operators.MetaPopulation(meta_pop,
                                         self.meta_pop_sample_sizes,
                                         at=sampling_generations),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp'],
                         at=sampling_generations),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp'],
                         at=sampling_generations),
                operators.StoreStatistics(at=sampling_generations),
                sim.MergeSubPops(),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=sizes, randomize=False),
            ],
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(pc.recursive_pairwise_parent_chooser),
                sim.OffspringGenerator(
                    ops=[sim.IdTagger(), sim.PedigreeTagger(),
                         sim.Recombinator(
                             rates=recombination_rates)],
                    numOffspring=1),
                subPopSize=offspring_pops,
                subPops=list(range(1, self.number_of_breeding_subpops, 1))
            ),
            postOps=[
                sim.MergeSubPops(),
                operators.DiscardRandomOffspring(
                    self.number_of_offspring_discarded),
            ],
            finalOps=[
                sim.InfoExec('generation=gen'),
                operators.GenoAdditive(qtl, aes),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved),
                operators.MetaPopulation(meta_pop, self.meta_pop_sample_sizes),
                sim.PyEval(
                    r'"Final: Sampled %d individuals from generation %d\n" '
                    r'% (ss, gen_sampled_from)'),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                operators.Sorter('p'),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp']),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp']),
                operators.StoreStatistics(),
                sim.MergeSubPops(),
                operators.Sorter('p'),
            ],
            gen=self.generations_of_selection)


    def replicate_selection(self, multi_pop, multi_meta_pop, qtl, aes,
                            recombination_rates):
        """
        Runs recurrent truncation selection on a multi-replicate population.

        :param multi_pop: Simulator object of full-sized population
        :param multi_meta_pop: Simulator object of meta-populations
        :param qtl: Loci whose alleles have effects
        :param aes: Allele effect container
        :param recombination_rates: Probabilities for recombination at each locus
        """


        for pop_rep in multi_pop.populations():
            pop_rep.dvars().gen = 0
        for meta_pop_rep in multi_meta_pop.populations():
            meta_pop_rep.dvars().gen = 0

        sizes = [self.individuals_per_breeding_subpop] \
                * self.number_of_breeding_subpops + \
                [self.number_of_nonbreeding_individuals]
        offspring_pops = [self.offspring_per_breeding_subpop] \
                         * self.number_of_breeding_subpops + [0]

        assert len(sizes) == len(offspring_pops), "Number of parental " \
                                                  "subpopulations must equal " \
                                                  "the number of offspring " \
                                                  "subpopulations"

        sampling_generations = [i for i in range(2,
                                                 self.generations_of_selection,
                                                 2)]

        pc = breed.HalfSibBulkBalanceChooser(
            self.individuals_per_breeding_subpop, self.offspring_per_female)

        multi_pop.evolve(
            initOps=[
                sim.InitInfo(0, infoFields=['generation']),
                operators.GenoAdditive(qtl, aes),
                operators.CalculateErrorVariance(self.heritability),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved),
                operators.ReplicateMetaPopulation(multi_meta_pop,
                                         self.meta_pop_sample_sizes),
                sim.PyEval(r'"Initial: Sampled %d individuals from generation '
                           r'%d Replicate: %d.\n" % (ss, gen_sampled_from, '
                           r'rep)'),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp']),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp']),
                operators.StoreStatistics(),
                sim.MergeSubPops(),
                operators.Sorter('p'),
            ],
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
                operators.GenoAdditive(qtl, aes, begin=1),
                sim.InfoExec('generation=gen'),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved, begin=1),
                operators.ReplicateMetaPopulation(multi_meta_pop,
                                         self.meta_pop_sample_sizes,
                                         at=sampling_generations),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp'],
                         at=sampling_generations),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp'],
                         at=sampling_generations),
                operators.StoreStatistics(at=sampling_generations),
                sim.MergeSubPops(),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=sizes, randomize=False),
            ],
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(pc.recursive_pairwise_parent_chooser),
                sim.OffspringGenerator(
                    ops=[sim.IdTagger(), sim.PedigreeTagger(),
                         sim.Recombinator(
                             rates=recombination_rates)],
                    numOffspring=1),
                subPopSize=offspring_pops,
                subPops=list(range(1, self.number_of_breeding_subpops, 1))
            ),
            postOps=[
                sim.MergeSubPops(),
                operators.DiscardRandomOffspring(
                    self.number_of_offspring_discarded),
            ],
            finalOps=[
                sim.InfoExec('generation=gen'),
                operators.GenoAdditive(qtl, aes),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved),
                operators.ReplicateMetaPopulation(multi_meta_pop,
                                          self.meta_pop_sample_sizes),
                sim.PyEval(
                    r'"Final: Sampled %d individuals from generation %d\n" '
                    r'% (ss, gen_sampled_from)'),
                operators.Sorter('p'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=False),
                operators.Sorter('p'),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp']),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp']),
                operators.StoreStatistics(),
                sim.MergeSubPops(),
                operators.Sorter('p'),
            ],
            gen=self.generations_of_selection)





class Drift(object):
    """
    Class which has functions similar to :class Truncation: but simulates
    conditions of genetic drift instead of recurrent truncation selection.
    """

    """
    ``Truncation`` is a class which encapsulates all parameters and
    functions to perform recurrent truncation selection on a quantitative
    trait. The trait is assumed to have a simple additive basis so all
    phenotypes are calculated by adding contribution of ``n`` loci plus error.
    """

    def __init__(self, generations_of_drift=1,
                 generations_of_random_mating=1,
                 operating_population_size=2000,
                 proportion_of_individuals_saved=0.05,
                 overshoot_as_proportion=0.50,
                 individuals_per_breeding_subpop=5,
                 heritability=0.7,
                 meta_pop_sample_sizes=100,
                 number_of_replicates=1,
                 prefounders_file_name='',
                 ae_file_name='',
                 ):
        self.generations_of_drift = generations_of_drift
        self.generations_of_random_mating = generations_of_random_mating
        self.operating_population_size = operating_population_size
        self.proportion_of_individuals_saved = proportion_of_individuals_saved
        self.overshoot_as_proportion = overshoot_as_proportion
        self.individuals_per_breeding_subpop = individuals_per_breeding_subpop
        self.heritability = heritability
        self.meta_pop_sample_sizes = meta_pop_sample_sizes
        self.number_of_replicates = number_of_replicates

        self.prefounders_file_name = prefounders_file_name
        self.ae_file_name = ae_file_name

        self.breeding_parameters = col.OrderedDict()
        self.number_of_breeding_individuals = int(
            proportion_of_individuals_saved * operating_population_size)
        self.breeding_parameters['number_of_breeding_individuals'] = \
            self.number_of_breeding_individuals
        self.number_of_breeding_subpops = \
            int(
                self.number_of_breeding_individuals / self.individuals_per_breeding_subpop)
        self.breeding_parameters['number_of_breeding_subpops'] = \
            self.number_of_breeding_subpops
        self.total_number_of_offspring_per_generation = \
            int(operating_population_size * (1 + self.overshoot_as_proportion))
        self.breeding_parameters['total_number_of_offspring_per_generation'] = \
            self.total_number_of_offspring_per_generation
        self.offspring_per_breeding_subpop = \
            int(self.total_number_of_offspring_per_generation
                / self.number_of_breeding_subpops)
        self.breeding_parameters['offspring_per_breeding_subpop'] = \
            self.offspring_per_breeding_subpop
        self.offspring_per_female = \
            int(self.offspring_per_breeding_subpop /
                self.individuals_per_breeding_subpop)
        self.breeding_parameters['offspring_per_female'] = \
            self.offspring_per_female
        self.number_of_nonbreeding_individuals = \
            int(self.operating_population_size -
                self.number_of_breeding_individuals)
        self.breeding_parameters['number_of_nonbreeding_individuals'] = \
            self.number_of_nonbreeding_individuals
        self.number_of_offspring_discarded = int(
            self.overshoot_as_proportion * self.operating_population_size)
        self.breeding_parameters['number_of_offspring_discarded'] = \
            self.number_of_offspring_discarded

    def determine_breeding_parameters(self, operating_population_size,
                                      proportion_of_individuals_saved,
                                      overshoot_as_proportion,
                                      individuals_per_breeding_subpop):
        """
        The parameters which determine the breeding population and how
        it is structured into sub-populations are determine by the four
        parameters:
        :param operating_population_size:
        :type operating_population_size:
        :param proportion_of_individuals_saved:
        :type proportion_of_individuals_saved:
        :param overshoot_as_proportion:
        :type overshoot_as_proportion:
        :param individuals_per_breeding_subpop:
        :type individuals_per_breeding_subpop:
        :return:
        :rtype:
        """
        self.breeding_parameters = col.OrderedDict()
        self.number_of_breeding_individuals = \
            int(proportion_of_individuals_saved * operating_population_size)
        self.breeding_parameters['number_of_breeding_individuals'] \
            = self.number_of_breeding_individuals
        self.number_of_breeding_subpops = \
            int(self.number_of_breeding_individuals /
                self.individuals_per_breeding_subpop)
        self.breeding_parameters['number_of_breeding_subpops'] = \
            self.number_of_breeding_subpops
        self.total_number_of_offspring_per_generation = \
            int(operating_population_size * (1 + overshoot_as_proportion))
        self.breeding_parameters['total_number_of_offspring_per_generation'] = \
            self.total_number_of_offspring_per_generation
        self.offspring_per_breeding_subpop = \
            int(self.total_number_of_offspring_per_generation /
                self.number_of_breeding_subpops)
        self.breeding_parameters['offspring_per_breeding_subpop'] = \
            self.offspring_per_breeding_subpop
        self.offspring_per_female = int(self.offspring_per_breeding_subpop /
                                        self.individuals_per_breeding_subpop)
        self.breeding_parameters['offspring_per_female'] = \
            self.offspring_per_female
        self.number_of_nonbreeding_individuals = \
            int(self.operating_population_size -
                self.number_of_breeding_individuals)
        self.breeding_parameters['number_of_nonbreeding_individuals'] = \
            self.number_of_nonbreeding_individuals
        self.number_of_offspring_discarded = \
            int(self.overshoot_as_proportion * self.operating_population_size)
        self.breeding_parameters['number_of_offspring_discarded'] = \
            self.number_of_offspring_discarded

    @staticmethod
    def pairwise_merge_protocol(pop: sim.Population):
        """
        Convenience function to merge half-populations into breeding subpopulations. Does not test for an even number
        of subpopulations. Testing for an even number of subpopulations is carried out by a different method.
        :param pop: sim.Population split into multiple subpopulations
        :return: Population with half the number of subpopulations
        """
        k = pop.numSubPop() - 1
        j = k - 1
        while j > -1:
            pop.mergeSubPops([j, k], toSubPop=j)
            k -= 2
            j -= 2

    def pop_halver(self, pop):
        """
        Discards half of each subpopulation at random and merges the results.
        """
        num_sub_pops = pop.numSubPop()
        size_of_each_sub_pop = pop.subPopSize(subPop=0)
        number_discarded = int(size_of_each_sub_pop / 2)
        for i in range(num_sub_pops):
            removed_ids = random.sample(pop.indInfo('ind_id', subPop=i),
                                        number_discarded)
            pop.removeIndividuals(IDs=removed_ids)
            # self.pairwise_merge_protocol(pop)

    @staticmethod
    def odd_to_even(pop: sim.Population):
        """
        Tests number of subpopulations in 'pop'. If odd then subpop 0 is paired with subpop -1.
        :param pop: sim.Population with more than one subpopulation
        :return: sim.Population with even number of subpopulations.
        """
        if pop.numSubPop() % 2 != 0 and pop.numSubPop() != 1:
            dummy_pop = pop.clone()
            sub_pop_indexes = list(range(pop.numSubPop()))
            dummy_pop.removeSubPops(sub_pop_indexes[1:])
            pop.addIndFrom(dummy_pop)

    def generate_f_one(self, pop, recombination_rates, parental_id_pairs):
        """
        Crosses pairs of founders as they are listed in founder indices
        using breed.PairwiseIDChooser
        """

        founder_chooser = breed.PairwiseIDChooser(parental_id_pairs)
        if len(parental_id_pairs) % 2 != 0:
            parental_id_pairs.append(random.choice(parental_id_pairs))
        offspring_pop_sizes = int(len(parental_id_pairs) / 2)

        print("Creating the F_one population from selected founders.")
        pop.evolve(
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
            ],
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(founder_chooser.by_id_pairs),
                sim.OffspringGenerator(ops=[
                    sim.IdTagger(), sim.ParentsTagger(), sim.PedigreeTagger(),
                    sim.Recombinator(rates=recombination_rates)],
                    numOffspring=1),
                subPopSize=offspring_pop_sizes),
            postOps=sim.SplitSubPops(sizes=[1] * offspring_pop_sizes),
            gen=1,
        )

    def expand_by_selfing(self, pop, recombination_rates):
        """
        Specific for plant populations capable of selfing.
        Creates an F2 subpopulations generation by selfing the individuals of
        'pop'. Works on a population with one or more subpopulations.
        """
        # self.odd_to_even(pop)
        num_sub_pops = pop.numSubPop()
        progeny_per_individual = int(self.operating_population_size / 2)
        print("Creating the F_two population.")
        return pop.evolve(
            preOps=[
                sim.MergeSubPops(),
                sim.PyEval(r'"Generation: %d\n" % gen'),
                sim.SplitSubPops(sizes=[1] * num_sub_pops, randomize=False),
            ],
            matingScheme=sim.SelfMating(subPopSize=[
                                                       progeny_per_individual] * num_sub_pops,
                                        numOffspring=progeny_per_individual,
                                        ops=[
                                            sim.Recombinator(
                                                rates=recombination_rates),
                                            sim.IdTagger(),
                                            sim.PedigreeTagger()],
                                        ),
            gen=1,
        )

    def mate_and_merge(self, pop, recombination_rates):
        """
        mate_and_merge was designed to handle scenarios where there is more
        than a single pair of prefounders used from the NAM population and
        where the number of prefounders used is not a power of two.

        :param pop:
        :type pop:
        :param recombination_rates:
        :type recombination_rates:
        :return:
        :rtype:
        """
        starting_gen = pop.vars()['gen']
        print(
            "Initiating recombinatorial convergence at generation: %d" % pop.dvars().gen)
        while pop.numSubPop() > 1:
            # self.pop_halver(pop)
            self.odd_to_even(pop)
            self.pairwise_merge_protocol(pop)
            sub_pop_sizes = list(pop.subPopSizes())
            pop.evolve(
                preOps=[
                    sim.MergeSubPops(),
                    sim.PyEval(r'"Generation: %d\n" % gen'),
                    sim.SplitSubPops(sizes=sub_pop_sizes, randomize=False),
                ],
                matingScheme=sim.HomoMating(sim.SequentialParentChooser(),
                                            ops=[
                                                sim.Recombinator(
                                                    rates=recombination_rates),
                                                sim.IdTagger(),
                                                sim.PedigreeTagger()],
                                            numOffspring=1,
                                            subPopSize=sub_pop_sizes),
                gen=1,
            )

    def interim_random_mating(self, pop, recombination_rates):
        """
        Randomly mates 'pop' for 'gens_of_random_mating' generations to further recombine founder genomes and dissolve
        population structure.
        :param pop: Founder population after mate_and_merge procedure
        :return: Population ready to be subjected to selection
        """
        print("Initiating interim random mating for {} generations.".format(
            self.generations_of_random_mating))
        pop.evolve(
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
            ],
            matingScheme=sim.RandomMating(
                subPopSize=self.operating_population_size,
                ops=[sim.IdTagger(), sim.PedigreeTagger(),
                     sim.Recombinator(
                         rates=recombination_rates)]),
            gen=self.generations_of_random_mating,
        )

    def recurrent_drift_selection(self, pop, meta_pop, qtl, aes,
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

        sizes = [self.individuals_per_breeding_subpop] \
                * self.number_of_breeding_subpops + \
                [self.number_of_nonbreeding_individuals]
        offspring_pops = [self.offspring_per_breeding_subpop] \
                         * self.number_of_breeding_subpops + [0]

        assert len(sizes) == len(offspring_pops), "Number of parental " \
                                                  "subpopulations must equal " \
                                                  "the number of offspring " \
                                                  "subpopulations"

        sampling_generations = [i for i in range(2,
                                                 self.generations_of_drift, 2)]

        pc = breed.HalfSibBulkBalanceChooser(
            self.individuals_per_breeding_subpop, self.offspring_per_female)

        pop.evolve(
            initOps=[
                sim.InitInfo(0, infoFields=['generation']),
                operators.GenoAdditive(qtl, aes),
                operators.CalculateErrorVariance(self.heritability),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved),
                operators.MetaPopulation(meta_pop,
                                         self.meta_pop_sample_sizes),
                sim.PyEval(r'"Initial: Sampled %d individuals from generation '
                           r'%d Replicate: %d.\n" % (ss, gen_sampled_from, '
                           r'rep)'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=True),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp']),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp']),
                operators.StoreStatistics(),
                sim.MergeSubPops(),
            ],
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
                operators.GenoAdditive(qtl, aes, begin=1),
                sim.InfoExec('generation=gen'),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved, begin=1),
                operators.MetaPopulation(meta_pop,
                                         self.meta_pop_sample_sizes,
                                         at=sampling_generations),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=True),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp'],
                         at=sampling_generations),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp'],
                         at=sampling_generations),
                operators.StoreStatistics(at=sampling_generations),
                sim.MergeSubPops(),
                sim.SplitSubPops(sizes=sizes, randomize=True),
            ],
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(pc.recursive_pairwise_parent_chooser),
                sim.OffspringGenerator(
                    ops=[sim.IdTagger(), sim.PedigreeTagger(),
                         sim.Recombinator(
                             rates=recombination_rates)],
                    numOffspring=1),
                subPopSize=offspring_pops,
                subPops=list(range(1, self.number_of_breeding_subpops, 1))
            ),
            postOps=[
                sim.MergeSubPops(),
                operators.DiscardRandomOffspring(
                    self.number_of_offspring_discarded),
            ],
            finalOps=[
                sim.InfoExec('generation=gen'),
                operators.GenoAdditive(qtl, aes),
                operators.PhenotypeCalculator(
                    self.proportion_of_individuals_saved),
                operators.MetaPopulation(meta_pop, self.meta_pop_sample_sizes),
                sim.PyEval(
                    r'"Final: Sampled %d individuals from generation %d\n" '
                    r'% (ss, gen_sampled_from)'),
                sim.SplitSubPops(sizes=[self.number_of_breeding_individuals,
                                        self.number_of_nonbreeding_individuals],
                                 randomize=True),
                sim.Stat(meanOfInfo=['g', 'p'], vars=['meanOfInfo',
                                                      'meanOfInfo_sp']),
                sim.Stat(varOfInfo=['g', 'p'], vars=['varOfInfo',
                                                     'varOfInfo_sp']),
                operators.StoreStatistics(),
                sim.MergeSubPops(),
            ],
            gen=self.generations_of_drift)
