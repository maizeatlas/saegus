import random

import simuPOP as sim

from . import operators


class Cross(object):

    def generate_f_one(self, pop: sim.Population):
        """
        A very basic implementation of the F_1 cross between pairs of individuals. Relies on
        pre-formatting of desired mating pairs into an ordered list.
        [1, 3, 5, 11, 12, 9, 22, 2]
        The mating pattern would be:
        1x3, 5x11, 12x9, 22x2.
        Rearranging the order of the indices would
        change which pairs of individuals mate.
        :param pop:
        :type pop:
        :return:
        :rtype:
        """
        pop.dvars().generations[1] = 'F_1'
        pop.dvars().gen = 1
        pairs_of_founders = int(pop.popSize() / 2)
        self.odd_to_even(pop)
        print("Creating the F_one population from selected founders.")
        return pop.evolve(
            preOps=[
                sim.PyEval(r'"Generation: %d\n" % gen'),
                operators.CalcTripletFrequencies(),
                sim.PyExec('triplet_freq[gen]=tripletFreq'),
                sim.SplitSubPops(sizes=[2] * pairs_of_founders, randomize=False),
            ],
            matingScheme=sim.RandomMating(subPopSize=[1] * pairs_of_founders,
                                          ops=[sim.Recombinator(rates=0.01), sim.IdTagger(), sim.PedigreeTagger()]),
            gen=1,
        )

    def generate_f_two(self, pop: sim.Population) -> sim.Population:
        """
        Creates an F2 subpopulations generation by selfing the individuals
        of 'pop'. Works on a population with one or more subpopulations.
        """
        pop.vars()['generations'][2] = 'F_2'
        self.odd_to_even(pop)
        num_sub_pops = pop.numSubPop()
        progeny_per_individual = int(self.selected_population_size/2)
        print("Creating the F_two population.")
        return pop.evolve(
            preOps=[
                sim.MergeSubPops(),
                sim.PyEval(r'"Generation: %d\n" % gen'),
                operators.CalcTripletFreq(),
                sim.PyExec('triplet_freq[gen]=tripletFreq'),
                sim.SplitSubPops(sizes=[1]*num_sub_pops, randomize=False),
            ],
            matingScheme=sim.SelfMating(subPopSize=[progeny_per_individual] * num_sub_pops,
                                        numOffspring=progeny_per_individual,
                                        ops=[sim.Recombinator(rates=0.01), sim.IdTagger(), sim.PedigreeTagger()],
                                        ),
            gen=1,
        )

    def _mate_and_merge(self, pop: sim.Population):
        print("Initiating recombinatorial convergence at generation: %d" % pop.dvars().gen)
        while pop.numSubPop() > 1:
            #self.odd_to_even(pop)
            #self.pairwise_merge_protocol(pop)
            #sub_pop_sizes = list(pop.subPopSizes())
            pop.evolve(
                preOps=[
                    sim.MergeSubPops(),
                    sim.PyEval(r'"Generation: %d\n" % gen'),
                    sim.PyExec('triplet_freq[gen]=tripletFreq'),
                    sim.SplitSubPops(sizes=sub_pop_sizes, randomize=False),
                ],
                matingScheme=sim.RandomMating(ops=[sim.Recombinator(rates=0.01),
                                                   sim.IdTagger(), sim.PedigreeTagger()]),
                gen=1,
            )


class ForcedPopulationStructureParentChooser(object):
    """
    For use in expanding the 105 Tuson founders into a population of
    ``expanded_population_size`` individuals while maintaining the empirical
    population structure.
    """
    def __init__(self, expanded_population_size, mating_probabilities):
        self.expanded_population_size = expanded_population_size
        self.mating_probabilities = mating_probabilities

    def forced_structure_parent_chooser(self, pop):

        for i in range(self.expanded_population_size):
            first_random_id = random.choice(list(pop.indInfo('ind_id')))
            first_parent = pop.indByID(first_random_id)
            compatible_mating_subpopulation = \
                self.mating_probabilities[first_random_id].rvs()

            second_random_id = random.choice(list(
                pop.indInfo(
                    'ind_id', subPop=[0, compatible_mating_subpopulation])))
            second_parent = pop.indByID(second_random_id)

            if first_parent.sex() == second_parent.sex():
                if first_parent.sex() == 1:
                    second_parent.setSex(2)
                elif first_parent.sex() == 2:
                    second_parent.setSex(1)
            yield pop.indByID(first_random_id), pop.indByID(second_random_id)

class PairwiseIDChooser(object):
    """
    Designed for simulations involving the NAM prefounders. This chooser
    provides a generator function which is used as input for a
    PyParentsChooser object during mating. PairwiseIDChooser allows the user to
    control the mating by providing tuples of ID pairs.

    The parental IDs are provided in a list of tuples. The tuple is a pair
    integers corresponding to the ind_id infoField.

    Observes genetics convention of writing the female as the leftmost
    member of the pair.
    """

    def __init__(self, pairs_of_parents):
        """
        **Example**: *Mate selected pairs of parents*

        ::

           >>> founders = [(0, 1), (10, 11)]

        :param pairs_of_parents:
        :type pairs_of_parents:
        """
        self.pairs_of_parents = pairs_of_parents

    def by_id_pairs(self, pop):
        for pair in self.pairs_of_parents:
            female_id, male_id = pair
            female = pop.indByID(float(female_id))
            male = pop.indByID(float(male_id))
            female.setSex(2)
            male.setSex(1)
            yield pop.indByID(float(male_id)), pop.indByID(float(female_id))



class HalfSibBulkBalanceChooser(object):

    def __init__(self, individuals_per_breeding_subpop, offspring_per_female):
        self.individuals_per_breeding_subpop = individuals_per_breeding_subpop
        self.offspring_per_female = offspring_per_female

    def recursive_pairwise_parent_chooser(self, pop, subPop):
        """
        A parent chooser which is used to effect a half-sib mating scheme in
        conjunction with an additive model of a quantitative trait.
        Individuals within the population are sorted in descending order
        according to infoField ``p``. Then individuals are cut into breeding
        units and mated according to the half-sib design.
        :param pop:
        :type pop:
        :param subPop:
        :type subPop:
        :return:
        :rtype:
        """
        males = []
        females = []
        for i in range(self.individuals_per_breeding_subpop):
            pop.individual(i, subPop-1).setSex(1)
            males.append(pop.individual(i, subPop-1))
            pop.individual(i, subPop).setSex(2)
            females.append(pop.individual(i, subPop))
        for f in females:
            for i in range(self.offspring_per_female):
                yield random.choice(males), f


class MaximumRecombinatorialConvergenceChooser(object):

    @staticmethod
    def convergent_parent_chooser(pop: sim.Population, subPop: int):
        """
        Currently not in a working state. Do not use.
        This method is used as the chooser function of sim.PyParentsChooser for a breeding structure which maximizes
        the number of recombinations between descendants.
        :param pop: sim.Population with multiple subpopulations of prefounder descendants
        :param subPop: sim.PyParentsChooser operates subpopulation by subpopulation
        """
        female_chooser = sim.RandomParentChooser(True, sexChoice=sim.FEMALE_ONLY)
        female_chooser.initialize(pop, subPop)

        male_chooser = sim.RandomParentChooser(True, sexChoice=sim.MALE_ONLY)
        male_chooser.initialize(pop, subPop)

        while True:
            f = female_chooser.chooseParents()[0]
            m = male_chooser.chooseParents()[0]
            yield m, f

