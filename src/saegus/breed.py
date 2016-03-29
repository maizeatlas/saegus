import random
import numpy as np
import simuPOP as sim

from . import operators

class MAGIC(object):
    """
    MAGIC: Multi-parent Advanced Generation Inter Crosses
    MAGIC is a cross-design which incorporates a large amount of geneticn
    diversity. MAGIC uses a 'funneling' strategy whereby pairs of lines are
    crossed with each other until only a single ``line`` remains.
    """

    def __init__(self, pop, recombination_rates):
        """
        An instance of MAGIC is intended to use in a particular population
        assuming recombination rates stay constant throughout breeding
        procedure.
        """
        self.pop = pop
        self. recombination_rates = recombination_rates

    def generate_f_one(self, parental_id_pairs,
                       offspring_per_pair):
        """
        Crosses pairs of *founders* specified by ``parental_id_pairs``.

         :param parental_id_pairs: Nested lists of founder IDs.
         :param offspring_per_pair: How many offspring per pair.

         :note: If there are an uneven number ``parental_id_pairs appends a
         random choice to end of list.
        """

        founder_chooser = PairwiseIDChooser(parental_id_pairs,
                                                  offspring_per_pair)


        if len(parental_id_pairs) % 2 != 0:
            parental_id_pairs.append(random.choice(parental_id_pairs))

        number_of_pairs = len(parental_id_pairs)
        self.pop.evolve(
            preOps=[],
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(founder_chooser.by_id_pairs),
                sim.OffspringGenerator(ops=[
                    sim.IdTagger(),
                    sim.ParentsTagger(),
                    sim.PedigreeTagger(),
                    sim.Recombinator(rates=self.recombination_rates)
                                        ],
                    numOffspring=1),
                subPopSize=[offspring_per_pair * number_of_pairs],
            ),
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

    def restructure_offspring(self, offspring_per_subpop, number_subpops):
        """
        Rearranges offspring after the F_1 mating scenario. F_1 is merged into
        a single population. This function splits single aggregate population into
        uniformly sized sub-populations to easily choose mating pairs.

        """
        self.pop.splitSubPop(0, offspring_per_subpop * number_subpops)
        subpop_list = list(range(self.pop.numSubPop()))

        breeding_groups = []
        for pair in zip(subpop_list[0::2], subpop_list[1::2]):
            first_maters = random.sample(self.pop.indInfo('ind_id', pair[0]), offspring_per_subpop)
            second_maters = random.sample(self.pop.indInfo('ind_id', pair[1]), offspring_per_subpop)
            breeding_groups.append([first_maters, second_maters])
        breeding_array = np.array(breeding_groups)

        return breeding_array

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

    def __init__(self, pairs_of_parents, offspring_per_parental_pair=1):
        """
        **Example**: *Mate selected pairs of parents*

        ::

           >>> founders = [[0, 1], [10, 11]]
           >>> offspring_per_parental_pair = 1
           >>> parent_chooser = PairwiseIDChooser(founders, offspring_per_parental_pair)
           >>> for pair in founders:
           >>>     print(pair)
           ...     [0, 1]
           ...     [10, 11]


        :param pairs_of_parents:
        :type pairs_of_parents:
        :param offspring_per_parental_pair:
        """
        self.pairs_of_parents = pairs_of_parents
        self.offspring_per_parental_pair = offspring_per_parental_pair

    def by_id_pairs(self, pop):
        for pair in self.pairs_of_parents:
            copied_parental_pairs = [pair]*self.offspring_per_parental_pair
            for copied_pair in copied_parental_pairs:
                female_id, male_id = copied_pair
                female = pop.indByID(float(female_id))
                male = pop.indByID(float(male_id))
                female.setSex(2)
                male.setSex(1)
                yield pop.indByID(float(male_id)), pop.indByID(float(female_id))


class SecondOrderPairIDChooser(object):
    """
    MultiSubGroupChooser is a generalization of the PairwiseIDChooser to scenarios
    where multiple sub-populations need to be merged. It is an elegant way
    to solve the mating scheme issue while preserving predictable
    and easily testable code.
    """

    def __init__(self,
                 female_parent_ids,
                 male_parent_ids,
                 offspring_per_parental_pair=1):
        """
        PairwiseIDChooser which allows separate lists for females and separate
        list for males. Instead of providing pairs of parents i.e
        ::

            founders = [[1, 2], [3, 4], [5, 6], [7, 8]]

        This chooser uses a separate list for each parent i.e.
        ::

            female_parent_ids = [1, 3, 5, 7]
             male_parent_ids = [2, 4, 6, 8]


        :parameter female_parent_ids: List of individual IDs (selfing allowed)
        :parameter male_parent_ids: List of individual IDS (selfing allowed)
        :parameter offspring_per_parental_pair: Family size per pair of parents.

        """
        self.female_parent_ids = female_parent_ids
        self.male_parent_ids = male_parent_ids
        self.offspring_per_parental_pair = offspring_per_parental_pair

    def snd_ord_id_pairs(self, pop):
        for female_id, male_id in zip(self.female_parent_ids, self.male_parent_ids):
            female = pop.indByID(female_id)
            male = pop.indByID(male_id)
            female.setSex(2)
            male.setSex(1)
            yield pop.indByID(male_id), pop.indByID(female_id)



class ListsOfIDsChooser(object):
    """
    The user provides a pair of lists of ``float``s corresponding to the female
    and male ``ind_id``s.
    """

    def __init__(self, female_parent_ids, male_parent_ids):
        self.female_parent_ids = female_parent_ids
        self.male_parent_ids = male_parent_ids

    def by_lists(self, pop):
        for female_id, male_id in zip(self.female_parent_ids,
                                      self.male_parent_ids):
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

