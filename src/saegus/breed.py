import random
import numpy as np
import simuPOP as sim

from . import operators


# todo Create documentation for DirectionalSelection

class DirectionalSelection(object):

    def __init__(self, pop, subPop=[]):
        self.pop = pop
        self.sub_pop = subPop
        self.selected_inds = self.pop.indInfo('ind_id', self.sub_pop)


    def dir_sel_chooser(self, pop):
        while True:
            selected_parent_id = random.choice(self.selected_inds)
            yield self.pop.indByID(selected_parent_id), self.pop.indByID(second_parent_id)

class MAGIC(object):
    """
    MAGIC: Multi-parent Advanced Generation Inter Crosses
    MAGIC is a cross-design which incorporates a large amount of geneticn
    diversity. MAGIC uses a 'funneling' strategy whereby pairs of lines are
    crossed with each other until only a single ``line`` remains.
    """

    def __init__(self, pop, founders, recombination_rates):
        """
        An instance of MAGIC is intended to use in a particular population
        assuming recombination rates stay constant throughout breeding
        procedure.
        """
        self.pop = pop
        self.founders = founders
        self.recombination_rates = recombination_rates
        self._convergence = False


    def generate_f_one(self, parental_id_pairs, offspring_per_pair):
        """
        Crosses pairs of founders as they are listed in founder indices.
        using breed.PairwiseIDChooser

        :note: Data is specified as pairs. Testing for even-number unnecessary.
        """

        founder_chooser = PairwiseIDChooser(parental_id_pairs, offspring_per_pair)
        number_of_pairs = len(parental_id_pairs)
        self.pop.evolve(
            preOps=[
#                sim.PyEval(r'"Generation: %d\n" % gen',),
                ],
            matingScheme=sim.HomoMating(
                sim.PyParentsChooser(founder_chooser.by_id_pairs),
                sim.OffspringGenerator(ops=[
                    sim.IdTagger(),
                    sim.PedigreeTagger(),
                    sim.Recombinator(rates=self.recombination_rates)],
                    numOffspring=1),
                subPopSize=[offspring_per_pair * number_of_pairs],
            ),
            gen=1,
        )

    def recombinatorial_convergence(self, multi_replicate_populations,
                                    number_sub_populations,
                                    offspring_per_pair):
        """
        Breeds individuals from different sub-populations together until a
        single hybrid sub-population is created.
        :note:`number_sub_populations*offspring_per_pair should equal operating_population_size.`
        :note:`For the time being only works with powers of 2.`


        :param sim.Simulator multi_replicate_populations:
        :param int number_sub_populations:
        :param int offspring_per_pair:
        :return:
        """
        print("Start of recombinatorial convergence.")
        while number_sub_populations > 1:
            mrc = MultiRandomCross(multi_replicate_populations,
                                   number_sub_populations, offspring_per_pair)
            mothers, fathers = mrc.determine_random_cross()
            multi_snd_order_chooser = MultiSecondOrderPairIDChooser(mothers, fathers)
            print("Prior to convergence: {}".format(number_sub_populations))
            multi_replicate_populations.evolve(
                matingScheme=sim.HomoMating(
                    sim.PyParentsChooser(multi_snd_order_chooser.snd_ord_id_pairs),
                    sim.OffspringGenerator(ops=[sim.IdTagger(), sim.PedigreeTagger(), sim.Recombinator(rates=self.recombination_rates)],
                                           numOffspring=1),
                                           subPopSize=[int(number_sub_populations*offspring_per_pair)]
                                           ),
                gen=1,
            )
            number_sub_populations = int(number_sub_populations/2)
            offspring_per_pair = int(2*offspring_per_pair)
        self._convergence = True


    def random_mating(self, generations_of_random_mating, pop_size):
        """
        Randomly mates 'pop' for 'gens_of_random_mating' generations to further
        recombine founder genomes and dissolve population structure.
        """
        print("Initiating random mating for {} generations.".format(generations_of_random_mating))
        self.pop.evolve(
            matingScheme=sim.RandomMating(
                subPopSize=pop_size,
                ops=[sim.IdTagger(), sim.PedigreeTagger(),
                     sim.Recombinator(rates=self.recombination_rates)]),
            gen=generations_of_random_mating,
        )


class SelfCross(object):
    """
    Contains functions to implement top crossing or self mating to produce
    families for further analysis. Parental ID can be mapped to offspring via
    father_id information field.

    """

    def __init__(self, recombination_rates):
        self.recombination_rates = recombination_rates

    def create_self_crosses(self, existing_pop, offspring_per_individual):
        new_pop_size = offspring_per_individual * existing_pop.popSize()
        existing_pop.evolve(
            matingScheme=sim.SelfMating(replacement=False,
                                        numOffspring=offspring_per_individual,
                                        subPopSize=new_pop_size,
                                        ops=[
                                            sim.IdTagger(),
                                            sim.PedigreeTagger(),
                                            sim.Recombinator(rates=0.01)], ),
            gen=1,
        )

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


        :parameter pairs_of_parents: Pairs of IDs corresponding to individuals
        :parameter offspring_per_parental_pair: Number of offspring to generate for each pair of parents.
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
        This chooser uses a separate list for each parent i.e.
        ::

        >>> female_parent_ids = [1, 3, 5, 7]
        >>> male_parent_ids = [2, 4, 6, 8]
        >>> snd_order_chooser = breed.SecondOrderPaidIDChooser(female_parent_ids, male_parent_ids)
        >>> example_population.evolve(
        ...     matingScheme=sim.HomoMating(
        ...     sim.PyParentsChooser(snd_order_chooser.snd_ord_id_pairs),
        ...     sim.OffspringGenerator(ops=[
        ...         sim.IdTagger(),
        ...         sim.PedigreeTagger(),
        ...         sim.Recombinator(rates=0.01)],
        ...         numOffspring=1),
        ...     subPopSize=2000),
        ...     gen=1,
        ...     )

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


class RandomCross(object):
    """
    A class intended for use with simulations using a single replicate. This
    class contains functions to pre-determine mating partners. The
    pre-determined pairs are stored in arrays and used as parameters for
    ParentChooser functions.
    """

    def __init__(self, single_replicate_population,
                 number_sub_pops, sub_pop_size):
        """
        :parameter single_replicate_population: simuPOP.Population object
        :parameter int number_sub_pops: Determines number of sub-populations of each replicate
        :parameter int sub_pop_size: Uniform size of sub-populations
        """

        self.single_replicate_population = single_replicate_population
        self.number_sub_pops = number_sub_pops
        self.sub_pop_size = sub_pop_size

    def __str__(self):
        return "Number of Sub-Populations: {nbr_sps}\n" \
               "Sub-Population Size: {sp_size}\n"\
            .format(nbr_sps=self.number_sub_pops, sp_size=self.sub_pop_size)

    def converging_random_cross(self):
        """
        Creates a pair of arrays which contain the IDs of mating pairs. The
        order of the IDs in each array determines the mate of the `mother`
        or `father`.
        """

        self.single_replicate_population.splitSubPop(0,
             sizes=[self.sub_pop_size]*self.number_sub_pops)

        cross_choices = np.zeros((self.number_sub_pops,
                                  2 * self.sub_pop_size))

        for sp in range(self.single_replicate_population.numSubPop()):
            cross_choices[sp] = [random.choice(
                self.single_replicate_population.indInfo('ind_id', sp))
                                     for k in range(2*self.sub_pop_size)]

            mother_idxs = list(range(self.number_sub_pops))[::2]
            father_idxs = list(range(self.number_sub_pops))[1::2]

            mother_ids = np.concatenate([cross_choices[m_idx]
                                      for m_idx in mother_idxs])
            father_ids = np.concatenate([cross_choices[f_idx]
                                      for f_idx in father_idxs])

        self.single_replicate_population.mergeSubPops()

        return mother_ids, father_ids



class MultiRandomCross(object):
    """
    A class intended for use with simulations using multiple replicates which
    requires predictable mating among sub-populations of each replicate.
    """

    def __init__(self, multiple_replicate_population,
                 number_sub_pops, sub_pop_size):
        """
        :parameter multiple_replicate_population: simuPOP.Simulator object
        :parameter int number_sub_pops: Determines number of sub-populations of each replicate
        :parameter int sub_pop_size: Uniform size of sub-populations
        """
        self.multiple_replicate_population = multiple_replicate_population
        self.number_sub_pops = number_sub_pops
        self.sub_pop_size = sub_pop_size

    def __str__(self):
        return "Number of Sub-Populations: {nbr_sps}\n" \
               "Sub-Population Size: {sp_size}\n"\
            .format(nbr_sps=self.number_sub_pops, sp_size=self.sub_pop_size)

    def determine_random_cross(self):
        """
        Creates separate dictionaries for IDs of mothers and fathers respectively.
        Entries are keyed corresponding to ``rep`` of the replicate the
        IDs are taken from.
        """

        multi_mothers = {}
        multi_fathers = {}

        for rep in self.multiple_replicate_population.populations():
            rep.splitSubPop(0, sizes=[self.sub_pop_size] * self.number_sub_pops)

            cross_choices = np.zeros((self.number_sub_pops,
                                      2 * self.sub_pop_size))

            for sp in range(rep.numSubPop()):
                cross_choices[sp] = [random.choice(rep.indInfo('ind_id', sp))
                                     for k in range(2*self.sub_pop_size)]

            mother_idxs = list(range(self.number_sub_pops))[::2]
            father_idxs = list(range(self.number_sub_pops))[1::2]

            mothers = np.concatenate([cross_choices[m_idx]
                                      for m_idx in mother_idxs])
            fathers = np.concatenate([cross_choices[f_idx]
                                      for f_idx in father_idxs])
            multi_mothers[rep.dvars().rep] = mothers
            multi_fathers[rep.dvars().rep] = fathers

            rep.mergeSubPops()

        return multi_mothers, multi_fathers


class MultiSecondOrderPairIDChooser(object):
    """
    MultiSecondOrderPairIDChooser supports controlled mating of multiple
    replicate populations.

    """

    def __init__(self,
                 multi_mother_ids,
                 multi_father_ids,
                 offspring_per_parental_pair=1):
        """
        This chooser uses a separate list for each parent i.e.


        female_parent_ids = [1, 3, 5, 7]
         male_parent_ids = [2, 4, 6, 8]


        :parameter female_parent_ids: Dict by rep of lists of individual IDs (selfing allowed)
        :parameter male_parent_ids: Dict by rep of lists of individual IDS (selfing allowed)
        :parameter offspring_per_parental_pair: Family size per pair of parents.

        """
        self.multi_mother_ids = multi_mother_ids
        self.multi_father_ids = multi_father_ids
        self.offspring_per_parental_pair = offspring_per_parental_pair

    def snd_ord_id_pairs(self, pop):
        rep = pop.dvars().rep
        for female_id, male_id in zip(self.multi_mother_ids[rep],
                                      self.multi_father_ids[rep]):
            female = pop.indByID(female_id)
            male = pop.indByID(male_id)
            female.setSex(2)
            male.setSex(1)
            yield pop.indByID(male_id), pop.indByID(female_id)


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


class ForcedPopulationStructureParentChooser(object):
    """
    For use in expanding the 105 Tuson founders into a population of
    ``expanded_population_size`` individuals while maintaining the empirical
    population structure.
    """
    def __init__(self, expanded_population_size, mating_probabilities):
        self.expanded_population_size = expanded_population_size
        self.mating_probabilities = mating_probabilities

    def __str__(self):
        return "Expanded Population Size: {}".format(self.expanded_population_size)

    def forced_structure_parent_chooser(self, pop):

        for i in range(self.expanded_population_size):
            first_random_id = random.choice(list(pop.indInfo('ind_id')))
            first_parent = pop.indByID(first_random_id)
            compatible_mating_subpopulation = \
                self.mating_probabilities[first_random_id].rvs()

            second_random_id = random.choice(list(
                pop.indInfo('ind_id',
                            subPop=[0, compatible_mating_subpopulation])))
            second_parent = pop.indByID(second_random_id)

            if first_parent.sex() == second_parent.sex():
                if first_parent.sex() == 1:
                    second_parent.setSex(2)
                elif first_parent.sex() == 2:
                    second_parent.setSex(1)
            yield pop.indByID(first_random_id), pop.indByID(second_random_id)

