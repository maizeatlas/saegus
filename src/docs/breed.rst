============
:mod:`breed`
============


A module which has functions and classes related to customized mating schemes implemented in simuPOP.
``saegus`` is made unique largely because of the way individuals are mated with each other.



:class:`PairwiseIDChooser`
==========================

        :param pairs_of_founders: List of lists with pairs of integers.

        Contains a :class:`PyParentsChooser` which allows the user to select pairs of individuals
        by their ``ind_id`` to be mated together.

        Tests
        -----

                I ran the following tests::

                assert pop.popSize == offspring_per_pair*num_pairs, "Incorrect number of individuals in population."
                assert pop.numSubPop() == num_pairs, "Incorrect number of subpopulations."

