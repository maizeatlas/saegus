===============
:mod:`simulate`
===============

.. module:: simulate





:class: Truncation

    ``Truncation`` is a class which encapsulates all parameters and
    functions to perform recurrent truncation selection on a quantitative
    trait. The trait is assumed to have a simple additive basis so all
    phenotypes are calculated by adding contribution of ``n`` loci plus error.



    .. py:function:: generate_f_one(self, pop, recombination_rates, parental_id_pairs)

        Crosses pairs of founders as they are listed in founder indices
        using breed.PairwiseIDChooser

        .. code-block::

           founder_chooser = breed.PairwiseIDChooser(parental_id_pairs)
           if len(parental_id_pairs) % 2 != 0:
               parental_id_pairs.append(random.choice(parental_id_pairs))
           os_size = len(parental_id_pairs)

           logging.info("Creating the F_one population from selected "
                        "founders.")
           # while pop.popSize() > 1:
           pop.evolve(
               preOps=[
                   sim.PyEval(r'"Generation: %d\n" % gen',
                              ),
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

