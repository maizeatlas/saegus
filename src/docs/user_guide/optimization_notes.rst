.. _general_notes_about_optimization:

################################
General Notes About Optimization
################################

A collection of optimization examples of :py:mod:`saegus` functions.


.. _optimizing_qtl_and_allele_effects:

Optimizing Allele Effect Assignment and Genotypic Effect Calculation
====================================================================

The current operators which take an individuals's genotype and look up
the corresponding allele effects then add them to get ``g`` is
:py:class:`operators.GenoAdditive`. The input is a dictionary of allele
effects and the function used to perform the operation uses list comprehension
to go through each genotype and take out the corresponding QTL. For a small
number of QTL this probably doesn't hurt anything; however, when the number
of QTL is large the time spent in simulation rises substantially.


.. _list_comprehension_to_fancy_indexing:

Changing List Comprehension to Numpy Fancy Indexing
===================================================

Numpy provides a method to lookup many values simultaneously in an array. This
method is called **fancy indexing**.

.. code-block:: py
   :caption: Example of fancy indexing

   >>> y = np.arange(35).reshape(5,7)
   >>> y[np.array([0, 2, 4], [0, 1, 2])]
   array([0, 15, 30])

I can take advantage of that to gain a significant speed-up in simulation time.

Old Way
^^^^^^^

.. code-block:: py

   >>> for ind in pop.individuals():
   ...      ind.g = sum((self.allele_effects[locus][ind.genotype(ploidy=0)[locus]]
   ...              for locus in self.qtl)) +\
   ...         sum((self.allele_effects[locus][ind.genotype(ploidy=1)[locus]]
   ...              for locus in self.qtl))

Old Way cProfile Results
^^^^^^^^^^^^^^^^^^^^^^^^

Initial: Sampled 500 individuals from generation 0 Replicate: 0.
Initial: Sampled 500 individuals from generation 0 Replicate: 1.
...
Final: Sampled 500 individuals from generation 10
Final: Sampled 500 individuals from generation 10

         1103271215 function calls in **1131.222** seconds

   Ordered by: standard name

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        5    0.000    0.000    0.000    0.000 <embed>:1(<module>)
        1    0.000    0.000 1131.261 1131.261 <string>:1(<module>)
   110190    0.028    0.000    0.028    0.000 __init__.py:502(__init__)
   110190    0.095    0.000    0.217    0.000 __init__.py:512(dvars)
   440000    0.345    0.000    1.587    0.000 __init__.py:534(ind_setInfo3)
   110000    0.050    0.000    0.468    0.000 __init__.py:540(ind_getInfo3)
        1    0.000    0.000    0.000    0.000 __init__.py:754(__init__)
        5    0.000    0.000    0.000    0.000 _methods.py:43(_count_reduce_items)
        5    0.000    0.000    0.001    0.000 _methods.py:76(_var)
       54    0.000    0.000    0.004    0.000 common.py:75(wake)
        5    0.000    0.000    0.001    0.000 fromnumeric.py:2988(var)
       55    0.001    0.000    0.006    0.000 ioloop.py:928(add_callback)
     1275    0.001    0.000    0.001    0.000 iostream.py:227(_is_master_process)
     1275    0.002    0.000    0.008    0.000 iostream.py:240(_schedule_flush)
     1275    0.003    0.000    0.013    0.000 iostream.py:308(write)
       10    0.000    0.000    0.000    0.000 numeric.py:476(asanyarray)
        1    0.000    0.000    0.000    0.000 operators.py:11(__init__)
        3    0.000    0.000    0.000    0.000 operators.py:148(__init__)
       30    0.001    0.000    0.091    0.003 operators.py:153(add_to_meta_pop)
        5    0.000    0.000    0.002    0.000 operators.py:16(calculate_error_variance)
        3    0.000    0.000    0.000    0.000 operators.py:29(__init__)
       55    0.923    0.017 1124.262   20.441 operators.py:36(additive_model)
110110000  115.413    0.000  545.162    0.000 operators.py:42(<genexpr>)
110110000  115.012    0.000  543.262    0.000 operators.py:44(<genexpr>)
        3    0.000    0.000    0.000    0.000 operators.py:78(__init__)
       55    0.466    0.008    1.988    0.036 operators.py:83(phenotypic_effect_calculator)
       30    0.037    0.001    0.045    0.001 random.py:258(shuffle)
   110000    0.212    0.000    0.267    0.000 random.py:370(normalvariate)
       30    0.000    0.000    0.000    0.000 sampling.py:140(__init__)
       30    0.000    0.000    0.000    0.000 sampling.py:150(prepareSample)
       30    0.000    0.000    0.000    0.000 sampling.py:183(__init__)
       30    0.002    0.000    0.088    0.003 sampling.py:189(drawSample)
       30    0.001    0.000    0.089    0.003 sampling.py:218(drawRandomSample)
       30    0.000    0.000    0.045    0.002 sampling.py:95(random_shuffle)
       60    0.000    0.000    0.000    0.000 sampling.py:98(isSequence)
       30    0.000    0.000    0.000    0.000 simuPOP_std.py:1287(getRNG)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:2766(__init__)
        4    0.000    0.000    0.000    0.000 simuPOP_std.py:4923(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:529(__init__)
       10    0.000    0.000    0.000    0.000 simuPOP_std.py:5297(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:5359(__init__)
220000030  109.268    0.000  361.001    0.000 simuPOP_std.py:553(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:5712(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:5954(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:644(__init__)
        3    0.000    0.000    0.000    0.000 simuPOP_std.py:6487(__init__)
        5    0.000    0.000    0.002    0.000 simuPOP_std.py:6648(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:668(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:7401(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:7574(__init__)
   330005    0.237    0.000    0.891    0.000 simuPOP_std.py:765(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:7824(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:835(__init__)
        3    0.000    0.000    0.000    0.000 simuPOP_std.py:859(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:9196(__init__)
        1    0.001    0.001 1131.261 1131.261 simulate.py:287(replicate_random_mating)
        1    0.000    0.000    0.000    0.000 simulate.py:313(<listcomp>)
       55    0.001    0.000    0.001    0.000 stack_context.py:253(wrap)
        1    0.000    0.000    0.000    0.000 {built-in method HomoMating_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method IdTagger_swiginit}
220000000  496.997    0.000  857.998    0.000 {built-in method Individual_genotype}
   110000    0.200    0.000    0.418    0.000 {built-in method Individual_info}
   220000    0.568    0.000    1.241    0.000 {built-in method Individual_setInfo}
        5    0.000    0.000    0.000    0.000 {built-in method InfoExec_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method InitInfo_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method OffspringGenerator_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method PedigreeTagger_swiginit}
       30    0.040    0.001    0.041    0.001 {built-in method Population_extractIndividuals}
        5    0.001    0.000    0.001    0.000 {built-in method Population_indInfo}
      110    0.001    0.000    0.001    0.000 {built-in method Population_individuals}
       60    0.000    0.000    0.000    0.000 {built-in method Population_popSize}
   110190    0.094    0.000    0.094    0.000 {built-in method Population_vars}
        3    0.000    0.000    0.000    0.000 {built-in method PyEval_swiginit}
       10    0.000    0.000    0.000    0.000 {built-in method PyOperator_swiginit}
    59970    0.008    0.000    0.008    0.000 {built-in method RNG_randUniform}
        1    0.000    0.000    0.000    0.000 {built-in method RandomParentsChooser_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method Recombinator_swiginit}
        1    4.863    4.863 1131.257 1131.257 {built-in method Simulator_evolve}
        1    0.000    0.000    0.000    0.000 {built-in method Simulator_populations}
       10    0.000    0.000    0.000    0.000 {built-in method array}
        1    0.000    0.000 1131.261 1131.261 {built-in method exec}
        3    0.000    0.000    0.000    0.000 {built-in method floatListFunc_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method floatList_swiginit}
       30    0.000    0.000    0.000    0.000 {built-in method getRNG}
       55    0.000    0.000    0.000    0.000 {built-in method get_ident}
     1275    0.000    0.000    0.000    0.000 {built-in method getpid}
      120    0.000    0.000    0.000    0.000 {built-in method hasattr}
        1    0.000    0.000    0.000    0.000 {built-in method intList_swiginit}
     1290    0.000    0.000    0.000    0.000 {built-in method isinstance}
       10    0.000    0.000    0.000    0.000 {built-in method issubclass}
       32    0.000    0.000    0.000    0.000 {built-in method len}
   150507    0.032    0.000    0.032    0.000 {built-in method log}
        5    0.000    0.000    0.000    0.000 {built-in method max}
        1    0.000    0.000    0.000    0.000 {built-in method new_HomoMating}
        1    0.000    0.000    0.000    0.000 {built-in method new_IdTagger}
        5    0.002    0.000    0.002    0.000 {built-in method new_InfoExec}
        1    0.000    0.000    0.000    0.000 {built-in method new_InitInfo}
        1    0.000    0.000    0.000    0.000 {built-in method new_OffspringGenerator}
        1    0.000    0.000    0.000    0.000 {built-in method new_PedigreeTagger}
        3    0.000    0.000    0.000    0.000 {built-in method new_PyEval}
       10    0.000    0.000    0.000    0.000 {built-in method new_PyOperator}
        1    0.000    0.000    0.000    0.000 {built-in method new_RandomParentsChooser}
        1    0.000    0.000    0.000    0.000 {built-in method new_Recombinator}
        3    0.000    0.000    0.000    0.000 {built-in method new_floatListFunc}
        1    0.000    0.000    0.000    0.000 {built-in method new_floatList}
        1    0.000    0.000    0.000    0.000 {built-in method new_intList}
        4    0.000    0.000    0.000    0.000 {built-in method new_opList}
        1    0.000    0.000    0.000    0.000 {built-in method new_stringList}
        1    0.000    0.000    0.000    0.000 {built-in method new_subPopList}
        1    0.000    0.000    0.000    0.000 {built-in method new_uintListFunc}
220000030   65.924    0.000   65.924    0.000 {built-in method new_uintList}
   330005    0.306    0.000    0.306    0.000 {built-in method new_uintString}
        4    0.000    0.000    0.000    0.000 {built-in method opList_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method stringList_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method subPopList_swiginit}
   220000   33.899    0.000 1122.323    0.005 {built-in method sum}
        1    0.000    0.000    0.000    0.000 {built-in method uintListFunc_swiginit}
220000030  185.810    0.000  185.810    0.000 {built-in method uintList_swiginit}
   330005    0.348    0.000    0.348    0.000 {built-in method uintString_swiginit}
       85    0.000    0.000    0.000    0.000 {method 'append' of 'list' objects}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
   301014    0.022    0.000    0.022    0.000 {method 'random' of '_random.Random' objects}
       10    0.000    0.000    0.000    0.000 {method 'reduce' of 'numpy.ufunc' objects}
       54    0.004    0.000    0.004    0.000 {method 'send' of '_socket.socket' objects}
     1275    0.000    0.000    0.000    0.000 {method 'write' of '_io.StringIO' objects}




New Way
^^^^^^^

.. code-block:: py

   >>> for ind in pop.individuals():
   ...         alpha_genotype, beta_genotype = np.asarray(ind.genotype(ploidy=0)), np.asarray(ind.genotype(ploidy=1))
   ...         ind.g = sum(self.allele_effects[range(1478), alpha_genotype]) + sum(self.allele_effects[range(1478), beta_genotype])


New Way cProfile Results
^^^^^^^^^^^^^^^^^^^^^^^^

Initial: Sampled 500 individuals from generation 0 Replicate: 0.
Initial: Sampled 500 individuals from generation 0 Replicate: 1.
...
Final: Sampled 500 individuals from generation 10
Final: Sampled 500 individuals from generation 10
         4370725 function calls in **140.042** seconds

   Ordered by: standard name

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        5    0.000    0.000    0.000    0.000 <embed>:1(<module>)
        1    0.000    0.000  140.075  140.075 <string>:1(<module>)
   110190    0.024    0.000    0.024    0.000 __init__.py:502(__init__)
   110190    0.085    0.000    0.193    0.000 __init__.py:512(dvars)
   440000    0.246    0.000    1.330    0.000 __init__.py:534(ind_setInfo3)
   110000    0.044    0.000    0.415    0.000 __init__.py:540(ind_getInfo3)
        1    0.000    0.000    0.000    0.000 __init__.py:754(__init__)
        5    0.000    0.000    0.000    0.000 _methods.py:43(_count_reduce_items)
        5    0.000    0.000    0.001    0.000 _methods.py:76(_var)
       55    0.000    0.000    0.004    0.000 common.py:75(wake)
        5    0.000    0.000    0.001    0.000 fromnumeric.py:2988(var)
       55    0.001    0.000    0.005    0.000 ioloop.py:928(add_callback)
     1275    0.001    0.000    0.001    0.000 iostream.py:227(_is_master_process)
     1275    0.002    0.000    0.007    0.000 iostream.py:240(_schedule_flush)
     1275    0.003    0.000    0.011    0.000 iostream.py:308(write)
   220000    0.147    0.000   39.368    0.000 numeric.py:406(asarray)
       10    0.000    0.000    0.000    0.000 numeric.py:476(asanyarray)
        1    0.000    0.000    0.000    0.000 operators.py:11(__init__)
        3    0.000    0.000    0.000    0.000 operators.py:148(__init__)
       30    0.000    0.000    0.089    0.003 operators.py:153(add_to_meta_pop)
        5    0.000    0.000    0.001    0.000 operators.py:16(calculate_error_variance)
        3    0.000    0.000    0.000    0.000 operators.py:49(__init__)
       55   44.759    0.814  133.606    2.429 operators.py:56(additive_model)
        3    0.000    0.000    0.000    0.000 operators.py:78(__init__)
       55    0.409    0.007    1.730    0.031 operators.py:83(phenotypic_effect_calculator)
       30    0.035    0.001    0.043    0.001 random.py:258(shuffle)
   110000    0.185    0.000    0.232    0.000 random.py:370(normalvariate)
       30    0.000    0.000    0.000    0.000 sampling.py:140(__init__)
       30    0.000    0.000    0.000    0.000 sampling.py:150(prepareSample)
       30    0.000    0.000    0.000    0.000 sampling.py:183(__init__)
       30    0.001    0.000    0.087    0.003 sampling.py:189(drawSample)
       30    0.001    0.000    0.088    0.003 sampling.py:218(drawRandomSample)
       30    0.000    0.000    0.043    0.001 sampling.py:95(random_shuffle)
       60    0.000    0.000    0.000    0.000 sampling.py:98(isSequence)
       30    0.000    0.000    0.000    0.000 simuPOP_std.py:1287(getRNG)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:2766(__init__)
        4    0.000    0.000    0.000    0.000 simuPOP_std.py:4923(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:529(__init__)
       10    0.000    0.000    0.000    0.000 simuPOP_std.py:5297(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:5359(__init__)
   220030    0.146    0.000    0.463    0.000 simuPOP_std.py:553(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:5712(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:5954(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:644(__init__)
        3    0.000    0.000    0.000    0.000 simuPOP_std.py:6487(__init__)
        5    0.000    0.000    0.001    0.000 simuPOP_std.py:6648(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:668(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:7401(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:7574(__init__)
   330005    0.210    0.000    0.766    0.000 simuPOP_std.py:765(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:7824(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:835(__init__)
        3    0.000    0.000    0.000    0.000 simuPOP_std.py:859(__init__)
        1    0.000    0.000    0.000    0.000 simuPOP_std.py:9196(__init__)
        1    0.001    0.001  140.075  140.075 simulate.py:287(replicate_random_mating)
        1    0.000    0.000    0.000    0.000 simulate.py:313(<listcomp>)
       55    0.000    0.000    0.001    0.000 stack_context.py:253(wrap)
        1    0.000    0.000    0.000    0.000 {built-in method HomoMating_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method IdTagger_swiginit}
   220000    0.683    0.000    1.146    0.000 {built-in method Individual_genotype}
   110000    0.175    0.000    0.372    0.000 {built-in method Individual_info}
   220000    0.514    0.000    1.084    0.000 {built-in method Individual_setInfo}
        5    0.000    0.000    0.000    0.000 {built-in method InfoExec_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method InitInfo_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method OffspringGenerator_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method PedigreeTagger_swiginit}
       30    0.043    0.001    0.043    0.001 {built-in method Population_extractIndividuals}
        5    0.001    0.000    0.001    0.000 {built-in method Population_indInfo}
      110    0.001    0.000    0.001    0.000 {built-in method Population_individuals}
       60    0.000    0.000    0.000    0.000 {built-in method Population_popSize}
   110190    0.084    0.000    0.084    0.000 {built-in method Population_vars}
        3    0.000    0.000    0.000    0.000 {built-in method PyEval_swiginit}
       10    0.000    0.000    0.000    0.000 {built-in method PyOperator_swiginit}
    59970    0.008    0.000    0.008    0.000 {built-in method RNG_randUniform}
        1    0.000    0.000    0.000    0.000 {built-in method RandomParentsChooser_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method Recombinator_swiginit}
        1    4.602    4.602  140.073  140.073 {built-in method Simulator_evolve}
        1    0.000    0.000    0.000    0.000 {built-in method Simulator_populations}
   220010   39.221    0.000   39.221    0.000 {built-in method array}
        1    0.000    0.000  140.075  140.075 {built-in method exec}
        3    0.000    0.000    0.000    0.000 {built-in method floatListFunc_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method floatList_swiginit}
       30    0.000    0.000    0.000    0.000 {built-in method getRNG}
       55    0.000    0.000    0.000    0.000 {built-in method get_ident}
     1275    0.000    0.000    0.000    0.000 {built-in method getpid}
      120    0.000    0.000    0.000    0.000 {built-in method hasattr}
        1    0.000    0.000    0.000    0.000 {built-in method intList_swiginit}
     1290    0.000    0.000    0.000    0.000 {built-in method isinstance}
       10    0.000    0.000    0.000    0.000 {built-in method issubclass}
       32    0.000    0.000    0.000    0.000 {built-in method len}
   150343    0.027    0.000    0.027    0.000 {built-in method log}
        5    0.000    0.000    0.000    0.000 {built-in method max}
        1    0.000    0.000    0.000    0.000 {built-in method new_HomoMating}
        1    0.000    0.000    0.000    0.000 {built-in method new_IdTagger}
        5    0.001    0.000    0.001    0.000 {built-in method new_InfoExec}
        1    0.000    0.000    0.000    0.000 {built-in method new_InitInfo}
        1    0.000    0.000    0.000    0.000 {built-in method new_OffspringGenerator}
        1    0.000    0.000    0.000    0.000 {built-in method new_PedigreeTagger}
        3    0.000    0.000    0.000    0.000 {built-in method new_PyEval}
       10    0.000    0.000    0.000    0.000 {built-in method new_PyOperator}
        1    0.000    0.000    0.000    0.000 {built-in method new_RandomParentsChooser}
        1    0.000    0.000    0.000    0.000 {built-in method new_Recombinator}
        3    0.000    0.000    0.000    0.000 {built-in method new_floatListFunc}
        1    0.000    0.000    0.000    0.000 {built-in method new_floatList}
        1    0.000    0.000    0.000    0.000 {built-in method new_intList}
        4    0.000    0.000    0.000    0.000 {built-in method new_opList}
        1    0.000    0.000    0.000    0.000 {built-in method new_stringList}
        1    0.000    0.000    0.000    0.000 {built-in method new_subPopList}
        1    0.000    0.000    0.000    0.000 {built-in method new_uintListFunc}
   220030    0.106    0.000    0.106    0.000 {built-in method new_uintList}
   330005    0.222    0.000    0.222    0.000 {built-in method new_uintString}
        4    0.000    0.000    0.000    0.000 {built-in method opList_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method stringList_swiginit}
        1    0.000    0.000    0.000    0.000 {built-in method subPopList_swiginit}
   220000   47.483    0.000   47.483    0.000 {built-in method sum}
        1    0.000    0.000    0.000    0.000 {built-in method uintListFunc_swiginit}
   220030    0.211    0.000    0.211    0.000 {built-in method uintList_swiginit}
   330005    0.334    0.000    0.334    0.000 {built-in method uintString_swiginit}
       85    0.000    0.000    0.000    0.000 {method 'append' of 'list' objects}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
   300686    0.021    0.000    0.021    0.000 {method 'random' of '_random.Random' objects}
       10    0.000    0.000    0.000    0.000 {method 'reduce' of 'numpy.ufunc' objects}
       55    0.003    0.000    0.003    0.000 {method 'send' of '_socket.socket' objects}
     1275    0.000    0.000    0.000    0.000 {method 'write' of '_io.StringIO' objects}


Speed Difference of Using Fancy Indexing for Lookups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Order of magnitude difference for using arrays to get allele effects.
Assumes 1478 loci, 1000 qtl, 5 replicates, operating size 2000, sample size 500,
10 generations and random mating. I obtain similar performance for non-random
mating.


.. _optimizing_sample_analyzer_functions:

Op