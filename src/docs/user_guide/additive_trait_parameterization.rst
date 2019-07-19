.. _additive_trait_parameterization:

###############################
Additive Trait Parameterization
###############################

.. code-block:: python
   :caption: Modules we will need for this example

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='short', quiet=True)
   >>> import simuPOP as sim
   >>> import pandas as pd
   >>> import numpy as np
   >>> import random
   >>> from saegus import analyze, operators, parameters
   >>> np.set_printoptions(suppress=True, precision=5)

.. _overview_of_additive_trait_example:

Our goal is to show how to parameterize an additive trait in :mod`saegus` code.
We will make use of :mod:`simuPOP` ``infoFields`` to store the information about
the additive trait for each individual. Also we will add a unique identifier
for each individual which will help us verify that the trait is being calculated
correctly. In this example we will randomly choose quantitative trait loci from
among loci which are segregating in ``example_pop``.

.. _preparing_the_population:

Preparing the Population
########################

We will re-load the population and add information fields.
``example_pop`` needs the information fields ``ind_id``, ``g`` and ``p``.

.. _load_the_population:

Load the Population
===================

We will use the population we created in the last step instead of creating
a new population.

.. code-block:: python
   :caption: Loading our example population from a file

   >>> example_pop = sim.loadPopulation('example_pop.pop')


.. _add_information_fields:

Add Information Fields
======================

At present our population has no special information assigned to its members.
Each individual is only a genotype. By default :mod:`simuPOP` uses ``ind_id``,
``father_id``, ``mother_id``, ``father_idx`` and ``mother_idx`` for its very
useful _Tagger functions. We can save some hassle by using these for
identifying individuals.

.. _Tagger: http://simupop.sourceforge.net/manual_svn/build/refManual_ch3_sec10.html

.. code-block:: python
   :caption: Add the ``infoFields`` ``ind_id``, ``g`` and ``p``

   >>> example_pop.infoFields()
   ()
   >>> example_pop.addInfoFields(['ind_id', 'g', 'p'])
   >>> example_pop.infoFields()
   ('ind_id', 'g', 'p')

By default information fields are set to ``0.0``. We can initialize the
``ind_id`` field using a :mod:`simuPOP` function.

.. code-block:: python
   :caption: Initialize individual identifiers

   >>> print(np.array(example_pop.indInfo('ind_id')))
   [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
   >>> sim.tagID(example_pop)
   >>> print(np.array(example_pop.indInfo('ind_id')))
   [   1.    2.    3.    4.    5.    6.    7.    8.    9.   10.   11.   12.
      13.   14.   15.   16.   17.   18.   19.   20.   21.   22.   23.   24.
      25.   26.   27.   28.   29.   30.   31.   32.   33.   34.   35.   36.
      37.   38.   39.   40.   41.   42.   43.   44.   45.   46.   47.   48.
      49.   50.   51.   52.   53.   54.   55.   56.   57.   58.   59.   60.
      61.   62.   63.   64.   65.   66.   67.   68.   69.   70.   71.   72.
      73.   74.   75.   76.   77.   78.   79.   80.   81.   82.   83.   84.
      85.   86.   87.   88.   89.   90.   91.   92.   93.   94.   95.   96.
      97.   98.   99.  100.  101.  102.  103.  104.  105.]

.. note::
   ::
   In this step we converted the output into a np.array for aesthetics
   
.. calculate allele frequencies:

Calculate Allele Frequencies
==========================

.. code-block:: python
   :caption: Using :mod:`simuPOP` to find compute allele frequencies
   
   >>> sim.stat(example_pop, alleleFreq=sim.ALL_AVAIL)

.. _determine_segregating_loci:

Determine Segregating Loci
==========================

For simplicity we will use loci which have more than one allele i.e.
segregating.

.. code-block:: python
   :caption: Using :mod:`simuPOP` to find segregating loci

   >>> sim.stat(example_pop, numOfSegSites=sim.ALL_AVAIL,
   ...              vars=['numOfSegSites', 'segSites', 'fixedSites'])
   >>> example_pop.dvars().numOfSegSites
   42837
   >>> print(example_pop.dvars().segSites[::1000]) # every 1000th segregating locus
   [0, 1040, 2072, 3098, 4124, 5156, 6199, 7217, 8248, 9282, 10338, 11361,
   12392, 13407, 14468, 15502, 16562, 17599, 18637, 19665, 20700, 21766, 22805,
   23813, 24837, 25882, 26910, 27923, 28955, 30026, 31057, 32103, 33142,
   34173, 35185, 36207, 37223, 38243, 39351, 40419, 41477, 42537, 43578]

There are 42,837 segregating loci in this population. ``saegus`` has a function
to put the alleles into an array and assign the alleles at ``qtl`` an effect as
a draw from a specified distribution.

.. _additive_trait:

Additive Trait
##############

We have all the information we need from the previous steps. We will randomly
choose ``20`` QTL from the segregating loci. Both alleles at each QTL are
assigned an effect as a random draw with an exponential distribution.

.. _choose_QTL:

Choosing QTL and Assign Effects
===============================

For this example we will pick 5 loci to designate as quantitative trait loci.

.. code-block:: python
   :caption: Choosing QTL and assigning allele effects

   >>> segregating_loci = example_pop.dvars().segSites
   >>> qtl = sorted(random.sample(segregating_loci, 5))
   >>> print(qtl)
   [6943, 14327, 16868, 17119, 35312]

Every allele is initially assigned an effect of ``0``. Now alleles only at each QTL 
will be assigned a non-zero effect drawn from the Exponential distribution.

.. code-block:: python
   :caption: Assign allele effects using an exponential distribution

   >>> example_run = analyze.Study('example_pop')
   >>> allele_states = example_run.gather_allele_data(example_pop)
   >>> alleles = np.array([allele_states[:, 1], allele_states[:, 2]]).T
   >>> trait = parameters.Trait()
   >>> ae_table = trait.construct_allele_effects_table(alleles, qtl, random.expovariate, 1)
   >>> print(ae_table[qtl]) # qtl only
   [[ 6943.        1.        0.938     3.        0.315]
    [14327.        1.        0.436     2.        2.439]
    [16868.        1.        1.3       3.        0.99 ]
    [17119.        1.        0.28      3.        0.702]
    [35312.        1.        0.449     3.        0.281]]
   >>> print(ae_table) # all loci
   [[     0.      1.      0.      2.      0.]
    [     1.      2.      0.      3.      0.]
    [     2.      2.      0.      3.      0.]
    ...,
    [ 44442.      1.      0.      2.      0.]
    [ 44443.      1.      0.      3.      0.]
    [ 44444.      1.      0.      3.      0.]]

Alternatively, we could use another distribution, such as the Normal.
This overwrites the previously assigned effects.

.. code-block:: python
   :caption: Assign allele effects using a normal distribution

   >>> ae_table = trait.construct_allele_effects_table(alleles, qtl, random.normalvariate, 0, 1)
   >>> print(ae_table[qtl]) # qtl only
   [[ 6943.        1.        1.927     3.       -0.827]
    [14327.        1.       -0.51      2.       -0.649]
    [16868.        1.       -0.863     3.        4.06 ]
    [17119.        1.       -0.292     3.       -0.763]
    [35312.        1.       -0.388     3.        0.148]]

For speed of computation we construct an array of allele effects where the row
of the array corresponds to the locus and the column corresponds to the integer
representing the allele state.

.. code-block:: python
   :caption: Putting the allele effects in an array for speed of computation

   >>> ae_array = trait.construct_ae_array(ae_table, qtl)
   >>> print(ae_array[qtl])
   [[ 0.     1.927  0.    -0.827  0.     0.   ]
    [ 0.    -0.51  -0.649  0.     0.     0.   ]
    [ 0.    -0.863  0.     4.06   0.     0.   ]
    [ 0.    -0.292  0.    -0.763  0.     0.   ]
    [ 0.    -0.388  0.     0.148  0.     0.   ]]

.. _definition_of_g:

Definition of ``g``
===================

``g`` is the sum of the allele effects of an individual's genotype. There is
no noise or error in ``g`` because we have *a priori* determined the allele
effects.

.. code-block:: python
   :caption: Calculating g values

   >>> operators.calculate_g(example_pop, ae_array)
   >>> print(np.array(example_pop.indInfo('g')))
   [-6.702 -5.63  -6.702 -1.195 -5.63  -6.702 -5.695 -4.227 -6.166 -0.658
    -6.702  0.284 -1.334 -3.412 -6.702 -6.166 -1.242  4.736 -5.695  0.349
    -0.724 -5.695 -6.702 -0.658 -3.006 -6.166 -6.702 -6.231 -0.658 -6.702
    -1.195 -5.63  -5.695 -1.195 -3.412 -3.948 -2.405 -5.159 -1.195 -3.948
    -5.63  -1.195 -0.658 -6.166 -6.702 -0.122 -3.412 -6.231 -5.63  -1.195
    -1.195 -0.658 -1.195 -2.47  -0.658 -3.948 -1.195 -5.695 -6.702 -3.412
    -1.195 -3.412 -6.231 -6.702 -6.702 -3.412 -6.166 -6.702 -1.195 -6.841
    -1.473 -6.166 -0.658 -0.658 -1.242 -6.166 -0.122 -5.63  -5.63  -6.231
    -6.166 -0.658 -6.166 -0.724 -5.695 -3.756 -0.122 -6.702 -0.724 -6.166
     1.446 -1.195 -1.334 -3.948 -1.195 -2.47  -0.658 -2.941 -6.166 -0.658
    -0.798 -6.702 -0.724 -5.76  -5.695]

.. _calculating_error:

Calculation of Error Term
=========================

To simulate the experimental noise a term :math:`\epsilon` is added to each
individual's ``g`` value.
:math:`\epsilon` is a random variable with a normal distribution given by
mean :math:`0` and variance given by:

.. math::

   \sigma^2_\epsilon = \frac{V_g - (h^2 * V_g)}{h^2}

where :math:`V_g` is the variance of ``g`` and :math:`h^2` is the
narrow sense heritability.

.. math::

   \varepsilon \sim \mathcal{N} (0,\sigma^2_\epsilon)

Hence an individual's value of ``p`` is calculated by

.. math::

   p = g + \epsilon

.. _calculating_p:

Calculating ``p``
=================

It is straightforward to calculate ``p`` for the population but we already
have a function to make it even easier for ourselves.

.. code-block:: python
   :caption: Computing ``p`` for each individual

   >>> heritability = 0.7
   >>> operators.calculate_error_variance(example_pop, heritability)
   >>> operators.calculate_p(example_pop)
   >>> print(np.array(example_pop.indInfo('p')))
   [ -5.914  -9.91   -3.431   1.705  -6.228  -8.191  -6.529  -7.639  -7.721
      3.255  -4.3    -2.664   4.516  -6.719  -8.796  -7.221  -6.325   8.624
     -9.833   0.183   6.828  -8.382  -9.432  -4.054  -7.174  -1.427  -4.127
     -3.868  -2.605  -4.644   3.063  -4.283  -7.793  -2.364  -7.342  -3.79
      1.104  -3.4    -3.65   -2.31   -9.322  -1.742   1.93  -10.422  -5.688
     -3.107  -4.476  -5.138  -2.316   3.798   0.795  -6.71   -3.408  -9.865
     -0.851  -2.047  -4.579  -0.868  -9.356   0.209   2.896   3.036  -4.482
     -2.621  -4.892  -3.376  -5.189  -7.666  -2.429 -10.87   -0.462  -4.31
     -0.184  -1.023  -0.967  -7.608   0.757  -2.729  -6.495  -7.947  -5.493
      2.968  -1.071  -5.463  -8.298  -7.276  -4.92   -8.31   -3.426 -10.872
     -0.99   -6.105  -0.051  -5.167  -0.692  -0.158  -3.649   0.146  -6.078
     -0.35    0.057  -5.414  -7.03    1.191  -4.866]

.. _validating_the_calculate_g_function:

Validating the ``calculate_g`` Function
=======================================

Let's make sure that our function is correctly matching allele to its effect and
summing the effects correctly. We will look at the alleles individual ``1`` of
``example_pop`` at the QTL. Then we will sum the effects and compare the result
with our function :func:`calculate_g`.

.. code-block:: python
   :caption: Validating the calculation of ``g``

   >>> example_ind = example_pop.individual(0)
   >>> alpha_qtl_alleles = np.array(example_ind.genotype(ploidy=0))[qtl]
   >>> omega_qtl_alleles = np.array(example_ind.genotype(ploidy=1))[qtl]
   >>> example_g = [[], []]
   >>> for locus, alpha, omega in zip(qtl, alpha_qtl_alleles, omega_qtl_alleles):
   ...  print(locus, alpha, ae_array[locus, alpha], omega, ae_array[locus, omega])
   ...  example_g[0].append(ae_array[locus, alpha])
   ...  example_g[1].append(ae_array[locus, omega])
   6943 3 -0.8270650481760465 3 -0.8270650481760465
   14327 1 -0.5096871785660831 1 -0.5096871785660831
   16868 1 -0.8631600556024023 1 -0.8631600556024023
   17119 3 -0.7627590609820143 3 -0.7627590609820143
   35312 1 -0.388313399999193 1 -0.388313399999193
   >>> sum(example_g[0]) + sum(example_g[1])
   -6.701969486651478
   >>> example_pop.indByID(1).g
   -6.701969486651478
