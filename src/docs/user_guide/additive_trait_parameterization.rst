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

   >>> print(np.array(example_pop.indInfo('example_pop')))
   [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
     0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
   >>> sim.tagID(example_pop)
   >>> print(np.array(example_pop.dvars().ind_id))
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
   >>> print(example_pop.dvars().segSites[::1000] # every 1000th segregating locus
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

For this example we will pick 20 loci to designate as quantitative trait loci.
The alleles at each chosen QTL will be assigned a non-zero effect via a draw
from an exponential distribution.

.. code-block:: python
   :caption: Choosing QTL and assigning allele effects

   >>> segregating_loci = example_pop.dvars().segSites
   >>> qtl = sorted(random.sample(segregating_loci, 20))
   >>> qtl
   [248,
    2609,
    4351,
    5444,
    5467,
    7902,
    8008,
    8951,
    10983,
    17571,
    17573,
    24130,
    25900,
    26640,
    30553,
    36841,
    38387,
    40501,
    42632,
    44217]

Every allele is assigned an effect of ``0``. Only the alleles at QTL have
non-zero effects.

.. code-block:: python
   :caption: Assign allele effects using an exponential distribution

   >>> test_run = analyze.Study('test')
   >>> allele_states = test_run.gather_allele_data(example_pop)
   >>> alleles = np.array([astates[:, 1], astates[:, 2]]).T
   >>> trait = parameters.Trait()
   >>> ae_table = trait.construct_allele_effects_table(alleles, qtl, random.expovariate, 1)
   >>> ae_table[qtl]
   [[   248.         2.         1.293      3.         2.876]
    [  2609.         1.         0.578      3.         1.497]
    [  4351.         1.         0.326      3.         0.024]
    [  5444.         1.         2.481      3.         0.247]
    [  5467.         2.         0.148      3.         1.209]
    [  7902.         1.         0.649      3.         0.868]
    [  8008.         1.         0.346      3.         0.341]
    [  8951.         1.         3.         3.         2.026]
    [ 10983.         1.         0.072      3.         1.709]
    [ 17571.         1.         0.586      3.         0.804]
    [ 17573.         1.         5.986      3.         1.062]
    [ 24130.         1.         1.623      2.         0.244]
    [ 25900.         1.         0.029      3.         0.074]
    [ 26640.         1.         1.266      3.         3.342]
    [ 30553.         1.         0.616      3.         0.278]
    [ 36841.         1.         0.14       3.         1.247]
    [ 38387.         1.         0.882      3.         0.669]
    [ 40501.         1.         2.083      3.         1.123]
    [ 42632.         1.         0.527      3.         0.029]
    [ 44217.         1.         0.703      2.         0.337]]
   >>> print(ae_table) # non-qtl
   [[     0.      1.      0.      2.      0.]
    [     1.      2.      0.      3.      0.]
    [     2.      2.      0.      3.      0.]
    ...,
    [ 44442.      1.      0.      2.      0.]
    [ 44443.      1.      0.      3.      0.]
    [ 44444.      1.      0.      3.      0.]]


For speed of computation we construct an array of allele effects where the row
of the array corresponds to the locus and the column corresponds to the integer
representing the allele state.

.. code-block:: python
   :caption: Putting the allele effects in an array for speed of computation

   >>> ae_array = trait.construct_ae_array(ae_table, qtl)
   >>> print(ae_array[qtl])
   [[ 0.     0.     1.293  2.876  0.   ]
    [ 0.     0.578  0.     1.497  0.   ]
    [ 0.     0.326  0.     0.024  0.   ]
    [ 0.     2.481  0.     0.247  0.   ]
    [ 0.     0.     0.148  1.209  0.   ]
    [ 0.     0.649  0.     0.868  0.   ]
    [ 0.     0.346  0.     0.341  0.   ]
    [ 0.     3.     0.     2.026  0.   ]
    [ 0.     0.072  0.     1.709  0.   ]
    [ 0.     0.586  0.     0.804  0.   ]
    [ 0.     5.986  0.     1.062  0.   ]
    [ 0.     1.623  0.244  0.     0.   ]
    [ 0.     0.029  0.     0.074  0.   ]
    [ 0.     1.266  0.     3.342  0.   ]
    [ 0.     0.616  0.     0.278  0.   ]
    [ 0.     0.14   0.     1.247  0.   ]
    [ 0.     0.882  0.     0.669  0.   ]
    [ 0.     2.083  0.     1.123  0.   ]
    [ 0.     0.527  0.     0.029  0.   ]
    [ 0.     0.703  0.337  0.     0.   ]]

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
   [ 47.551  43.782  41.252  45.229  43.477  44.363  46.361  45.871  38.205
     44.067  40.832  48.246  35.99   48.896  43.381  44.006  40.275  39.4
     42.295  42.297  49.207  47.154  40.323  44.564  45.493  46.109  48.308
     48.786  37.21   40.446  43.844  37.579  45.431  38.584  48.469  44.869
     43.905  45.335  46.453  40.682  42.834  42.491  47.074  49.875  44.902
     47.802  42.042  43.454  43.898  43.955  34.524  47.52   42.267  44.827
     45.142  47.934  43.975  46.465  47.375  40.209  40.356  45.553  50.139
     49.649  36.133  41.16   36.637  47.069  45.871  45.299  37.93   41.483
     40.249  47.552  43.699  41.168  44.15   48.072  40.277  42.204  43.874
     44.48   39.389  45.467  44.215  45.45   46.685  43.067  34.726  45.275
     42.136  42.337  36.428  46.922  42.695  44.359  49.863  43.432  41.794
     40.63   43.844  39.743  44.324  44.548  44.843]


.. _calculating_error:

Calculation of Error Term
=========================

To simulate the experimental noise a term :math:`\epsilon` is added to each
individual's ``g`` value.
:math:`\epsilon` is a random variable with a normal distribution given by
mean :math:`0` and variance given by:

.. math::

   \sigma^2_g = V_g * (\frac{1}{h^2} - 1)

Where :math:`V_g` is the variance of ``g`` and :math:`h^2` is the
narrow sense heritability.


.. math::

   \varepsilon \sim \mathcal{N} (0, \sigma^2_g)

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
   [ 51.912  43.698  31.867  39.772  42.367  55.394  45.751  36.806  30.67
     49.342  39.44   53.184  34.946  47.489  51.119  45.355  49.242  42.428
     43.989  41.591  40.911  51.125  51.154  44.06   38.491  45.406  48.175
     44.659  31.458  35.266  46.712  41.525  53.201  45.575  54.09   48.758
     37.973  49.883  49.092  39.151  45.827  47.702  46.029  54.084  42.357
     52.86   37.865  49.702  39.409  36.099  34.894  43.194  44.701  41.302
     48.899  51.28   41.661  44.914  41.055  47.     38.409  44.145  49.102
     43.193  39.99   41.011  41.165  56.536  52.146  41.167  41.825  29.432
     40.996  46.574  41.725  31.878  47.394  48.373  35.343  46.933  44.161
     37.527  40.506  44.     47.073  42.077  39.873  36.111  33.092  43.489
     39.986  46.277  36.09   39.514  38.24   47.853  50.019  49.436  43.239
     40.166  45.012  43.531  43.878  53.475  45.172]


Using a Normal Distribution Instead of Exponential
==================================================

Suppose we wanted to use a normal distribution for allele effects instead of
an exponential. All we need to do is change the parameter in the
``construct_allele_effects_table`` function.

.. code-block:: python
   :caption: Allele effects drawn from a normal distribution

   >>> normal_ae_table = trait.construct_allele_effects_table(example_pop, qtl, random.normalvariate, 0, 1)
   >>> print(normal_ae_table[qtl])
   [[   248.         2.         0.414      3.        -0.983]
    [  2609.         1.         0.252      3.        -1.557]
    [  4351.         1.        -0.312      3.        -0.314]
    [  5444.         1.        -0.998      3.         1.352]
    [  5467.         2.         0.939      3.        -1.758]
    [  7902.         1.         1.075      3.        -1.841]
    [  8008.         1.         0.554      3.         1.327]
    [  8951.         1.         0.051      3.         0.1  ]
    [ 10983.         1.         0.976      3.         1.786]
    [ 17571.         1.        -0.754      3.         1.578]
    [ 17573.         1.         0.561      3.        -0.436]
    [ 24130.         1.         1.027      2.        -0.745]
    [ 25900.         1.        -0.363      3.        -2.082]
    [ 26640.         1.         0.503      3.        -1.156]
    [ 30553.         1.         0.605      3.        -0.605]
    [ 36841.         1.        -0.568      3.         2.009]
    [ 38387.         1.        -0.605      3.        -0.625]
    [ 40501.         1.        -1.32       3.         0.672]
    [ 42632.         1.         0.715      3.        -0.453]
    [ 44217.         1.         0.346      2.         1.045]]

Recomputing Using Normally Distributed Allele Effects
-----------------------------------------------------

.. code-block:: python
   :caption: Recalculate ``g``

   >>> normal_ae_array = trait.construct_ae_array(normal_ae_table, qtl)
   >>> operators.calculate_g(example_pop, normal_ae_array)
   >>> print(np.array(example_pop.indInfo('g')))
   [  9.391   8.02    4.225   0.09   11.742   7.926  11.427  15.576  -0.911
      6.048   9.484   6.613  -3.472   3.755  14.837   0.041   9.211   8.871
      6.226   0.924   6.576  15.055   1.922   9.408   1.357   0.158   5.3
     10.485   8.384  11.894   0.574   6.921   5.443  -3.902  10.89    3.441
     10.557   6.855  13.957   6.83    0.855   7.258   5.761  11.502   9.072
     10.946   5.282   6.848  10.172   6.145   0.627   2.356   3.401  10.021
     -0.208   8.668   2.678   2.451   4.81    3.761   3.486   5.8    11.652
      3.193   5.9     2.796   6.423   5.589  10.268  14.335   4.298   2.229
      6.918   3.191   4.046   2.462  -0.218   6.304   3.416   7.002   7.478
      4.093   7.007  -1.809   5.602   6.77    8.306  14.293   0.037   7.616
      9.862   3.128   8.275   4.544   8.11    0.443   2.44    8.743   5.935
      0.777  12.231   9.319   2.053   1.085  13.271]

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
   >>> sum(example_g[0]) + sum(example_g[1])
   40.500306681374511
   >>> example_pop.indByID(1).g
   40.500306681374504
