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
   >>> from saegus import operators, parameters
   >>> np.set_printoptions(suppress=True, precision=5)


.. _overview_of_additive_trait_example:

Our goal is to show how to parameterize an additive trait in :mod`saegus` code.
We will make use of :mod:`simuPOP` ``infoFields`` to store the information about
the additive trait for each individual. Also we will add a unique identifier
for each individual which will help us verify that the trait is being calculated
correctly. In this example we will randomly choose quantitative trait loci from
among loci which are segregating in ``example_pop``.

Steps:

   + Load Population
   + Add Information Fields
   + Determine Segregating Loci
   + Choose QTL
   + Assign Allele Effects
   + Calculate G
   + Calculate error
   + Calculate P

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
   [1812,
    1905,
    4802,
    6092,
    7776,
    9225,
    11426,
    17994,
    18169,
    19480,
    21206,
    22754,
    27998,
    28313,
    29297,
    31358,
    36316,
    36354,
    40565,
    44143]

Every allele is assigned an effect of ``0``. Only the alleles at QTL have
non-zero effects.

.. code-block:: python
   :caption: Assign allele effects using an exponential distribution

   >>> trait = parameters.Trait()
   >>> ae_table = trait.construct_allele_effects_table(example_pop, qtl, random.expovariate, 1)
   >>> ae_table[qtl]
   [[  1812.         1.         0.069      3.         1.832]
    [  1905.         1.         0.192      3.         2.812]
    [  4802.         1.         0.009      3.         0.935]
    [  6092.         1.         3.329      2.         0.274]
    [  7776.         1.         0.885      3.         0.349]
    [  9225.         1.         0.018      2.         1.521]
    [ 11426.         1.         1.026      3.         0.223]
    [ 17994.         1.         0.374      2.         0.618]
    [ 18169.         1.         1.141      3.         0.688]
    [ 19480.         1.         6.983      3.         1.049]
    [ 21206.         1.         2.583      2.         0.173]
    [ 22754.         1.         1.162      3.         2.465]
    [ 27998.         1.         0.535      2.         1.631]
    [ 28313.         1.         4.603      3.         0.686]
    [ 29297.         1.         1.071      3.         0.001]
    [ 31358.         1.         2.123      3.         2.785]
    [ 36316.         1.         0.138      3.         0.951]
    [ 36354.         1.         0.465      2.         0.853]
    [ 40565.         1.         5.387      3.         0.006]
    [ 44143.         1.         1.22       3.         0.039]]
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
   [[ 0.     0.069  0.     1.832  0.   ]
    [ 0.     0.192  0.     2.812  0.   ]
    [ 0.     0.009  0.     0.935  0.   ]
    [ 0.     3.329  0.274  0.     0.   ]
    [ 0.     0.885  0.     0.349  0.   ]
    [ 0.     0.018  1.521  0.     0.   ]
    [ 0.     1.026  0.     0.223  0.   ]
    [ 0.     0.374  0.618  0.     0.   ]
    [ 0.     1.141  0.     0.688  0.   ]
    [ 0.     6.983  0.     1.049  0.   ]
    [ 0.     2.583  0.173  0.     0.   ]
    [ 0.     1.162  0.     2.465  0.   ]
    [ 0.     0.535  1.631  0.     0.   ]
    [ 0.     4.603  0.     0.686  0.   ]
    [ 0.     1.071  0.     0.001  0.   ]
    [ 0.     2.123  0.     2.785  0.   ]
    [ 0.     0.138  0.     0.951  0.   ]
    [ 0.     0.465  0.853  0.     0.   ]
    [ 0.     5.387  0.     0.006  0.   ]
    [ 0.     1.22   0.     0.039  0.   ]]

.. _definition_of_g:

Definition of ``g``
===================

``g`` is the sum of the allele effects of an individual's genotype. There is
no noise or error in ``g`` because we have *a priori* determined the allele
effects.

.. code-block:: python
   :caption: Calculating g values

   >>> operators.calculate_g(example_pop)
   >>> print(np.array(example_pop.indInfo('g')))
   array([ 40.5  ,  57.516,  42.954,  44.655,  58.748,  45.196,  44.301,
        37.803,  42.125,  48.263,  59.79 ,  46.791,  44.018,  40.228,
        46.464,  54.358,  50.271,  48.995,  49.538,  34.851,  43.836,
        47.706,  54.652,  40.614,  47.126,  48.786,  42.837,  42.593,
        54.974,  45.717,  44.98 ,  41.022,  47.093,  42.612,  47.278,
        46.156,  49.569,  45.891,  43.185,  46.977,  40.895,  39.624,
        46.451,  40.221,  41.131,  44.719,  46.342,  49.455,  42.355,
        49.107,  37.983,  46.371,  45.825,  49.369,  40.751,  42.464,
        48.045,  49.075,  47.905,  49.164,  46.342,  41.702,  41.419,
        45.088,  47.784,  48.206,  42.946,  46.279,  41.376,  48.122,
        40.604,  53.401,  43.177,  42.734,  40.98 ,  44.888,  46.668,
        43.456,  55.55 ,  43.821,  45.745,  40.688,  46.057,  44.673,
        49.514,  38.059,  40.034,  42.149,  40.867,  42.66 ,  49.946,
        44.809,  39.963,  46.583,  43.055,  49.495,  41.973,  46.353,
        43.615,  46.172,  39.211,  44.044,  44.618,  42.06 ,  43.291])

.. _plotting_the_distribution_of_g:

Plotting the Distribution of ``g``
----------------------------------

Let's visually inspect the distribution of ``g`` for this population.




.. _validating_the_calculate_g_function:

Validating the ``calculate_g`` Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

   >>> operators.calculate_p(example_pop)
   >>> print(np.array(example_pop.indInfo('p')))



Using a Normal Distribution Instead of Exponential
==================================================

Suppose we wanted to use a normal distribution for allele effects instead of
an exponential. All we need to do is change the parameter in the
``construct_allele_effects_table`` function.

.. code-block:: python
   :caption: Allele effects drawn from a normal distribution

   >>> normal_ae_table = trait.construct_allele_effects_table(example_pop, qtl, random.normalvariate, 0, 1)
   >>> print(normal_ae_table[qtl])
   [[  1812.         1.        -1.081      3.         0.317]
    [  1905.         1.         0.675      3.        -1.652]
    [  4802.         1.         0.307      3.        -1.259]
    [  6092.         1.         0.695      2.        -0.429]
    [  7776.         1.        -0.141      3.        -1.2  ]
    [  9225.         1.        -0.754      2.        -0.253]
    [ 11426.         1.        -0.499      3.        -1.067]
    [ 17994.         1.         0.804      2.         2.749]
    [ 18169.         1.        -0.354      3.         0.079]
    [ 19480.         1.         0.112      3.        -0.726]
    [ 21206.         1.        -0.812      2.         0.74 ]
    [ 22754.         1.        -0.125      3.         0.314]
    [ 27998.         1.        -1.239      2.         0.172]
    [ 28313.         1.         0.49       3.         1.02 ]
    [ 29297.         1.         1.022      3.         0.763]
    [ 31358.         1.         0.525      3.         0.563]
    [ 36316.         1.        -0.803      3.         0.73 ]
    [ 36354.         1.         0.266      2.        -2.607]
    [ 40565.         1.        -1.582      3.        -0.679]
    [ 44143.         1.         0.046      3.         1.264]]

Recomputing Using Normal Values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python
   :caption: Recalculate ``g``

   >>> normal_ae_array = trait.construct_ae_array(normal_ae_table, qtl)
   >>> operators.calculate_g(example_pop, normal_ae_array)
   >>> print(np.array(example_pop.indInfo('g')))
   [ -3.553  -9.525  -4.702  -4.797  -8.954   0.677  -0.047  -4.165  -6.304
     -1.938  -4.17    0.239  -5.376  -0.775  -3.369  -3.671  -4.242  -0.578
     -6.075  -6.511   0.25   -2.213  -2.302  -7.594  -3.914  -6.419  -3.559
      0.92  -10.755  -4.719   1.3    -1.734  -2.431  -4.007  -8.386   0.575
      0.719  -5.358  -3.105  -4.266  -5.877  -1.723  -3.222   2.485  -6.532
     -3.478  -5.369   1.964  -1.525  -0.737  -3.519  -8.021  -1.33   -2.929
     -0.985  -7.34   -4.304  -2.914  -1.826  -2.955  -2.134  -2.592  -7.036
     -4.123   0.51   -3.507   0.668   0.327  -2.461  -0.584   1.26   -6.559
     -7.789  -2.213  -6.319  -0.808  -4.924   0.751 -11.156  -5.651   0.903
      1.676  -1.173  -4.805  -0.773   4.606  -7.018   1.822  -0.15   -3.242
     -2.086  -1.359  -5.043   2.78   -2.491  -4.629  -3.859   2.17   -1.853
      1.854  -3.509  -3.715  -2.368   0.242   4.075]

