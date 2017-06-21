.. _example_of_additive_quantitative_trait:

##########################################
Example of Additive Trait Parameterization
##########################################

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

Overview
========

Our goal is to show how to parameterize an additive trait in ``saegus`` code.
We will make use of ``simuPOP`` ``infoFields`` to store the information about
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

.. _load_population:

Load Population
^^^^^^^^^^^^^^^

We will use the population we created in the last step instead of creating
a new population.

.. code-block:::: python
   :caption: Loading our example population from a file

   >>> sim.loadPopulation('example_pop.pop')

.. _add_information_fields:

Add Information Fields
~~~~~~~~~~~~~~~~~~~~~~

At present our population has no special information assigned to its members.
Each individual is only a genotype. By default ``simuPOP`` uses ``ind_id``,
``father_id``, ``mother_id``, ``father_idx`` and ``mother_idx`` for its very
useful ``IdTagger`` functions. We can save some hassle by using these for
identifying individuals.

.. code-block:: python
   :caption: Add the ``infoFields`` ``ind_id``, ``g`` and ``p``

   >>> example_pop.infoFields()
   ()
   >>> example_pop.addInfoFields(['ind_id', 'g', 'p'])
   >>> example_pop.infoFields()
   ('ind_id', 'g', 'p')

By default information fields are set to ``0.0``. We can initialize the
``ind_id`` field using a ``simuPOP`` function.

.. code-block:: python
   :caption: Initialize individual identifiers

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

.. note:: In this step we converted the output into a np.array for aesthetics

.. _determine_segregating_loci:

Determine Segregating Loci
~~~~~~~~~~~~~~~~~~~~~~~~~~

For simplicity we will loci which have more than one allele i.e. segregating.
It will be useful to extract the alleles from each locus for later use.

.. code-block:: python
   :caption: Using ``simuPOP`` to find segregating loci

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




.. code-block:: python
   :caption: Gather the alleles at each segregating site

   >>> sim.stat(example_pop, alleleFreq=sim.ALL_AVAIL)
   >>> segregating_loci = example_pop.dvars().segSites
   >>> alpha_alleles = []
   >>> beta_alleles = []
   >>> for locus in segregating_loci:
   ...      alpha_alleles.append(list(example_pop.dvars().alleleFreq[locus])[0])
   ...      beta_alleles.append(list(example_pop.dvars().alleleFreq[locus])[1])

We have the alleles at each segregating site in two separate Python lists
i.e. ``alpha_alleles`` and ``beta_alleles``. Let's check to make sure that all
entries in ``alpha_alleles`` are different from ``beta_alleles``.

.. code-block:: python
   :caption: A quick check to see if our code is semantically correct

   >>> alpha_allele_array = np.array(alpha_alleles)
   >>> beta_allele_array = np.array(beta_alleles)
   >>> sum(alpha_allele_array == beta_allele_array)
   0

Because the result is ``0`` that means that every entry of ``alpha_alleles`` is
different from ``beta_alleles``.

.. _choose_QTL:

Choosing QTL and Assign Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For this example we will pick 20 loci to designate as quantitative trait loci.
The alleles at each chosen QTL will be assigned a non-zero effect via a draw
from an exponential distribution. We are choosing QTL from only
segregating loci.

.. code-block:: python
   :caption: Choosing QTL and assigning allele effects

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
   :caption: Assign allele effects as an exponential distribution

   >>> trait = parameters.Trait()
   >>> ae_table = trait.construct_allele_effects_table(example_pop, qtl, random.expovariate, 1)
   >>> ae_table[qtl]
   array([[  1812.   ,      1.   ,      2.559,      3.   ,      1.962],
          [  1905.   ,      1.   ,      0.169,      3.   ,      0.199],
          [  4802.   ,      1.   ,      0.533,      3.   ,      0.523],
          [  6092.   ,      1.   ,      0.5  ,      2.   ,      4.702],
          [  7776.   ,      1.   ,      1.825,      3.   ,      0.156],
          [  9225.   ,      1.   ,      0.793,      2.   ,      1.657],
          [ 11426.   ,      1.   ,      1.064,      3.   ,      0.228],
          [ 17994.   ,      1.   ,      0.221,      2.   ,      0.015],
          [ 18169.   ,      1.   ,      1.011,      3.   ,      1.45 ],
          [ 19480.   ,      1.   ,      1.443,      3.   ,      0.046],
          [ 21206.   ,      1.   ,      0.554,      2.   ,      1.086],
          [ 22754.   ,      1.   ,      0.904,      3.   ,      0.628],
          [ 27998.   ,      1.   ,      0.361,      2.   ,      0.023],
          [ 28313.   ,      1.   ,      1.953,      3.   ,      0.033],
          [ 29297.   ,      1.   ,      2.737,      3.   ,      3.567],
          [ 31358.   ,      1.   ,      0.778,      3.   ,      1.601],
          [ 36316.   ,      1.   ,      6.54 ,      3.   ,      2.131],
          [ 36354.   ,      1.   ,      0.573,      2.   ,      1.766],
          [ 40565.   ,      1.   ,      0.137,      3.   ,      0.351],
          [ 44143.   ,      1.   ,      0.338,      3.   ,      0.719]])

For speed of computation we construct an array of allele effects where the row
of the array corresponds to the locus and the column corresponds to the integer
representing the allele state.

.. code-block:: python
   :caption: Putting the allele effects in an array for speed of computation

   >>> ae_array = trait.construct_ae_array(ae_table, qtl)
   >>> ae_array[qtl]
   array([[ 0.   ,  2.559,  0.   ,  1.962,  0.   ],
       [ 0.   ,  0.169,  0.   ,  0.199,  0.   ],
       [ 0.   ,  0.533,  0.   ,  0.523,  0.   ],
       [ 0.   ,  0.5  ,  4.702,  0.   ,  0.   ],
       [ 0.   ,  1.825,  0.   ,  0.156,  0.   ],
       [ 0.   ,  0.793,  1.657,  0.   ,  0.   ],
       [ 0.   ,  1.064,  0.   ,  0.228,  0.   ],
       [ 0.   ,  0.221,  0.015,  0.   ,  0.   ],
       [ 0.   ,  1.011,  0.   ,  1.45 ,  0.   ],
       [ 0.   ,  1.443,  0.   ,  0.046,  0.   ],
       [ 0.   ,  0.554,  1.086,  0.   ,  0.   ],
       [ 0.   ,  0.904,  0.   ,  0.628,  0.   ],
       [ 0.   ,  0.361,  0.023,  0.   ,  0.   ],
       [ 0.   ,  1.953,  0.   ,  0.033,  0.   ],
       [ 0.   ,  2.737,  0.   ,  3.567,  0.   ],
       [ 0.   ,  0.778,  0.   ,  1.601,  0.   ],
       [ 0.   ,  6.54 ,  0.   ,  2.131,  0.   ],
       [ 0.   ,  0.573,  1.766,  0.   ,  0.   ],
       [ 0.   ,  0.137,  0.   ,  0.351,  0.   ],
       [ 0.   ,  0.338,  0.   ,  0.719,  0.   ]])

Then we calculate ``g``: the value corresponding to the alleles of an individual
without any noise or error.

.. code-block:: python
   :caption: Calculating g values

   >>> operators.calculate_g(example_pop)
   >>> np.array(example_pop.indInfo('g'))
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

.. _validating_the_calculate_g_function:

Validating the ``calculate_g`` Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's make sure that our function is correctly matching allele to its effect and
summing the effects correctly. We will look at the alleles individual ``1`` of
``example_pop`` at the QTL. Then we will sum the effects and compare the result
with our function ``calculate_g``.

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