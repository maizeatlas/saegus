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
   [ 0.  0.  0.  0.  0. ... 0.]
   >>> sim.tagID(example_pop)
   >>> print(np.array(example_pop.indInfo('ind_id')))
   [   1.    2.    3.    4.    5. ... 105.]

.. note::
   ::
   In this step we converted the output into a np.array for aesthetics
   
.. _calculate_allele_frequencies:

Calculate Allele Frequencies
==========================

.. code-block:: python
   :caption: Using :mod:`simuPOP` to compute allele frequencies
   
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
   [0, 1040, 2072, 3098, 4124, ... 43578]

There are 42,837 segregating loci in this population. ``saegus`` has a function
to put the alleles into an array and assign the alleles at ``qtl`` an effect as
a draw from a specified distribution.

.. _additive_trait:

Additive Trait
##############

We have all the information we need from the previous steps. We will randomly
choose ``5`` QTL from the segregating loci. Both alleles at each QTL are
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
   [7790, 21801, 22978, 29480, 30705]

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
   [[ 7790.          1.          1.69685     3.          0.02152]
    [21801.          1.          0.9653      3.          0.30436]
    [22978.          1.          0.21302     3.          0.22749]
    [29480.          1.          1.58062     3.          0.57265]
    [30705.          1.          0.18288     3.          0.04947]]
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
   [[ 7790.          1.         -0.41871     3.         -0.49912]
    [21801.          1.          0.5188      3.         -1.04793]
    [22978.          1.          0.9032      3.          0.49079]
    [29480.          1.          0.24459     3.          0.64607]
    [30705.          1.         -0.88107     3.         -0.80171]]

For speed of computation we construct an array of allele effects where the row
of the array corresponds to the locus and the column corresponds to the integer
representing the allele state.

.. code-block:: python
   :caption: Putting the allele effects in an array for speed of computation

   >>> ae_array = trait.construct_ae_array(ae_table, qtl)
   >>> print(ae_array[qtl])
   [[ 0.      -0.41871  0.      -0.49912  0.       0.     ]
    [ 0.       0.5188   0.      -1.04793  0.       0.     ]
    [ 0.       0.9032   0.       0.49079  0.       0.     ]
    [ 0.       0.24459  0.       0.64607  0.       0.     ]
    [ 0.      -0.88107  0.      -0.80171  0.       0.     ]]

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
   [-0.83311  0.40057 -0.84405 -1.65794 -2.73393 ...  -2.82318]

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
   [-0.35134  0.27554 -1.03809 -1.68467 -3.88023 ... -1.95743]

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
   7790 1 -0.41871386917093 1 -0.41871386917093
   21801 3 -1.047928786709147 1 0.5187969288611954
   22978 1 0.9032044301593078 1 0.9032044301593078
   29480 1 0.24459159812004574 1 0.24459159812004574
   30705 1 -0.8810726262417609 1 -0.8810726262417609
   >>> sum(example_g[0]) + sum(example_g[1])
   -0.8331127921146263
   >>> example_pop.indByID(1).g
   -0.833112792114626
