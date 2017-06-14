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
   >>> from saegus import breed, operators, simulate, analyze, parse, parameters
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
   + Calculate G and P

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

There are 42,837 segregating loci in this population. Next we will gather the
alleles which are present at each segregating locus.

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

Choosing QTL
~~~~~~~~~~~~

For this example we will pick 20 loci to designate as quantitative trait loci.
The alleles at each chosen QTL will be assigned a non-zero effect via a draw
from an exponential distribution.

.. code-block:: python
   :caption: Choosing QTL and assigning allele effects

   >>> qtl = sorted(random.sample(segregating_loci, 20))
   >>> qtl
   [5027,
    7313,
    11571,
    13436,
    15145,
    15615,
    17727,
    17946,
    18912,
    22551,
    23076,
    26364,
    30497,
    31261,
    34355,
    34668,
    37124,
    37753,
    37920,
    40366]
