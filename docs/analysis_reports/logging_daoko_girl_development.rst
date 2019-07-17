
.. _run-daoko-girl:

======================================
Log of Development for Run: daoko_girl
======================================

.. _daoko-girl-parameters:

Location and Definition of Parameters
=====================================

Genetic Map
~~~~~~~~~~~

A ``genetic_map`` determines how many chromosomes and how many loci are on each
chromosome. The ``genetic_map`` also defines the rates of recombination between
sequential loci as defined by simuPOP_.

.. _simuPOP: http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec5.html#recombination-operator-recombinator


Synthesis Parameters
~~~~~~~~~~~~~~~~~~~~

The purpose of the ``synthesis_parameters`` is to tell the user exactly how the population
was created. The contents of ``synthesis_parameters`` should be able to tell you
where the genotypes of the founders came from and how the founders were mated with each
other. The following items are defined in ``synthesis_parameters``:

   * ``prefounder_file_name``
   * ``prefounder_names``
   * ``founders``
   * ``integer_to_snp``
   * ``snp_to_integer``
   * ``mating_scheme``
   * ``development_population_size``

Trait Parameters
~~~~~~~~~~~~~~~~

``saegus`` was developed with the intention of simulating quantitative phenotypes.
``trait_parameters`` should tell the user exactly how the phenotypes, ``p``, were
calculated. At present ``trait_parameters`` has the following entries:

   * ``number_of_qtl``
   * ``qtl``
   * ``allele_effect_distribution``
   * ``allele_effect_parameters``
   * ``allele_effects``
   * ``heritability``
   * ``epsilon``

Alleles
~~~~~~~

The alleles which are present at each locus for the ``magic_1478.pop`` population
can be found in``alleles_at_1478_loci.hdf``.

Analysis Parameters
~~~~~~~~~~~~~~~~~~~

The ability to map the output from TASSEL to the input from ``saegus`` is crucial.
The ``analysis_parameters`` record information which allows the user to make the
connection between input and output.

   * ``population_name``
   * ``scenario``
   * ``generations``
   * ``operating_population_size``
   * ``output_prefix``
   * ``sample_size``
   * ``sample_ind_ids``
   * ``sample_allele_frequencies``
   * ``segregating_loci``
   * ``saegus_to_tassel_loci``
   * ``tassel_to_saegus_loci``

Intermediate Data
~~~~~~~~~~~~~~~~~

There is a great deal of data and parameters moving around inside even a single
simulation of ``saegus``. It is wasteful to generate new data each time I work on
``saegus``. ``intermediate_data`` records parameters or data that in the case
I take a break from developing. ``intermediate_data`` can also be used for debugging
and metadata.

.. code-block:: python

   intermediate_data = shelve.open('daoko_girl_debug_data')
   intermediate_data['allele_frequencies'] = af
   intermediate_data['segregating_allele_frequencies'] = segregating_frame
   intermediate_data['g'] = np.array(magic1478.indInfo('g'))
   intermediate_data['p'] = np.array(magic1478.indInfo('p'))
   intermediate_data['segregating_loci'] = segregating_loci
   intermediate_data['run_name'] = 'daoko_girl'
   intermediate_data.close()



Running the Simulation and Formatting Output for TASSEL
=======================================================

.. code-block:: shell

   $ pwd
   C:\Users\DoubleDanks\wisser\code\rjwlab-scripts\saegus_project\devel\magic\1478

I loaded the **standard** population with 1478 loci made from prefounders 1 through 8.

.. code-block:: python

   magic1478 = sim.loadPopulation('populations\\magic_1478.pop')
   magic1478.dvars()

   {'rep': 0, 'gen': 3}


I loaded the standard genetic map of 1478 loci with recombination rates of 0.01
between sequential loci. The genetic map is saved in a ``.hdf`` file and is loaded
using :py:func:`pd.read_hdf`

.. code-block:: python

   genetic_map = pd.read_hdf('parameters\\genetic_map_1478.hdf)
   genetic_map

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>chromosome</th>
          <th>relative_locus</th>
          <th>recom_rate</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>1</td>
          <td>0</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>1</th>
          <td>1</td>
          <td>1</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>2</th>
          <td>1</td>
          <td>2</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>3</th>
          <td>1</td>
          <td>3</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>4</th>
          <td>1</td>
          <td>4</td>
          <td>0.01</td>
        </tr>
        <tr>
          <th>5</th>
          <td>1</td>
          <td>5</td>
          <td>0.01</td>
        </tr>
        <tr>
        <td>...</td>
        <td>...</td>
        <td>...</td>
        <td>...</td>
        </tr>
    </table>
    <p>1478 rows × 3 columns</p>
 </div>


Creating the Population to Analyze
==================================

Overview
~~~~~~~~

The population from ``magic_1478.pop`` undergoes random mating for three generations
with a fixed population size. After random mating each individual is assigned ``g``
and ``p``. The error term for each individual is a random draw from a normal distribution
with mean :math:`0` and variance equal to:

   Heritability, :math:`h^2`, is set to :math:`0.7` and the genotypic variance is :math:`V_g`

   .. math::

      e = V_g * (1 / h^2 - 1)

Once the phenotype data is assigned we take a random sample from our analysis population
and use the methods of the class :class:`GWAS` in :mod:`analyze` to generate
the input for TASSEL.

Executing Random Mating
~~~~~~~~~~~~~~~~~~~~~~~

I created a class in :mod:`breed` to handle mating. :class:`MAGIC` is initialized
with two parameters: population and recombination rates. The following block shows
an example of how to use the :class:`MAGIC` to perform random mating.

.. code-block:: python

   sim.tagID(base_population, reset=False)
   random_mater = breed.MAGIC(base_population, recombination_rates)
   random_mater.interim_random_mating(analysis_parameters['generations'],
                              analysis_parameters['operating_population_size'])

.. code-block:: none

   Initiating interim random mating for 3 generations.
   Generation: 3
   Generation: 4
   Generation: 5

.. note:: In this run we are assuming a simple additive trait model.

QTL and Allele Effects
~~~~~~~~~~~~~~~~~~~~~~

Once we have our randomly intermated population I use a ``simuPOP`` operator
to determine which loci are segregating. Only the loci that are segregating after
random mating are eligible to be defined as a quantitative trait locus.

.. code-block:: python

   sim.stat(base_population, numOfSegSites=sim.ALL_AVAIL, vars=['segSites'])
   qtl = sorted(random.sample(base_population.dvars().segSites, number_of_qtl))

For this run we assign an allele effect as three independent draws from
an exponential distribution with :math:`theta` equal to one. The parameter ``alleles``
is a list of the alleles present at each locus in the prefounders. This guarantees
that every allele is assigned an effect even in the case we decide to analyze a fixed
locus.

.. code-block:: python

   additive_trait = parameters.Trait()
   allele_effects = additive_trait.assign_allele_effects(alleles, qtl, random.expovariate, 1, multiplicity=3)

.. note:: ``saegus`` was developed with two alleles at each locus. Later versions will generalize these functions.

Given a ``dict`` of ``allele_effects`` we calculate ``g`` and ``p`` for each
individual using the function form of one of the ``saegus`` operators.

.. code-block:: python

   operators.assign_additive_g(base_population, qtl, allele_effects)
   operators.calculate_error_variance(base_population, heritability)
   operators.phenotypic_effect_calculator(base_population)

Now we have everything we need to create the input for TASSEL.

Generating TASSEL Input
=======================

The class :class:`GWAS` has methods to handle all of the formatting, calculations and
output for the four files we use for TASSEL's mixed linear model. An external
function calls the appropriate methods in the appropriate order



.. _analyze-magic1478-rdm-mating-results:

Analyzing Results of GWAS with TASSEL
=====================================

We use the mixed-linear modeling method implemented in TASSEL.
The MLM in TASSEL requires three files at minimum with the option for a fourth.

   * genotypes in *hapmap* format
   * kinship matrix
   * phenotypes
   * population structure matrix (optional)












QVALUES in R
============

We will follow Jim's tutorial to use the :mod:`qvalue` package in R; however, I
have found that the function we want to use :func:`qvalue` does not handle
missing data i.e. ``NaN``. Because I am more proficient with ``python`` than
``R`` I used the ``python`` :mod:`pandas` package to convert all ``NaN`` p-values
into values of :math:`1.0`

For example a sample of the P-values of ``gwas_out_2.txt`` are:

+ 0.4968
+ 5.6091E-28
+ NaN
+ 0.6236
+ 0.16525


If we use the :func:`qvalue` function directly it will result in an error.
Instead I use the values:

+ 0.4968
+ 5.6091E-28
+ 1.0
+ 0.6236
+ 0.16525

The edited file name is ``edited_gwas_out_2.txt``. I use these commands to
obtain the q-values.

.. code-block::

   results_header = scan("edited_gwas_out_2.txt", what="character", nlines=1, sep="\t")
   gwas_results = read.table("edited_gwas_out_2.txt", header=F, row.names=NULL, skip=2)
   colnames(gwas_results) = results_header

   pvalues = gwas_results$p
   library(qvalue)
   qobj = qvalue(p = pvalues)
   qobj$qvalues
   qvalues_of_magic1478_results = data.frame(qobj$qvalues)
   write.table(qvalues_of_magic1478_results, "qvalues_of_magic1478.txt", sep="\t")

We use the :func:`qqunif` function in R to produce the quantile-quantile plot of
the p-values.

.. figure:: qqplot.png

   Quantile-Quantile plot.

TASSEL detected two of the three QTL: 2 and 20.
The Q-values are 4.1451249e-25 and 1.606586e-73 respectively.

Two observations that are immediately obvious:

   1) there is almost no difference between allele effects at locus 10
   2) locus 10 is very close to 2 and 20 which might mute its already small effect


