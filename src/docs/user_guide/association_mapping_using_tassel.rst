.. _analysis_using_tassel:

#####################################################
Using Standalone Tassel for Multi-Locus Trait Mapping
#####################################################

This example will show how to perform association mapping of a quantitative
trait using tassel-5-standalone. tassel_5_standalone_ TASSEL requires certain types of input for
its mixed linear model mode.

.. _tassel_5_standalone: http://www.maizegenetics.net/tassel

* Hapmap formatted data
* Kinship Matrix
* Population Structure
* Trait

:py:mod:`saegus` has functions to calculate the (kinship) relationship matrix
using marker data [VanRaden2008]_. Population structure is computed using
[Patterson2006]_. We will give an example of hapmap formatted output. The
trait data is formatted by adding specific tags to the top of a trait vector.
:py:mod:`saegus` also writes an ``xml`` file for TASSEL as a config file. The
config file keeps the data grouped together appropriately and keeps you from
having to re-run the same commands over and over again.

.. code-block:: python
   :caption: Module imports

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='short', quiet=True)
   >>> import simuPOP as sim
   >>> import pandas as pd
   >>> import numpy as np
   >>> import random
   >>> from saegus import analyze, operators, parameters
   >>> np.set_printoptions(suppress=True, precision=5)



.. _setting_up_for_association_analysis:

Setting up the Population for Analysis
======================================

We will load the population and compute the required information as in previous
examples.

.. code-block:: python
   :caption: Preparing the population for TASSEL

   >>> example_pop = sim.loadPopulation('example_pop.pop')
   >>> example_pop.addInfoFields('ind_id', 'mother_id', 'father_id', 'g', 'p')
   >>> sim.tagID(example_pop)
   >>> sim.stat(example_pop, numOfSegSites=sim.ALL_AVAIL, vars=['segSites'])
   >>> segregating_loci = example_pop.dvars().segSites

We will use exponentially distributed mean ``1``) allele effects with 20 ``qtl``

.. code-block:: python
   :caption: Setting up additive parameterization

   >>> qtl = sorted(random.sample(segregating_loci, 20))
   >>> allele_effects_table = trait.construct_ae_table(example_pop, qtl, random.expovariate, 1)
   >>> allele_effects_array = trait.construct_ae_array(allele_effects_table, qtl)
   >>> heritability = 0.7
   >>> operators.calculate_g(example_pop, allele_effects_array)
   >>> operators.calculate_error_variance(example_pop, heritability)
   >>> operators.calculate_p(example_pop)
   >>> analyze.trait_formatter()


.. _hapmap_formatted_data:


Hapmap Formatted Data
=====================

Hapmap formatted data consists of 11 columns followed by the genotype data.

* rs#
* alleles
* chrom
* pos
* strand
* assembly#
* center
* protLSID
* assayLSID
* panelLSID
* QCcode

Each column corresponds to an individual. Each row corresponds to a locus.
Most of these columns we will fill with NA because we are only interested in
the association mapping aspect of TASSEL. In the hapmap format a genotype is
given as two upper case letters corresponding to a nucleotide or code.


A genotype is given in terms of two uppercase letters.

Code  Meaning
A  Adenine
C  Cytosine
G  Guanine
T  Thymine
R  A or G
Y  C or T
S  G or C
W  A or T
K  G or T
M  A or C
B  C or G or T
D  A or G or T
H  A or C or T
V  A or C or G
N  any base
. or -   gap

``saegus`` automatically does the conversion into hapmap format provided with a
``dict`` to code the allele of each individual into a letter.


.. _calculating_the_kinship_matrix:

Kinship Matrix
==============

The kinship matrix is calculated via the method given in VanRaden2008_. It
is the same method implemented in Synbreed. The marker allele is interpreted
as the minor allele. The elements of :math:`\vec{M}` are :math:`-1` for the minor
allele homozygote, :math:`0` for the heterozygote and :math:`1` for the
major allele homozygote.

:math:`n` The number of individuals
:math:`m` The number of loci


.. _calculating_population_structure:

Population Structure
====================

Population structure is used as a covariate. For the past examples the first
eigenvector explains the overwhelming majority of the variation. However,
you can check and compute a test statistic if desired.


.. _formatting_trait_data:

Trait
=====

Our simulated trait data is given a header and output as a text file.
TASSEL uses html style tags <tag> in the header to label the input of each
file.

.. code-block:: python
   :caption: Functions for handling trait data

   >>> allele_effects_array =
   >>> heritability = 0.7
   >>> operators.calculate_g(example_pop)
   >>> operators.calculate_error_variance(example_pop, heritability)
   >>> operators.calculate_p(example_pop)
   >>>

.. _tassel_config_file:

TASSEL Config File
==================

The final component is a ``xml`` file which specifies the protocol for
TASSEL to run. The config file can be in terms of relative or absolute paths
for its input. All TASSEL options can be specified in the config file.




.. [Patterson2006] Patterson, N, Price, A, Reich, D. (2006). Population Structure and Eigenanalysis. PLOS Genetics, 2(12). doi:10:1371/journal.pgen.0020190

.. [VanRaden2008] VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. Journal of Dairy Science, 91(11), 4414–23. doi:10.3168/jds.2007-0980
