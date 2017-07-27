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

+------+-------------+
| Code |   Meaning   |
+======+=============+
|   A  |   Thymine   |
+------+-------------+
|   C  |   Cytosine  |
+------+-------------+
|   G  |   Guanine   |
+------+-------------+
|   T  |   Thymine   |
+------+-------------+
|   R  |    A or G   |
+------+-------------+
|   Y  |    C or T   |
+------+-------------+
|   S  |    G or C   |
+------+-------------+
|   W  |    A or T   |
+------+-------------+
|   K  |    G or T   |
+------+-------------+
|   M  |    A or C   |
+------+-------------+
|   B  | C or G or T |
+------+-------------+
|   D  | A or G or T |
+------+-------------+
|   H  | A or C or T |
+------+-------------+
|   V  | A or C or G |
+------+-------------+
|   N  |   any base  |
+------+-------------+
|   .  |     gap     |
+------+-------------+
|   /  |     gap     |
+------+-------------+


``saegus`` automatically does the conversion into hapmap format provided with a
``dict`` to code the allele of each individual into a letter. For these
examples the allele states (``1, 2, 3, 4``) have nothing to do with the
corresponding nucleotides of the hapmap format. The :py:class:`GWAS`
automatically creates an attribute ``individual names`` which formats the
``ind_id`` for TASSEL.

.. code-block:: python
   :caption: Formatting genotype matrix in hapmap format

   >>> gwas = analyze.GWAS(example_pop, segregating_loci, 'example')
   >>> print(gwas.individual_names)
   ['I1' 'I2' 'I3' 'I4' 'I5' 'I6' 'I7' 'I8' 'I9' 'I10' 'I11' 'I12' 'I13' 'I14'
    'I15' 'I16' 'I17' 'I18' 'I19' 'I20' 'I21' 'I22' 'I23' 'I24' 'I25' 'I26'
    'I27' 'I28' 'I29' 'I30' 'I31' 'I32' 'I33' 'I34' 'I35' 'I36' 'I37' 'I38'
    'I39' 'I40' 'I41' 'I42' 'I43' 'I44' 'I45' 'I46' 'I47' 'I48' 'I49' 'I50'
    'I51' 'I52' 'I53' 'I54' 'I55' 'I56' 'I57' 'I58' 'I59' 'I60' 'I61' 'I62'
    'I63' 'I64' 'I65' 'I66' 'I67' 'I68' 'I69' 'I70' 'I71' 'I72' 'I73' 'I74'
    'I75' 'I76' 'I77' 'I78' 'I79' 'I80' 'I81' 'I82' 'I83' 'I84' 'I85' 'I86'
    'I87' 'I88' 'I89' 'I90' 'I91' 'I92' 'I93' 'I94' 'I95' 'I96' 'I97' 'I98'
    'I99' 'I100' 'I101' 'I102' 'I103' 'I104' 'I105']
   >>> hapmap_columns = ['rs', 'alleles', 'chrom', 'pos',
   ...              'strand', 'assembly', 'center',
   ...              'center', 'protLSID', 'assayLSID',
   ...              'panelLSID', 'QCode'] + list(gwas.individual_names)
   ...              ]
   >>> hapmap_matrix = pd.DataFrame(columns=hapmap_columns)
   >>> hapmap.rs = segregating_loci
   >>> hapmap.alleles = segregating_minor_alleles
   >>> chromosomes = np.array([example_pop.chromLocusPair(locus)[0] + 1
   ...                          for locus in segregating_loci], dtype=np.int8)
   >>> hapmap.chrom = chromosomes
   >>> hapmap_matrix.pos = np.arange(segregating_loci.shape)
   >>> hapmap_matrix.loc[:, 'strand':'QCode'] = np.core.defchararray.array(
   ...                          [['NA']*42837]*8).T
   >>> for i, ind in enumerate(example_pop.individuals()):
   ...       hapmap_matrix.loc[:, gwas.individual_names[i]] = [
   ...           ''.join(sorted(gwas.int_to_snp_conversions[a] +
   ...                     gwas.int_to_snp_conversions[b]))
   ...           for a, b, in zip(
   ...              np.array(ind.genotype(ploidy=0))[segregating_loci],
   ...              np.array(ind.genotype(ploidy=1))[segregating_loci])
   ...               ]
   >>> print(np.array(hapmap_matrix))
   [[0 1 1 ..., 'GG' 'GG' 'GG']
    [1 2 1 ..., 'GT' 'GG' 'TT']
    [2 3 1 ..., 'GG' 'GG' 'GG']
    ...,
    [44442 2 10 ..., 'GG' 'CG' 'CG']
    [44443 3 10 ..., 'CT' 'CC' 'CT']
    [44444 1 10 ..., 'CC' 'TT' 'TT']]

The final step is to write the hapmap to a ``txt`` file with the appropriate
rows and columns.

.. code-block:: python
   :caption: Writing the hapmap file for TASSEL

   >>> with open('example_hapmap.txt', 'w') as hapmap_file:
   ...         hapmap_matrix.to_csv(hapmap_file, sep='\t', index=False)

All of the operations to format the genotypes into hapmap format are
encapsulated inside of :func:`hapmap_formatter`.

.. _calculating_the_kinship_matrix:

Kinship Matrix
==============

The kinship matrix is calculated via the method given in VanRaden2008_. It
is the same method implemented in Synbreed. The marker allele is interpreted
as the minor allele. The elements of :math:`\mathbf{V_{n x m}}` are
:math:`-1` for the minor allele homozygote, :math:`0` for the heterozygote
and :math:`1` for the major allele homozygote.

* :math:`n` The number of individuals
* :math:`m` The number of loci
* :math:`p_i` Frequency of the minor allele at locus i

We can obtain :math:`\mathbf{V}` by multiplying :math:`\mathbf{C}` by
:math:`-1` and adding :math:`1`. Recall that
:math:`\mathbf{C}` is given by :func:`calculate_count_matrix`
in terms of minor alleles.

.. math::

   \left[
   \begin{array}
   ((-1)*2 + 1 = (-1) \\
   (-1)*1 + 1 = 0 \\
   (-1)*0 + 1 = 1
   \end{array}
   \right]

Which is exactly what is required. We use method (1) to compute
:math:`\mathbf{G}`. We only need to add ``individual_names`` to the header
of :math:`\mathbf{G}` for use in TASSEL.


.. math::

   \mathbf{P} = 2(p_i - \frac{1}{2})

.. math::

   \mathbf{Z}_{n \times m} = \mathbf{V_j} - \mathbf{P}

.. math::

   \mathbf{G} = \frac{\mathbf{Z}\mathbf{Z^T}}{2\sum_{i}p_i(1-p_i)}



.. _calculating_population_structure:

Population Structure
====================

Population structure is used as a covariate. For the past examples the first
eigenvector explains the overwhelming majority of the variation. However,
you :py:mod:`saegus` has functions to compute a test statistic. The value
of the test statistic must be compared against the *Tracy-Widom*
manually as there is not a distribution implemented in Python. We implement
the computation in Patterson2006_. Let:

* :math:`m` be the number of individuals
* :math:`n` the number of loci (called markers in the paper)
* :math:`a` is the minor allele at each locus
* :math:`\mathbf{C}` is a matrix whose entires are the counts of the minor allele for each individual for each locus

Hence :math:`\mathbf{C}(i, j)` is how many copies of the minor allele, :math:`a`
an individual, :math:`i`, has at locus :math:`j`. For example:

.. code-block:: python
   :caption: Example of count matrix

   >>> sim.stat(example_pop, alleleFreq=sim.ALL_AVAIL)
   >>> allele_states = analyze.gather_allele_data(example_pop)
   >>> minor_alleles = allele_states[:, 3]  # column corresponding to minor alleles
   >>> segregating_minor_alleles = minor_alleles[segregating_loci]
   >>> gwas = analyze.GWAS(example_pop, segregating_loci, 'example')
   >>> count_matrix = gwas.calculate_count_matrix(segregating_minor_alleles, segregating_loci)
   >>> count_matrix.shape
   (105, 42837)
   >>> print(count_matrix)
   [[1 1 1 ..., 1 1 1]
    [0 0 0 ..., 1 0 0]
    [1 0 0 ..., 1 0 0]
    ...,
    [0 1 0 ..., 2 1 2]
    [0 2 0 ..., 1 0 0]
    [0 0 0 ..., 1 1 0]]

Calculate the mean of each column: :math:`\mathbf{C_{\mu}}(j)`:

.. math::

   \mathbf{C}_{\mu}(j) = \frac{\sum_{i=1}^{m}\mathbf{C}(i, j)}{m}

We "correct" each entry of :math:`\mathbf{C}(i, j)` by the following process.
Let :math:`p(j)` be equal to :math:`\frac{\mathbf{C_{\mu}}(j)}{2}`. We obtain
the matrix :math:`\mathbf{M}` after the "correction" process.

.. math::

   \mathbf{M}(i,j) = \frac{\mathbf{C}(i, j) - \mathbf{C}_{\mu}(j)}{\sqrt{p(j)(1-p(j))}}

We are interested in the principal components of the matrix:

.. math::

   \mathbf{X} = (\frac{1}{n})\mathbf{M}\mathbf{M^T}

Hence we perform the eigenvalue decomposition of \mathbf{X}. We will use the
principal components as co-variates for TASSEL's mixed linear model. But before
that let us see if our results agree. Does this population seem to be descended
from six subpopulations?

.. image:: /images/PC1xPC2.svg
   :align: center

There appear to be six subgroups in the plot. Which is exactly what we obtained
from a different method. For the time being we will use the first two principal
components; however, the Patterson paper describes a test statistic to test
the significance of each principal component.

.. code-block:: python
   :caption: Population structure calculation

   >>> column_means = np.apply_along_axis(np.mean, axis=0, arr=count_matrix)
   >>> shifted = np.array(
   ...    [count_matrix[:, i] - column_means[i] for i in range(42837)]).T
   >>> P = column_means/2
   >>> scale = np.sqrt(P*(1-P))
   >>> M = np.matrix(
   ...     np.array([shifted[:, i] / scale[i] for i in range(42837)]).T)
   >>> X = (1/42837)*(M * M.T)
   >>> eigendata = linalg.eig(X)
   >>> eigenvalues = np.array(eigendata[0], dtype=np.float)
   >>> eigenvectors = np.array(eigendata[1], dtype=np.float)

All of the above calculations are performed by :func:`pop_struct_eigendecomp`.
The file for the population structure covariates is written by
:func:`population_structure_formatter`

.. code-block:: python
   :caption: Population structure file

   >>> gwas.population_structure_formatter(eigenvalues, eigenvectors,
   ...    number_of_pcs=2, pop_struct_file_name='example_structure.txt')

.. _formatting_trait_data:

Trait
=====

Our simulated trait data is given a header and output as a text file.
TASSEL uses html style tags <tag> in the header to label the input of each
file.

.. code-block:: python
   :caption: Functions for handling trait data

   >>> heritability = 0.7
   >>> operators.calculate_g(example_pop)
   >>> operators.calculate_error_variance(example_pop, heritability)
   >>> operators.calculate_p(example_pop)

.. _tassel_config_file:

TASSEL Config File
==================

The final component is a ``xml`` file which specifies the protocol for
TASSEL to run. The config file can be in terms of relative or absolute paths
for its input. All TASSEL options can be specified in the config file.



.. [Patterson2006] Patterson, N, Price, A, Reich, D. (2006). Population Structure and Eigenanalysis. PLOS Genetics, 2(12). doi:10:1371/journal.pgen.0020190

.. [VanRaden2008] VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. Journal of Dairy Science, 91(11), 4414–23. doi:10.3168/jds.2007-0980
