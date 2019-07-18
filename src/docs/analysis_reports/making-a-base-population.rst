========================
Making a Base Population
========================


Overview
--------

saegus is written with a specific strategy for simulating whole populations.
Typically the user starts with a relatively small number of individuals
called 'founders' or 'pre-founders'. In the case of the *maize NAM* population
26 lines were used to make all 5,000 lines by combinatorial crossing.
In the *Tuson* population we used a sample of 105 individuals to create
a population of 10,000. In both cases we have a file which contains the
genotype information of all founder individuals and a genetic map which
contains the number of loci, their distribution on the chromosomes and
recombination rates of each sequential pair of loci.

Parsing Genotype Data
---------------------

Imports to handle parsing:

.. code:: python

    import numpy as np
    import pandas as pd
    import collections as col

Genotype matrix and genetic map filenames:

.. code:: python

    raw_genotype_matrix_filename = 'raw_genotype_matrix.txt'
    genetic_map = 'genetic_map.txt'


Genotype File Formats
~~~~~~~~~~~~~~~~~~~~~

Here are some examples of file formats for genotype data that I have encountered.

.. table:: Example genotype file

+----------+--------+--------+
| label    | locus1 | locus2 |
+==========+========+========+
| C0_164_1 | 2/1    | 3/2    |
+----------+--------+--------+
| C0_164_2 | 2/1    | NA     |
+----------+--------+--------+
| C0_164_3 | 2/1    | 3/2    |
+----------+--------+--------+

.. table:: Another common format

+------------+----------------+----------------+
|  position  | individual_one | individual_two |
+============+================+================+
|   50486    |      0/0       |       0/1      |
+------------+----------------+----------------+
|   59185    |      0/2       |       ./.      |
+------------+----------------+----------------+
|   59188    |      0/1       |       0/0      |
+------------+----------------+----------------+

The *hapmap* format is a text file with certain information. Most of it is
extraneous for the purposes of creating a base population. However, a key
aspect of the *hapmap* format is that it provides the genetic map and the
genotypes.

.. table:: Hapmap format

+-----+--------+------+---------+----------------+----------------+
| chr | marker |  cM  | alleles | individual_one | individual_two |
+=====+========+======+=========+================+================+
|  1  |   m1   |  1.0 |   A/G   |      AA        |      AG        |
+-----+--------+------+---------+----------------+----------------+
|  1  |   m2   |  1.2 |   A/C   |      CC        |      NA        |
+-----+--------+------+---------+----------------+----------------+
|  2  |  m281  |  1.0 |   G/C   |      CG        |      GG        |
+-----+--------+------+---------+----------------+----------------+


Parsing Genotype File
~~~~~~~~~~~~~~~~~~~~~

We need to get the genotype out of the file and into a Python object.
The pandas_ package makes parsing the file very simple. For example
if we have a tab-delimited file named ``genotypes.txt`` we can assign
the file contents to ``genotype_data``.

.. _pandas: http://pandas.pydata.org/


.. code:: python

    >>> genotype_data = pd.read_csv(genotype_matrix_filename, sep='\t')
    >>> genotype_data

.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th><th>POS</th><th>ind_one</th><th>ind_two</th><th>ind_three</th><th>ind_four</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th><td>50486</td><td>0/0</td><td>./.</td><td>0/0</td><td>0/0</td>
        </tr>
        <tr>
          <th>1</th><td>59185</td><td>0/0</td><td>0/0</td><td>0/0</td><td>0/0</td>
        </tr>
        <tr>
          <th>2</th><td>59188</td><td>0/0</td><td>0/2</td><td>0/1</td><td>0/0</td>
        </tr>
        <tr>
          <th>3</th><td>59189</td><td>0/0</td><td>./.</td><td>0/1</td><td>0/0</td>
        </tr>
      </tbody>
    </table>
    </div>

The "./." in ``genotype_data`` represents missing data. We cannot generate a population
if some individuals have missing data; however, we do not want
to throw away individuals with missing data. All of the genotype
files I have worked with heretofore had missing data. We take
a simple approach to fill-in missing data.


Missing Genotype Data
~~~~~~~~~~~~~~~~~~~~~

Albeit there might be several ways of handling missing genotype data I
use the genotype frequencies at each locus to construct probability mass functions.
The random variable are the genotypes at that locus and the probabilities are simply
the frequencies of a genotype at that locus. Thus whenever an
individual is missing data we replace the missing data with a draw from the pmf
of the genotype frequencies for that locus. For example suppose that the
genotype frequencies for locus2 are: 3/2: 0.65, 2/2: 0.10, 3/3: 0.25
Assuming that we randomly draw 2/2 we would update our genotype matrix as:

Genotype matrix with missing data replaced by random draw:

+--------+----------+----------+----------+
| label  | C0_164_1 | C0_164_2 | C0_164_3 |
+========+==========+==========+==========+
| locus1 | 2/1      | 2/1      | 2/1      |
+--------+----------+----------+----------+
| locus2 | 3/2      | 2/2      | 3/2      |
+--------+----------+----------+----------+

We repeat the same procedure for every locus until all all loci are assigned an
appropriate genotype.


Genetic Map
-----------

A genetic map gives the positions of genes (in our case single nucleotide
polymorphisms) in terms of recombination between homologous chromosomes. Our
particular map has units of centiMorgan (cM). A genetic distance of 1 cM implies
that there is 1% recombination between a pair of loci. The Tuson genetic map
uses 44,445 markers and is structured like this:

+--------+-----+-------------+
| locus  | chr | cM          |
+========+=====+=============+
| locus1 | 1   | 2.279635824 |
+--------+-----+-------------+
| locus2 | 1   | 2.290091023 |
+--------+-----+-------------+
| locus3 | 1   | 2.290716944 |
+--------+-----+-------------+

simuPOP does not recognize any unit of measurement of recombination. simuPOP
requires the *probability* of recombination between sequential loci. Therefore,
we divide the absolute difference between sequential loci by 100 to obtain
the recombination rates which are supplied as input to the simuPOP Recombinator function.
