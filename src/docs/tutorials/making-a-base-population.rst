Making a Base Population
========================

saegus is written with a specific strategy for simulating whole populations.
Typically the user starts with a relatively small number of individuals
called 'founders' or 'pre-founders'. In the case of the maize NAM population
26 lines were used to make all 5,000 lines by combinatorial crossing.
In the Tuson population we used a sample of 105 individuals to create
a population of 10,000. In both cases we have a file which contains the
genotype information of all founder individuals and a genetic map which
gives us the crossover probabilities of each sequential pair of loci.

Parsing Genotype Data
---------------------

Imports to handle parsing:

.. code:: python

    import numpy as np
    import pandas as pd
    import collections as col
    from wgs import parser

Genotype matrix and genetic map filenames:

.. code:: python

    raw_genotype_matrix_filename = 'raw_genotype_matrix.txt'
    genetic_map = 'genetic_map.txt'

Genotype Matrix
%%%%%%%%%%%%%%%

wgs assumes that the genotype data will have the format rows=loci and
columns=individuals. For example the generation zero Tuson population has
105 individuals and 44445 markers. It is typical to have genotype file which
must be processed or edited before it can be used. The raw genotype matrix for
the Tuson population looks like this:

Example raw genotype matrix:

+----------+--------+--------+
| label    | locus1 | locus2 |
+==========+========+========+
| C0_164_1 | 2/1    | 3/2    |
+----------+--------+--------+
| C0_164_2 | 2/1    | NA     |
+----------+--------+--------+
| C0_164_3 | 2/1    | 3/2    |
+----------+--------+--------+


Notice that we have two issues:

- The data is transposed i.e. rows=individuals and columns=loci
- We have missing i.e. NA

First we read the genotype matrix in and transpose it to give us the required
row-column structure:

.. code:: python

    gmatrix = pd.read_csv(genotype_matrix_filename, sep='\t', index_col=0, low_memory=False)
    gmatrix = gmatrix.T

We then have:

+--------+----------+----------+----------+
| label  | C0_164_1 | C0_164_2 | C0_164_3 |
+========+==========+==========+==========+
| locus1 | 2/1      | 2/1      | 2/1      |
+--------+----------+----------+----------+
| locus2 | 3/2      | NA       | 3/2      |
+--------+----------+----------+----------+

Now we need to deal with the missing genotype data.

Missing Genotype Data
%%%%%%%%%%%%%%%%%%%%%

Albeit there might be several ways of handling missing genotype data I decided
to use the genotype frequencies of Tuson generation zero sample of 105
individuals as population probability mass functions. Thus whenever an
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
we take divide the absolute difference between sequential loci by 100 to obtain
the recombination rates which are supplied as input to the simuPOP Recombinator.
