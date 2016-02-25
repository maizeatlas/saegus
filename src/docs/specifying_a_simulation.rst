==========================
Parameters of a Simulation
==========================


Notes on Standard Simulation Parameters
=======================================


   I am writing this document to record the method i used when developing the
   typical parameters I use in a simulation.


   Recombination Rates

   Alleles
   ~~~~~~~

   The variable ``alleles`` refers to the alleles which are present in the NAM
   prefounders at each locus. The alleles present at each locus is conveniently
   given by the ``alleles`` column in ``hapmap3.txt``. For the sake of clarity I
   put the allele information into a separate file called ``alleles_by_locus.txt``
   which looks like:


   +---------+-----------------+-----------------+
   |  locus  | first(allele)   | second(allele)  |
   +=========+=================+=================+
   |    0    |       A         |      G          |
   +---------+-----------------+-----------------+
   |    1    |       G         |      A          |
   +---------+-----------------+-----------------+
   |    2    |       G         |      A          |
   +---------+-----------------+-----------------+
   |    3    |       C         |      T          |
   +---------+-----------------+-----------------+
   |    4    |       C         |      G          |
   +---------+-----------------+-----------------+
   |    5    |     ``+``       |    ``-``        |
   +---------+-----------------+-----------------+
   |    6    |       G         |      A          |
   +---------+-----------------+-----------------+


   Where `-` means deletion and `+` means insertion. For simuPOP each allele
   must be represented by an integer. So I mapped each allele onto an integer
   starting with 0 according to alphabetical order. The mapping is specified in
   a Python ``dict`` object. The following code snippet shows the code I used to
   map each nucleotide abbreviation onto an integer.


   .. code:: python

      snp_to_integer = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-':4, '+':5}
      alleles = {i: (snp_to_integer[hapmap.ix[i, 'alleles'][0]], snp_to_integer[
                        hapmap.ix[i, 'alleles'][-1]])
                        for i in range(len(hapmap))}

      print(alleles)

       {0: [0, 2],
        1: [2, 0],
        2: [2, 0],
        3: [1, 3],
        4: [1, 2],
        5: [5, 4],
        6: [2, 0],

        ...

        }

   Which would convert the above table into:

   +---------+-----------------+-----------------+
   |  locus  | first(allele)   | second(allele)  |
   +=========+=================+=================+
   |    0    |       0         |      2          |
   +---------+-----------------+-----------------+
   |    1    |       2         |      0          |
   +---------+-----------------+-----------------+
   |    2    |       2         |      0          |
   +---------+-----------------+-----------------+
   |    3    |       1         |      0          |
   +---------+-----------------+-----------------+
   |    4    |       1         |      2          |
   +---------+-----------------+-----------------+
   |    5    |       5         |      4          |
   +---------+-----------------+-----------------+
   |    6    |       2         |      0          |
   +---------+-----------------+-----------------+



   Recombination Rates
   ~~~~~~~~~~~~~~~~~~~

   The variable ``recombination_rates`` or ``recom_rates`` is a list of floats which specify the probability of a
   recombination at the locus immediately *after*. In other words the probability of recombination at the region between
   two loci.


   **Recombination Rates**

   +---------+--------+--------+----------+
   |  locus  |  rate  | after  |  before  |
   +=========+========+========+==========+
   |    0    |  0.00  |   NA   |    NA    |
   +---------+-----------------+----------+
   |    1    |  0.01  |   1    |    2     |
   +---------+-----------------+----------+
   |    2    |  0.00  |   NA   |    NA    |
   +---------+-----------------+----------+
   |    3    |  0.01  |    3   |    4     |
   +---------+-----------------+----------+
   |    4    |  0.00  |   NA   |    NA    |
   +---------+-----------------+----------+
   |    5    |  0.00  |   NA   |    NA    |
   +---------+-----------------+----------+
   |    6    |  0.01  |   6    |    7     |
   +---------+-----------------+----------+


For example in the above table at locus 0 the rate is 0.00. Because there is no recombination after that locus the after
and before columns are `NA`. At locus 1: the rate is 0.01 which means that recombination occurs after locus 1 but before
locus 2 with probability 0.01. So on and so forth for the other loci.





Universal Parameters
++++++++++++++++++++

   I created a .yaml file to hold the parameters which I re-use over and over again. Instead of recreating the
   parameters for every single run of the simulator I opted to write the constantly re-used parameters to a human
   readable format.


#chr_cM_positions = {}
#for i in range(1, 11):
#    chr_cM_positions[i] = []

#for idx in range(len(genetic_map)):
#    chrome = str(int())
#    chr_cM_positions[int(genetic_map.iloc[idx]['chr'])].append(genetic_map.iloc[idx]['cM_pos'])


#cM_positions = []
#for i in range(1, 11):
#    cM_positions.append(chr_cM_positions[i])

#integral_valued_loci = []
#relative_integral_valued_loci = {}
#for idx in range(len(genetic_map)):
#    if str(genetic_map.iloc[idx]['cM_pos'])[-2:] == '.0':
#        integral_valued_loci.append(idx)
#        relative_integral_valued_loci[idx] = (genetic_map.iloc[idx]['chr'], genetic_map.iloc[idx]['cM_pos'])


#alleles = {i: [snp_to_integer[hapmap.ix[i, 'alleles'][0]],
#               snp_to_integer[hapmap.ix[i, 'alleles'][-1]]] for i in
#          range(len(hapmap))}


#recombination_rates = []
#for chromosome in cM_positions:
#    for cM in chromosome:
#        if str(cM)[-2:] == '.6':
#            recombination_rates.append(0.01)
#        else:
#            recombination_rates.append(0.0)


#flat_cM_positions = []
#for cMs in cM_positions:
#    flat_cM_positions.extend(cMs)


#genetic_map_parameters = dict(
#    cM_positions=cM_positions, chr_cM_positions=chr_cM_positions, integral_valued_loci=integral_valued_loci,
#    relative_integral_valued_loci=relative_integral_valued_loci, alleles=alleles,
#    recombination_rates=recombination_rates
#)

.. code:: python

universal_parameters = dict(generations_of_selection=10,
                            generations_of_random_mating=3,
                            number_of_replicates=2,
                            operating_population_size=200,
                            proportion_of_individuals_saved=0.05,
                            overshoot_as_proportion=0.50,
                            individuals_per_breeding_subpop=5, heritability=0.7,
                            prefounder_file_name='nam_prefounders.pop',
                            founders=[[3, 18], [2, 13], [7, 14], [1, 19],
                                      [14, 17], [1, 20], [17, 21], [9, 22]],
                            allele_names=['A', 'C', 'T', 'G', 'D', 'I'],
                            snp_to_integer={'A': 0, 'C': 1, 'G': 2, 'T': 3,
                                            '-': 4, '+': 5},
                            integer_to_snp={0: 'A', 1: 'C', 2: 'G', 3: 'T',
                                            4: '-', 5: '+'},
                            number_of_qtl=10,
                            allele_effect_parameters=[1], pos_column=pos_column,
                            locus_names=locus_names,
                            )









