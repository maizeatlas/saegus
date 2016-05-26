Genotype Data
=============


Handling missing genotype data
------------------------------

In most cases genotype-by-sequencing will have some missing data. This is a
tutorial for handling missing genotype data.

.. code:: python

    import numpy as np
    import pandas as pd
    import h5py as hp
    import collections as col
    import scipy.stats as st


In this instance Randy wanted to include three more loci so I had to re-create a base population file to accommodate
the new genome structure. The new genotype matrix is in a file called 'raw_genotype_matrix.txt' and the corresponding
genetic map in 'genetic_map.txt'.


.. code:: python

    genotype_matrix_filename = 'raw_genotype_matrix.txt'
    genetic_map_filename = 'genetic_map.txt'



.. code:: python

    def parse_genotype_matrix(genotype_matrix_filename: str):
        genotype_matrix = pd.read_csv(genotype_matrix_filename, sep='\t', index_col=0, low_memory=False)
        genotype_matrix = genotype_matrix.drop('popdata', axis=1)
        return genotype_matrix

.. code:: python

    def parse_recombination_rates(genetic_map_filename):
        """
        Returns a list of crossover probabilities from a genetic map measured in centimorgans.
        """
        genetic_map = np.genfromtxt(genetic_map_filename, delimiter='\t', usecols=[1, 3])
        recombination_rates = col.OrderedDict()
        for i in range(1, len(genetic_map), 1):
            if genetic_map[i-1][0] == genetic_map[i][0]:
                recombination_rates[i] = np.divide(np.abs(genetic_map[i][1] - genetic_map[i-1][1]), 100)
            elif genetic_map[i-1][0] != genetic_map[i][0]:
                recombination_rates[i] = 0.0
        recombination_rates[len(genetic_map)] = 0.0
        return list(recombination_rates.values())

.. code:: python

    def end_pt_finder(genetic_map_filename, number_of_chromosomes=10):
        genetic_map = pd.read_csv(genetic_map_filename, sep='\t', index_col=None)
        genetic_map.drop(['locus', 'agpv2', 'namZmPRDA', 'namZmPRDS'], axis=1, inplace=True)
        genetic_map = np.array(genetic_map)
        end_pts = [0]
        for i in range(1, len(genetic_map)):
            if genetic_map[i-1][0] != genetic_map[i][0]:
                end_pts.append(i)
            else:
                pass
        end_pts.append(len(genetic_map))
        chrom_lengths = [(end_pts[i] - end_pts[i-1]) for i in range(1, number_of_chromosomes+1)]
        del end_pts[0]
        return end_pts, chrom_lengths

.. code:: python

    endpts, lens = end_pt_finder('genetic_map.txt')

.. code:: python

    gmatrix = parse_genotype_matrix('raw_genotype_matrix.txt')

.. code:: python

    garray = np.array(gmatrix)

.. code:: python

    missing_loci = [i for i in range(garray.shape[1]) if np.nan in geno_cnt[i]]

.. code:: python

    geno_cnt = {i:col.Counter(garray[:, i]) for i in range(garray.shape[1])}

.. code:: python

    def genotype_counts_to_frequencies(genotype_counts):
        """
        Converts a the dictionaries of genotype: count for each locus into their
        frequency equivalents by dropping and missing data and dividing by the adjusted
        total.
        """
        geno_frq = {}
        for mlocus in missing_loci:
            geno_frq[mlocus] = {}
            if np.nan in geno_cnt[mlocus]:
                del geno_cnt[mlocus][np.nan]
            inds_counted = sum(geno_cnt[mlocus].values())
            for genotype, cnt in geno_cnt[mlocus].items():
                geno_frq[mlocus][genotype] = cnt/inds_counted
        return geno_frq


.. code:: python

    def generate_genotype_pmfs(genotype_frequencies):
        """
        For the time being all of the information required to compute a custom
        probability mass function for each locus is stored a dictionary keyed by locus.
        The values are tuples:
        0: genotype: frequency
        1: integer: genotype
        2: density
        3: genotype: integer
        """
        genotype_pmfs = {}
        loci = list(genotype_frequencies.keys())
        for missinglocus in loci:
            genotype_pmfs[missinglocus] = {}
            genotype_pmfs[missinglocus]['states'] = [geno_state for geno_state in genotype_frequencies[missinglocus].keys()]
            genotype_pmfs[missinglocus]['frequencies'] = [frequency for frequency in genotype_frequencies[missinglocus].values()]
            # scipy.stats.rv_discrete only allows integer or float valued random variables.
            # So we have to map an integer to each genotype state value and vice versa.
            genotype_pmfs[missinglocus]['state_to_integer'] = {state:i for i, state 
                                                               in enumerate(genotype_pmfs[missinglocus]['states'])}
            genotype_pmfs[missinglocus]['integer_to_state'] = {i:state for i, state in enumerate(genotype_pmfs[missinglocus]['states'])}
            integer_states = list(genotype_pmfs[missinglocus]['state_to_integer'].values())
            genotype_pmfs[missinglocus]['pmf'] = st.rv_discrete(values=(integer_states, genotype_pmfs[missinglocus]['frequencies']))
        return genotype_pmfs

.. code:: python

    gf = genotype_counts_to_frequencies(geno_cnt)

.. code:: python

    gframe = pd.DataFrame(garray, index=gmatrix.index, columns=gmatrix.columns)

.. code:: python

    gframe.to_csv('genotype_matrix.txt', sep='\t')

