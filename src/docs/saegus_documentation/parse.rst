.. _parse:

=====
Parse
=====

A small module of miscellaneous code which does not fit in anywhere else.
Unsure if I am keeping it or merging it with another module.

*Update*
========

I have recently found a good use of parse. Aside from containing a very small
number of functions used to handle raw genotype data I will use it to derive
data from simulated output i.e. allele frequencies from genotype files.


.. py:function:: af_from_hapmap(hapmap_file_name, sep='\t')

   :parameter str hapmap_file_name: File name of simulated hapmap results.

   Determines the allele frequencies from a hapmap (.hmp or .txt) file. Allows
   the user to calculate the allele frequencies during the simulation or after
   the simulation.

   At present there are a many unnecessary columns in the hapmap file. Since the
   return value is a pandas.DataFrame we can just drop superfluous columns. We
   convert it into a numpy array for quicker calculation of frequencies.

   .. code-block:: python

      >>> af_from_hmp = pd.read_csv('R0_100_infinite_simulated_hapmap.txt', sep='\t')
      >>> af_from_hmp.drop(['rs', 'alleles', 'chrom', 'pos', 'strand',
      ...                  'assembly', 'center', 'protLSID', 'assayLSID',
      ...                  'panelLSID', 'QCcode'], axis=1, inplace=True)
      >>> genotypes = np.array(af_from_hmp)
      >>> genotypes
      array([['CC', 'CT', 'CC', ..., 'CC', 'CT', 'CC'],
             ['TT', 'TT', 'TT', ..., 'TT', 'CT', 'TT'],
             ['AG', 'GG', 'GG', ..., 'GG', 'GG', 'AG'],
             ...,
             ['GG', 'AA', 'AA', ..., 'AA', 'AG', 'AA'],
             ['AA', 'CC', 'CC', ..., 'AC', 'AC', 'CC'],
             ['CT', 'TT', 'TT', ..., 'TT', 'CT', 'TT']], dtype=object)