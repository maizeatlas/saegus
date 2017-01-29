import simuOpt
simuOpt.setOptions(alleleType='short', optimized=True, quiet=True, numThreads=4)
import simuPOP as sim
import numpy as np
import pandas as pd
from saegus import analyze

input_file_prefix = '/home/vakanas/tassel-5-standalone/output/'

mafrqs = pd.read_csv(input_file_prefix+'epsilon_0_maf_table.txt', sep='\t', index_col=0)
qtad = pd.read_csv(input_file_prefix+'epsilon_0_quant_allele_table.txt', sep='\t', index_col=0)

super_table = analyze.tassel_results_tables(input_file_prefix+'epsilon_0_out_2.txt', input_file_prefix+'epsilon_0_qvalues.txt',
                                            mafrqs, qtad)


super_table.to_csv(input_file_prefix+'epsilon_0_super_table.txt', sep='\t')