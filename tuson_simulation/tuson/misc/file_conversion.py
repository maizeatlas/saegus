#! /home/jjdoc/anaconda3/bin/python3.4
import numpy as np


batchfiles = ['batch_b'+str(batch)+'_af.txt' for batch in range(43, 48)]

for bfile in batchfiles:
    batch = np.loadtxt(bfile)
    np.savetxt(bfile, batch, fmt='%.3f', delimiter='\t')
