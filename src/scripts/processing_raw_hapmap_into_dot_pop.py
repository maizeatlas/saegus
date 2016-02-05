import pandas as pd
import collections as col
from wgs import parser
import random
import numpy as np

random.seed(1337)

'''
This is the exact sequence of commands which I used to process hapmap3.txt
into separate objects which I used to create nam_prefounders.pop
'''

hapmap = pd.read_csv('hapmap3.txt', sep='\t', index_col=0)
col_numbers = [0, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20]
droppable_columns = []
for number in col_numbers:
    droppable_columns.append(hapmap.columns[number])
hapmap.drop(droppable_columns, axis=1, inplace=True)
founder_indices = [i for i in range(3, 29)]
new_founder_names = {hapmap.columns[founder_index]: hapmap.columns[founder_index][:-3]
                    for founder_index in founder_indices}
hapmap.rename(columns=new_founder_names, inplace=True)
hapmap.to_csv('clean_hapmap.txt', index=None)
hapmap = pd.read_csv('clean_hapmap.txt')
genetic_map = hapmap.ix[:, :'cM_pos']
genetic_map.to_csv('nam_prefounders_genetic_map.txt', sep='\t', index=None)
genetic_map = pd.read_csv('nam_prefounders_genetic_map.txt', index_col=None,
                         sep='\t')
chr_cm_positions = col.OrderedDict()
for i in range(1, 11):
    chr_cm_positions[i] = []
for idx in range(len(genetic_map)):
    chr_cm_positions[int(genetic_map.iloc[idx]['chr'])].append(
    float(genetic_map.iloc[idx]['cM_pos']))
cM_positions = []
for k, v in chr_cm_positions.items():
    cM_positions.append(v)
cM_positions = []
for k, v in chr_cm_positions.items():
    cM_positions.append(v)
snp_to_integer = {'A':0, 'C':1, 'G':2, 'T':3, '-':4, '+':5}
all_genotypes_handle = [
    [snp_to_integer[genotypes[name][idx][0]] for idx in range(len(hapmap))] +
    [snp_to_integer[genotypes[name][idx][1]] for idx in range(len(hapmap))]
            for name in names]
integral_valued_loci = []
relative_integral_valued_loci = {}
for idx in range(len(genetic_map)):
    if str(genetic_map.iloc[idx]['cM_pos'])[-2:] == '.0':
        integral_valued_loci.append(idx)
        relative_integral_valued_loci[idx] = (genetic_map.iloc[idx]['chr'], genetic_map.iloc[idx]['cM_pos'])

        
triplet_qtl = []
for locus in qtl:
    triplet_qtl.append(locus-1)
    triplet_qtl.append(locus)
    triplet_qtl.append(locus+1)
triplet_qtl = sorted(triplet_qtl)
alleles = {i: (snp_to_integer[hapmap.ix[i, 'alleles'][0]], snp_to_integer[hapmap.ix[i, 'alleles'][-1]]) for i in
          range(len(hapmap))}

allele_effects = {}
for tqtl in triplet_qtl:
    for allele in alleles[tqtl]:
        allele_effects[tqtl, allele] = random.expovariate(1)

recombination_rates = []
for chromosome in cM_positions:
    for cM in chromosome:
        if str(cM)[-2:] == '.6':
            recombination_rates.append(0.01)
        else:
            recombination_rates.append(0.0)