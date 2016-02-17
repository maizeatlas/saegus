
# coding: utf-8

# In[1]:





import simuOpt
simuOpt.setOptions(alleleType='lineage', optimized=True, quiet=True)
import simuPOP as sim
import pandas as pd
import collections as col
from wgs import breed, operators, selection, helpers, parser, parameterizer, selection
import random
import numpy as np
random.seed(1337)


hapmap = pd.read_csv('clean_hapmap.txt')
genetic_map = hapmap.ix[:, :'cM_pos']
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


snp_to_integer = {'A':0, 'C':1, 'G':2, 'T':3, '-':4, '+':5}

integral_valued_loci = []
relative_integral_valued_loci = {}
for idx in range(len(genetic_map)):
    if str(genetic_map.iloc[idx]['cM_pos'])[-2:] == '.0':
        integral_valued_loci.append(idx)
        relative_integral_valued_loci[idx] = (genetic_map.iloc[idx]['chr'], genetic_map.iloc[idx]['cM_pos'])


qtl = sorted(random.sample(integral_valued_loci, 5))

triplet_qtl = []
for locus in qtl:
    triplet_qtl.append(locus-1)
    triplet_qtl.append(locus)
    triplet_qtl.append(locus+1)
triplet_qtl = sorted(triplet_qtl)


alleles = {i: (snp_to_integer[hapmap.ix[i, 'alleles'][0]], 
               snp_to_integer[hapmap.ix[i, 'alleles'][-1]]) for i in
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

allele_names = ['A', 'C', 'T', 'G', 'D', 'I']

flat_cM_positions = []
for cMs in cM_positions:
    flat_cM_positions.extend(cMs)


nam = sim.loadPopulation('nam_prefounders.pop')
nam.setSubPopName('prefounders', 0)
sample_sizes = {i: 100 for i in range(0, 11, 2)}

s = selection.Truncation(10,
                       3,
                       2000,
                       0.05,
                       0.50,
                       5,
                       0.7,
                       sample_sizes,
                       1)

nam.dvars().qtl = qtl
nam.dvars().triplet_qtl = triplet_qtl
nam.dvars().allele_effects = allele_effects
sim.tagID(nam, reset=True)

founders = [(1, 5), (7, 8)]
replicated_nam = sim.Simulator(nam, rep=2)
pop = replicated_nam.extract(0)
meta = replicated_nam.extract(0)
meta.removeSubPops(0)
s.generate_f_one(pop, recombination_rates, founders)
s.generate_f_two(pop, recombination_rates)
s.interim_random_mating(pop, recombination_rates)
s.mate_and_merge(pop, recombination_rates)
s.recurrent_truncation_selection(pop, meta, triplet_qtl, allele_effects,
                                recombination_rates)
                                


# In[12]:

pop.popSize()


# In[13]:

pop.save('testpop_recurrent_selection.pop')


# In[3]:

meta.popSize()


# In[4]:

len(meta.indInfo('ind_id'))


# In[5]:

meta.indInfo('generation')


# In[ ]:

pop.popSize()


# In[ ]:

pop.infoFields()


# In[ ]:

pop.__cmp__(meta)


# In[ ]:

pop == meta


# In[ ]:

s.generate_f_two(pop, recombination_rates)


# In[ ]:




# In[ ]:

replicated_nam.numRep()


# In[ ]:




# In[ ]:

parents = [(0,1), (10, 11)]


# In[ ]:

generator = (pair for pair in parents)


# In[ ]:

next(generator)


# In[ ]:

def parent_generator(pair_of_parents):
    for pair in pair_of_parents


# In[ ]:

pop.indInfo('generation')


# In[ ]:

pop.popSize()


# In[ ]:

pop.indInfo('ind_id')


# In[ ]:

s.generate_f_two(pop, recombination_rates)


# In[ ]:

pairsOfParents = [(0, 1), (10, 11)]


# In[ ]:

female, male = pairsOfParents[0]


# In[ ]:

male


# In[ ]:

len(pop.indInfo('ind_id'))


# In[ ]:

s.mate_and_merge(pop, recombination_rates)


# In[ ]:

pop.indInfo('ind_id')


# In[ ]:

len(pairsOfParents)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

s.recurrent_truncation_selection(pop, meta, triplet_qtl, allele_effects,
                                recombination_rates)


# In[7]:

allele_effects


# In[8]:

triplet_qtl


# In[9]:

sim.stat(pop, alleleFreq=triplet_qtl)


# In[10]:

pop.dvars().alleleFreq


# In[11]:

allele_effects[3750, 3]


# In[6]:

pop.indInfo('g')


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

import simuOpt
simuOpt.setOptions(alleleType='short', optimized=True, numThreads=4, quiet=True)
import simuPOP as sim
import pandas as pd
import collections as col
from wgs import breed, operators, selection, helpers, parser, parameterizer, selection
import random
import numpy as np
random.seed(1337)


# In[ ]:

hapmap = pd.read_csv('clean_hapmap.txt')


# In[ ]:

genetic_map = hapmap.ix[:, :'cM_pos']


# In[ ]:

genetic_map = pd.read_csv('nam_prefounders_genetic_map.txt', index_col=None,
                         sep='\t')


# In[ ]:

len(genetic_map)


# In[ ]:

chr_cm_positions = col.OrderedDict()
for i in range(1, 11):
    chr_cm_positions[i] = []


# In[ ]:

for idx in range(len(genetic_map)):
    chr_cm_positions[int(genetic_map.iloc[idx]['chr'])].append(
    float(genetic_map.iloc[idx]['cM_pos']))


# In[ ]:

cM_positions = []
for k, v in chr_cm_positions.items():
    cM_positions.append(v)


# In[ ]:

snp_to_integer = {'A':0, 'C':1, 'G':2, 'T':3, '-':4, '+':5}


# In[ ]:

integral_valued_loci = []
relative_integral_valued_loci = {}
for idx in range(len(genetic_map)):
    if str(genetic_map.iloc[idx]['cM_pos'])[-2:] == '.0':
        integral_valued_loci.append(idx)
        relative_integral_valued_loci[idx] = (genetic_map.iloc[idx]['chr'], genetic_map.iloc[idx]['cM_pos'])


# In[ ]:

qtl = sorted(random.sample(integral_valued_loci, 5))


# In[ ]:

triplet_qtl = []
for locus in qtl:
    triplet_qtl.append(locus-1)
    triplet_qtl.append(locus)
    triplet_qtl.append(locus+1)
triplet_qtl = sorted(triplet_qtl)


# In[ ]:

alleles = {i: (snp_to_integer[hapmap.ix[i, 'alleles'][0]], snp_to_integer[hapmap.ix[i, 'alleles'][-1]]) for i in
          range(len(hapmap))}


# In[ ]:

allele_effects = {}
for tqtl in triplet_qtl:
    for allele in alleles[tqtl]:
        allele_effects[tqtl, allele] = random.expovariate(1)


# In[ ]:

recombination_rates = []
for chromosome in cM_positions:
    for cM in chromosome:
        if str(cM)[-2:] == '.6':
            recombination_rates.append(0.01)
        else:
            recombination_rates.append(0.0)


# In[ ]:




# In[ ]:

allele_names = ['A', 'C', 'T', 'G', 'D', 'I']


# In[ ]:

flat_cM_positions = []
for cMs in cM_positions:
    flat_cM_positions.extend(cMs)


# In[ ]:

nam = sim.loadPopulation('nam_prefounders.pop')
sample_sizes = {i: 100 for i in range(0, 11, 2)}


# In[ ]:

sample_sizes


# In[ ]:

s = selection.Truncation(10,
                       3,
                       2000,
                       0.05,
                       0.50,
                       5,
                       0.7,
                       sample_sizes,
                       1)


# In[ ]:

s.breeding_parameters


# In[ ]:




# In[ ]:

nam.dvars().qtl = qtl
nam.dvars().triplet_qtl = triplet_qtl
nam.dvars().allele_effects = allele_effects
sim.tagID(nam, reset=True)


# In[ ]:

#nam.save('nam_prefounders.pop')
#nam = sim.loadPopulation('nam_prefounders.pop')
print(nam.vars().keys())


# In[ ]:

founders = [1, 5, 7, 8]
pop = nam.clone()
meta = nam.clone()
#meta.removeIndividuals(list(range(meta.popSize())))
s.arbitrary_ordering_of_founders(nam, pop, founders)
s.generate_f_one(pop, recombination_rates)


# In[ ]:

s.generate_f_two(pop, recombination_rates)


# In[ ]:

s.interim_random_mating(pop, triplet_qtl, allele_effects, recombination_rates)


# In[ ]:

s.mate_and_merge(pop, recombination_rates)


# In[ ]:

s.recurrent_truncation_selection(pop, meta, triplet_qtl, allele_effects,
                                recombination_rates)


# In[ ]:




# In[ ]:




# In[ ]:



