Overview of Tuson Simulation
============================


Tuson Population
----------------

Hallauer's Tuson is a landrace population of *Zea mays* that was adapted to a 
temperate environment. The base population of Tuson was created through a single 
generation of open pollination among multiple populations of the landrace Tuson. 

We designed a computer simulation using simuPOP_ and customized scripts in
Python. The computer simulation uses genotype data from a 55k SNP chip and
recombination frequencies from the empirically determined genetic map. The
computer simulation was run for 10,000 replicates under conditions of genetic
drift. The simulated allele frequency data was used as a null distribution
when testing for selection for each locus.

Genotype Data
-------------

   Genotype data was collected from samples of the Tuson population in
   G\ :sub:`0`, G\ :sub:`2`, G\ :sub:`4`, G\ :sub:`6`, G\ :sub:`8`,
   G\ :sub:`10`. Genotype data was collected with a 55k SNP chip. 
   The genotype data of the 105 individuals of G\ :sub:`0` were used as a
   base population from a computer simulation. Moreover, linkage information 
   from the genetic map is used to determine the recombination rate between 
   adjacent pairs of loci.

   We handled missing data by calculating the genotype frequencies at each locus.
   Missing genotypes were replaced with a random draw from the empirical
   genotype probability mass functions.

Population Expansion from G\ :sub:`0` Sample
--------------------------------------------

   The simulation uses the 105 individuals to create full sized *in silico*
   population of 10,000 individuals. The full sized *in silico* population is
   subjected to conditions of genetic drift for ten generations. We estimated
   the population structure of the 105 individuals by estimating the
   contribution of each founder.

   The genomes of the 105 individuals from G\ :sub:`0` reflected the open
   pollination of the founding populations. From the genomes of the sample 
   of the base population we estimated six hypothetical founder populations. 
   Each genome of the 105 individuals is a combination of one or more of 
   the theoretical founders.

Alleles Absent in G\ :sub:`0` Sample
------------------------------------

   There were some loci at fixation in G\ :sub:`0` and not at fixation in at
   least one of G\ :sub:`2` , G\ :sub:`4` , G\ :sub:`6` , G\ :sub:`8`,
   G\ :sub:`10` . Alleles not seen in the G\ :sub:`0` sample but present in a
   future generation sample implies there are low frequency alleles in
   G\ :sub:`0` We accounted for this by changing one allele of a single, randomly
   chosen individual at each fixed locus in the sample of G\ :sub:`0`.

Simulating Genetic Drift
------------------------

   The *in silico* population simulates drift by randomly selecting 400
   individuals as *females* and 400 distinct individuals as
   *males*. A male or female is capable of self-mating or out-crossing with
   the opposite sex. In each generation 10,000 offspring are generated.
   Each mating event produces a single individual. Male and female are
   chosen at random independent of each other. If the *male* parent and
   *female* parent are the same individual selfing occurs.











.. _simuPOP: http://simupop.sourceforge.net/
