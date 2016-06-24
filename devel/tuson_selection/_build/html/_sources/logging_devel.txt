.. _tuson_with_selection:

########################################
Tuson Revisited: Simulation of Selection
########################################

Many months ago we simulated genetic drift in the Tuson population. For comparison
we will now simulate truncation selection which mirrors the experimental protocol.
In other simulations of selection we have randomly assigned QTL and allele effects
were determined by random draws from some distribution (i.e. exponential). In
this case we will be using the experimentally estimated allele effects of every
single locus.

Initial Setup and Data Collection
=================================

We will use G:sub:`0` as founder population. In G:sub:`0` some sites are fixed
which are segregating in later generations. In other words the G:sub:`0` sample
failed to capture all of the segregating loci. We assume that this implies the
sites which are fixed in the sample have low frequency alleles.

Conversion of Fixed Sites to Segregating Sites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We determine which sites are fixed in the sample and then randomly assign a
single individual an allele which is **not** present at the fixed locus. This
converts a fixed site into a segregating site of low frequency.

.. code-block:: py
   :caption: Convert fixed sites to low frequency segregating sites

   >>> sim.stat(tuson, numOfSegSites=sim.ALL_AVAIL, vars=['numOfSegSites', 'segSites', 'numOfFixedSites', 'fixedSites'])
   >>> parameters.randomly_convert_fixed_sites(tuson, tuson.dvars().fixedSites)
   >>> tuson.dvars().numOfFixedSites
   0
   >>> tuson.dvars().numOfSegSites
   44445

Storage of Founder Data
~~~~~~~~~~~~~~~~~~~~~~~

Since each conversion of fixed sites to low frequency segregating sites is a
different random outcome I store each realization before carrying the simulation
forward. We store:

+ Alleles at each locus
+ Allele frequencies
+ Homozygote frequencies
+ Heterozygote frequencies

Examples of
