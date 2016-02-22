NAM Prefounders
===============


Features of NAM Prefounders
---------------------------

   The twenty-six lines which were used to make the maize nested association
   mapping population are used as *prefounders* in ``wgs`` to construct other
   populations. The genetic map of the NAM *prefounders* used in ``wgs``
   gives a SNP every 0.2 centiMorgans (cM). The SNPs are 'A', 'C', 'G', 'T',
   '--', '+' where "--", "+" stand for deletion, insertion respectively. The
   raw genetic map gives genotypes at 7368 loci across the ten chromosomes.

Triplet Scheme
~~~~~~~~~~~~~~

   In ``wgs`` we are using a special *triplet* structure of a locus.
   The newest version of ``wgs`` uses all 7368 loci for the genotypes
   of the prefounders; however, we only use the loci which are at integral
   valued positions in the genetic map are viable QTL (quantitative trait
   locus). There are 1478 loci which have integral-valued cM positions. QTL
   are designated by choosing *n* from these 1478 loci. Then the loci
   immediately flanking a QTL plus/minus 0.2 cM form a ``haplotype`` with the
   central locus. A ``haplotype`` in this case is a *triplet* of
   non-recombining loci i.e. there is no recombination between the triplet of
   alleles.

   However, since ``wgs``-0.1.1 there are extra loci in addition to the triplet
   loci. ``wgs``-0.1.1 accommodates this by extending the area of zero
   recombination to include loci which are plus/minus 0.4 cM away from the
   integral-valued locus. This allows ``wgs`` to make greater use of the native
   simuPOP operators such as Stat.haploFreq: a statistical operator which
   allows the user to calculate the frequency of specific haplotypes.

Quantitative Trait Loci
~~~~~~~~~~~~~~~~~~~~~~~

   QTL are chosen at random from the 1478 loci which have integral-valued cM
   positions in the genetic map. The triplet scheme is effected by
   adding the loci immediately flanking the random sample of ``n`` of the
   1478 loci. So for example:
   +-------------+----------------+
   |    Locus    |    Triplet     |
   +=============+================+
   |     100     | 99, 100, 101   |
   +-------------+----------------+
   |     278     | 277, 278, 279  |
   +-------------+----------------+
   |     900     | 899, 900, 901  |
   +-------------+----------------+
