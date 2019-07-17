
.. _magic_notes:

=====================================================
MAGIC: Multi-parent Advanced Generation Inter-crosses
=====================================================

Notes on MAGIC
==============

MAGIC begins with ``founders`` arranged into pairs. Pairs of parents
are crossed with each other to make hybrid offspring. The hybrid offspring
are crossed with a different group of hybrid offspring to make
double-hybrid offspring. This process continues until only a single
sub-population of individuals remains. Given the ``founders`` MAGIC
merges successive generations. The process is demonstrated visually with
a graph.

.. code-block:: python

   >>> founders = [[1, 2], [3, 4], [5, 6], [7, 8]]

.. figure:: reformed_graph.png
   :scale: 50%

   Successive rounds of mating merges together the genomes of distinct founders.



The result of these three rounds of crossing is a population with high
genetic diversity. High genetic diversity stemming from the combination of
many different parental lines greatly enhances mapping resolution.


Predicting Number of Rounds of Mating
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given :math:`n` founders to start with:

If :math:`n` is the :math:`m` a power of two, :math:`m`, then there are *m* rounds of mating.
That is if:

.. math::

   2^m = n

The number of rounds of mating is the largest power of :math:`2` which can
be subtracted from :math:`n` and then one more generation to combine the final pair
of individuals.
