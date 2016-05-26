.. _parse:

=====
Parse
=====

A small module of miscellaneous code which does not fit in anywhere else.
Unsure if I am keeping it or merging it with another module.


.. py:function:: af_from_hapmap(hapmap_file_name, sep='\t')

   :parameter str hapmap_file_name: File name of simulated hapmap results.

   Determines the allele frequencies from a hapmap (.hmp or .txt) file. Allows
   the user to calculate the allele frequencies during the simulation or after
   the simulation.

   At present there are a many unnecessary columns in the hapmap file. Since the
   return value is a pandas.DataFrame we can just drop superfluous columns

   .. code-block:: python

