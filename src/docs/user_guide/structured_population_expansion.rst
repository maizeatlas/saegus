.. _structured_population_expansion:

####################################################
Expanding A Population According to Structure Matrix
####################################################

This example will show how to take a small population and expand it to a larger
population using a population structure matrix. The population structure matrix
defines the proportion of the genome inherited by an individual from each
sub-population. Each individual is assigned a ``primary`` subpopulation. The
primary sub-population is the sub-population from which the individual inherited
the largest proportion of their genome. An individual is randomly chosen.
The proportions of inheritance are interpreted as probabilities for determining
which sub-population the mate will derive from.

.. code-block:: python
   :caption: Module imports

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='short', quiet=True)
   >>> import simuPOP as sim
   >>> import numpy as np, pandas as pd
   >>> from scipy import stats
   >>> from saegus import parameters, breed

We will continue to use the same population as the rest of our examples.

.. code-block:: python
   :caption: Load the population from ``example_pop.pop``

   >>> example_pop = sim.LoadPopulation('example_pop.pop')
   >>>