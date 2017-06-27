
.. _creating_multi_parental_populations:

###################################
Creating Multi-Parental Populations
###################################

This example will provide a walkthrough for how to create populations similar
to the MAGIC protocol. In this example we are not using inbred lines so it is
a slight misnomer to call it the MAGIC protocol; however, when I initially
wrote the code I was not aware of the difference.

.. _module_imports:

.. code-block:: python
   :caption: Required modules

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='short', quiet=True)
   >>> import simuPOP as sim
   >>> import pandas as pd, import numpy as np
   >>> from saegus import breed

