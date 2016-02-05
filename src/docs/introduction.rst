Introduction
============


   ``wgs`` is a Python package which utilizes the general purpose
   forward-time population genetics simulation environment simuPOP_.
   .. _simuPOP: http://simupop.sourceforge.net/
   simuPOP is an excellent package which allows the user to easily
   simulate a wide range of scenarios. Moreover, the user can
   design custom classes and functions to customize simuPOP to any
   population genetics scenario. ``wgs`` is a set of customized
   classes and functions which are tailored to simulate conditions
   of artificial selection in the commercially important *Zea mays*.
   ``wgs`` is primarily concerned with quantitative traits; however,
   support for non-quantitative traits is available in the native simuPOP
   package. **Note** simuPOP seems to be designed more for disease-gene
   associations; hence, it includes support for non-quantitative traits.

   The typical use of ``wgs`` is to take a fairly small number of individuals
   which we call *founders* or *prefounders* and generate a full sized
   population through a customized breeding scheme. The full sized population
   is subjected to a scenario such as recurrent selection or drift. Data
   from the simulation is then analyzed downstream for GWAS, SOSS or some other
   study. At present the simulator was used to generate a null distribution
   by in the *Tuson* population.