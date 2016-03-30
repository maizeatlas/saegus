============
Introduction
============

   The saegus package is built on top of the incredibly powerful package simuPOP_.
   The power of simuPOP_ comes at a price: intellectual burden to learn
   how to write scripts to utilize the ``simuPOP`` machinery. In this introduction
   we define some basic terminology and concepts required to understand the rest
   of the documentation.

   What makes saegus unique is the way in which simuPOP's basic components are
   put together as well as a heavy degree of customization for quantitative
   trait models. The :mod:`analyze` is an entirely independent set of classes and
   functions which are dedicated to analyzing, formatting and storing the results
   of simulations.

Evolutionary Process
====================

   simuPOP performs all population genetic processes such as mating, migration and selection
   inside of a :method:`evolve`. :method:`evolve` has five intervals:

      #) [initOps] Occurs before anything and only occurs once
      #) [preOps] Occurs before mating
      #) [matingScheme] Defines how to mate individuals and produce offspring
      #) [postOpts] Occurs after mating
      #) [finalOps] Absolute last step and only occurs once.

Operators, Functions and Customization
======================================



