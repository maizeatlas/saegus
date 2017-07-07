============
Introduction
============


A huge number of functions in this module are obsolete. Their existence reflects
the successive iterations of improvement. Many functions exist simply because I
was too sick to write proper documentation as I developed the code.

.. attention:: JJD 6/30/2017

Welcome to the `saegus` documentation! This page serves as an introduction for
some fairly basic information

Concepts of simuPOP
===================

The saegus package is built on top of the incredibly powerful package simuPOP
The power of simuPOP_ comes at a price: intellectual burden to learn
how to write scripts to utilize the simuPOP_ machinery. In this introduction
we define some basic terminology and concepts required to understand the rest
of the documentation.

What makes saegus unique is the way in which simuPOP's basic components are
put together as well as a heavy degree of customization for quantitative
traits.

.. _simuPOP: http://simupop.sourceforge.net/Main/HomePage

Evolutionary Process
====================

simuPOP performs all mating processes inside of a  call to evolve
evolve has five intervals:

#) [initOps] Occurs before anything and only occurs once
#) [preOps] Occurs before mating
#) [matingScheme] Defines how to mate individuals and produce offspring
#) [postOpts] Occurs after mating
#) [finalOps] Absolute last step and only occurs once.

At each step all of the functions listed are called to operate on the Population
or Population replicates in the case of a Simulator object.

Operators, Functions and Customization
======================================

There is almost no difference between a function and an Operator
A simuPOP Operator is a function applied during a call to evolve.
A simuPOP function is a function which is applied outside of a call to evolve.
All Operators can be turned into functions so long as the operations are legitimate
outside of a call to evolve.

simuPOP provides a large amount of highly memory optimized **Operator** and functions
with native support for parallel computation. Even better simuPOP allows the user
to make custom **Operator** written entirely in Python. Thus simuPOP can be tailored
to literally any situation in forward-time population genetic simulations.