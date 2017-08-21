.. _power_analysis:

###############################
Power Analysis: Start to Finish
###############################

This tutorial will show how to close the loop. From start in Python
to analysis in Java (TASSEL), calculation of Q-values in R back to Python
to interpret and analyze the results. These distinct phases are automated
via a shell script.

Simulation in Python
####################

The prior sections of the guide go into great detail about specifying
a simulation. We will mention how to put the results in a script i.e.
a ``.py`` file. The output from the simulation is fed into TASSEL via
the ``configFile``.

TASSEL Analysis
###############

After the simulation ceases the data files are written. All we need to do
is to run TASSEL using our configFile. If there are multiple replicates
we can use the shell to loop over each configuration file.

Calculation of Q-values in R
############################

We use the qvalue_ package from ``R`` to obtain the Q-values from the P-values
in the TASSEL results.

.. code-block:: R

   > library(qvalues)
   > results_header = scan("example_2.txt", what = "character", nlines=1, sep="\t")
   > colnames(gwas_results) = results_header
   > marker_column = matrix(gwas_results$Marker)
   > pvalues = gwas_results$p
   > add_pvalues = gwas_results$add_p
   > p_qobj = qvalues(p = pvalues)
   > add_p_qobj = qvalues(p = add_pvalues)
   > q_p = matrix(p_qobj$qvalues)
   > q_add_p = matrix(add_p_obj$qvalues)
   > rownames(q_p) = marker_column
   > rownames(q_add_p) = marker_column





Analysis in Python
##################




.. _qvalue: https://github.com/StoreyLab/qvalue