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
we can use the shell to loop over each configuration file. The shell script
is contained in ``simulated_mlm.sh``.

.. code-block:: bash
   :caption: Contents of the ``simulated_mlm.sh`` script

   #!/bin/bash


   echo "Run ID: $1, Number of Replicates $2"
   run_id=$1
   number_of_replicates=$2
   final_rep_index="$((number_of_replicates - 1))"

   echo "Beginning TASSEL analysis of Run ID: $run_id"
   echo "Number of Replicates: $number_of_replicates"
   echo "First configuration file: small_0_gwas_pipeline.xml"onca
   echo "Final configuration file: small_"$final_rep_index"_gwas_pipeline.xml"

   for i in `seq 0 $final_rep_index`
   do
       config_file_name=$run_id$i"_gwas_pipeline.xml"
       echo "$config_file_name"
       ./run_pipeline.pl -Xmx6g -configFile $config_file_name
   done


The input is the ``run_id`` followed by the number of replicates. For example:

.. code-block:: sh

   ~/tassel-5-standalone$ bash simulated_mlm.sh small 10

Which iterates through all of the TASSEL configuration files.

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

The result from this block of code is a text file which has two columns:
a column of the indexes corresponding to the loci and the Q values
corresponding to the ``add_p`` column in the TASSEL output file. The output
file is titled ``<run_id_prefix>_qvalues.txt``.


Analysis in Python
##################




.. _qvalue: https://github.com/StoreyLab/qvalue