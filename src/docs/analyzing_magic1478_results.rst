.. analyze-magic1478-rdm-mating-results:

=====================================
Analyzing Results of GWAS with TASSEL
=====================================

We use the mixed-linear modeling method implemented in TASSEL.
The MLM in TASSEL requires three files at minimum with the option for a fourth.

   * genotypes in *hapmap* format
   * kinship matrix
   * phenotypes
   * population structure matrix (optional)

:mod:`saegus` generates all four files. In this particular case the file names
are:

   * ``simulated_hapmap.txt``
   * ``phenotype_vector.txt``
   * ``kinship_matrix.txt``
   * ``structure_matrix.txt``

For the sake of simplicity I placed them in same directory as the
TASSEL executable file: ``sTASSEL.jar``. Running an analysis
with four files using the TASSEL command line interface requires typing
a very long string of code. Fortunately TASSEL allows the user to create
a ``config`` file in ``xml`` format. A ``config`` file provides a much clearer
description of how TASSEL is processing the input. Below is the exact file I used
to run the current analysis.


.. code-block:: xml

   <?xml version='1.0' encoding='UTF-8' standalone='no'?>
   <TasselPipeline>
       <fork1>
           <h>simulated_hapmap.txt</h>
       </fork1>
       <fork2>
           <t>phenotype_vector.txt</t>
       </fork2>
       <fork3>
           <q>structure_matrix.txt</q>
       </fork3>
       <fork4>
           <k>kinship_matrix.txt</k>
       </fork4>
       <combine5>
           <input1/>
           <input2/>
           <input3/>
           <intersect/>
       </combine5>
       <combine6>
           <input5/>
           <input4/>
           <mlm/>
           <mlmCompressionLevel>
               None
           </mlmCompressionLevel>
           <export>gwas_out_</export>
       </combine6>
       <runfork1/>
       <runfork2/>
       <runfork3/>
       <runfork4/>
   </TasselPipeline>

The output from this analysis is a set of three files. The file names
are determined by ``<export>gwas_out_</export>`` so that the exact file names
are:

   * ``gwas_out_1.txt``
   * ``gwas_out_2.txt``
   * ``gwas_out_3.txt``

Given this config file and the input files we use a ``.cmd`` file to run
the analysis in TASSEL:

.. code-block:: bat

   for %%f in (sim_gwas_pipeline.xml) do (
      echo %%~nf
      run_pipeline.bat -configFile "%%f"
      )


The file we are primarily interested in is ``gwas_out_2.txt``. This is the file
which records the values and p-values of all the statistical tests. The goal is
to import the ``gwas_out_2.txt`` into R and determine the Q values using the false
discovery rate to control for false-positives.


Dataset Specifics
=================

This is some useful information about the data-set used in the TASSEL MLM:

Population Size: 2000
QTL: 2, 10, 20

qvalues in R
------------



We will follow Jim's tutorial to use the :mod:`qvalue` package in R; however, I
have found that the function we want to use :func:`qvalue` does not handle
missing data i.e. ``NaN``. Because I am more proficient with ``python`` than
``R`` I used the ``python`` :mod:`pandas` package to convert all ``NaN`` p-values
into values of :math:`1.0`

For example a sample of the P-values of ``gwas_out_2.txt`` are:

+ 0.4968
+ 5.6091E-28
+ NaN
+ 0.6236
+ 0.16525


If we use the :func:`qvalue` function directly it will result in an error.
Instead I use the values:

+ 0.4968
+ 5.6091E-28
+ 1.0
+ 0.6236
+ 0.16525

The edited file name is ``edited_gwas_out_2.txt``. I use these commands to
obtain the q-values.

.. code-block:: R

   results_header = scan("edited_gwas_out_2.txt", what="character", nlines=1, sep="\t")
   gwas_results = read.table("edited_gwas_out_2.txt", header=F, row.names=NULL, skip=2)
   colnames(gwas_results) = results_header

   pvalues = gwas_results$p
   library(qvalue)
   qobj = qvalue(p = pvalues)
   qobj$qvalues
   qvalues_of_magic1478_results = data.frame(qobj$qvalues)
   write.table(qvalues_of_magic1478_results, "qvalues_of_magic1478.txt", sep="\t")