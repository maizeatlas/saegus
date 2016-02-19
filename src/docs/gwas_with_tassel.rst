Generating Data for GWAS with TASSEL
====================================

TASSEL requires three files for GWAS: phenotypes, kinship and a hapmap of
genotypes (population structure is optional). For most runs we will use all
four files unless otherwise noted.

Command Line Interface to TASSEL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The TASSEL command line interface requires a considerable number of
options to run GWAS. It is impractical to run the command line manually
for the number of replications in a simulated study. Here is an example of
the command required for a single run of TASSEL

.. code-block:: console

    run_pipeline.bat -fork1 -h C:\GWAS\rs_run_1_simulated_hapmap.txt -fork2
    -t C:\GWAS\rs_run_1_phenotype_vector.txt -fork3 -q
    C:\GWAS\rs_run_1_structure_matrix.txt -fork4 -k
    C:\GWAS\rs_run_1_kinship_matrix.txt -combine5 -input1 -input2 -input3
    -intersect -combine6 -input5 -input4 -mlm -mlmCompressionLevel None
    -export C:\GWAS\result\rs_run_1_ -runfork1 -runfork2 -runfork3 -runfork4

Fortunately TASSEL's command line interface is able to use a .xml file
containing the same information for input. Using an XML file only requires a
single option to the run_pipeline.bat script. For example if we have a file
called rs_run_1_gwas_pipeline.xml we can use the -configFile option.

.. code-block:: console

    run_pipeline.bat -configFile rs_run_1_gwas_pipeline.xml

Using a .xml file is much cleaner than specifying all of those options in the
terminal. Moreover, a .xml file also elucidates how TASSEL is processing the
input files.

Example XML File for MLM GWAS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: xml

    <?xml version="1.0" encoding="UTF-8" standalone="no" ?>
        <TasselPipeline>
            <fork1>
                <h>C:\GWAS\rs_run_1_simulated_hapmap.txt</h>
            </fork1>
            <fork2>
                <t>C:\GWAS\rs_run_1_phenotype_vector.txt</t>
            </fork2>
            <fork3>
                <q>C:\GWAS\rs_run_1_structure_matrix.txt</q>
            </fork3>
            <fork4>
                <k>C:\GWAS\rs_run_1_kinship_matrix.txt</k>
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
                <export>C:\GWAS\result\rs_run_1_</export>
            </combine6>
            <runfork1/>
            <runfork2/>
            <runfork3/>
            <runfork4/>
        </TasselPipeline>


Running GWAS With TASSEL
^^^^^^^^^^^^^^^^^^^^^^^^

Performing exploratory analyses with simulated populations is divided into
three sets of tasks:

    - simulation
    - GWAS
    - analysis

A script which uses saegus is performs the simulation. A single set of
founders is used to generate multiple replicates of a population under
selection or drift or both. The script also generates additional files
containing:

    - configuration files used to run GWAS with TASSEL
    - parameter information
    - mean and variance data
    - minor allele matrix


