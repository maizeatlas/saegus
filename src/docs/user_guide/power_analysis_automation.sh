#!/usr/bin/env bash

MACHINE=`uname -a`
PY_VERSION=`python - V`
R_VERSION=`Rscript --version`
JAVAC_VERSION=`javac --version`
TASSEL_LOC="/home/vakanas/tassel-5-standalone"
$CONFIG_FILE_NAME="example_gwas_pipeline.xml"


python example_analysis_automation.py
/run_pipeline.pl -Xmx4g -configFile
Rscript 

