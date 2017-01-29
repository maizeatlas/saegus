#! /usr/lib/R/bin/R

library(qvalue)
library(ggplot2)
library(gap)

setwd("/home/vakanas/tassel-5-standalone/output")

input_file_name = "epsilon_0_out_2.txt"
output_file_name = "epsilon_0_qvalues.txt"

results_header = scan(input_file_name, what="character", nlines=1, sep="\t")
gwas_results = read.table(input_file_name, header=F, row.names = NULL, skip=2)
colnames(gwas_results) = results_header
pvalues = gwas_results$p
qobj = qvalue(p = pvalues)
qvalues = data.frame(qobj$qvalues)
colnames(qvalues) = "q"
write.table(qvalues, output_file_name, sep="\t")

