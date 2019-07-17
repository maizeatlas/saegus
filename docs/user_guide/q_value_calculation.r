#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test to determine if the file name parameter is supplied to the script
if (length(args)==0) {
  stop("At least one argument must be suppled (input file).\n", call.=FALSE)
}
library(qvalue)
library(ggplot2)
library(gap)
#setwd("/home/vakanas/tassel-5-standalone/output")

run_id = args[1]
file_name_match_pattern = paste(run_id, "(.*)_2.txt", sep='')
file_names = list.files(pattern = file_name_match_pattern)

for(n in file_names) {
    print(n)
    input_file_name = n
    run_id_prefix_terminus = nchar(input_file_name) - 5
    run_id_prefix = substring(input_file_name, 1, run_id_prefix_terminus)
    output_file_name = paste(run_id_prefix, 'q_values.txt', sep='')
    print(output_file_name)
    results_header = scan(input_file_name, what="character", nlines=1, sep="\t")
    gwas_results = read.table(input_file_name, header=F, row.names = NULL, skip=2)
    colnames(gwas_results) = results_header
    pvalues = gwas_results$p
    qobj = qvalue(p = pvalues)
    qvalues = data.frame(qobj$qvalues)
    colnames(qvalues) = "q"
    rownames(qvalues) = gwas_results$Marker
    write.table(qvalues, output_file_name, sep="\t")
}