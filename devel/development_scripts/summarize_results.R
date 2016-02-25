# Prerequisite packages

library(gap)
library(qvalue)
library(dplyr)
library(tidyr)
library(ggplot2)


read_analysis_file = function(file_name, skip_rows){
    header = scan(file_name, what='character', nlines=1, sep='\t')
    analysis_data = read.table(file_name, header=F, row.names=NULL, skip=skip_rows)
    colnames(analysis_data) = header
    return(analysis_data)
}

simtrait_across = read_analysis_file('rs_R0_gwas_out_2.txt', 2)

snp_list = simtrait_across$Marker

qqunif(simtrait_across$p)
pdf("QQ_plot.pdf")



