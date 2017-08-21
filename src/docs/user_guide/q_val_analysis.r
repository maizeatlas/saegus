library(qvalue)

setwd("/home/vakanas/tassel-5-standalone")
results_header = scan("example_2.txt", what = "character", nlines=1, sep="\t")
gwas_results = read.table("example_2.txt", header = F, row.names = NULL, skip=2)
colnames(gwas_results) = results_header

add_pvalues = gwas_results$add_p)
padd_qobj = qvalue(p = add_pvalues)

allele_data = h5file("/home/vakanas/BISB/rjwlab-scripts/saegus_project/src/docs/user_guide/8517example_trait_data.hdf5")
allele_effects = allele_data["effects"]
allele_states = allele_data["states"]
allele_frequencies = allele_data["frequencies"]
segloci = allele_data["segregating_loci"]

write.table(padd_obj$qvalues, "example_q_matrix.txt", sep="\t", quote = FALSE, col.names = FALSE)
