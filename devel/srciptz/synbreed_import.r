rm(list=ls(all=TRUE))

library(synbreed)
setwd("C:\\Users\\DoubleDanks\\BISB\\wisser\\code")



geno <- read.table('syn_genotype.txt', header=TRUE, sep='\t', row.names=1)
pheno <- read.table('syn_phenotype.txt', header=TRUE, sep='\t', row.names=1)
map <- read.table('syn_marker_map.txt', header=TRUE, sep='\t', row.names=1)
gp <- create.gpData(pheno, geno, map)
gp2 <- codeGeno(gp, label.heter="1")
G <- kin(gp2, ret="realized")
plot(G)
write.csv(G, "first_wgs_relationship_matrix.txt")

# Individuals on the diagonal of the relationship matrix G are compared
# to themselves. Hence the diagonal of G is much darker when plotted
