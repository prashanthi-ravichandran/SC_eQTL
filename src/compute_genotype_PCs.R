rm(list = ls())
library(gdsfmt)
library(SNPRelate)
library(dplyr)

# Compute the genotype PCs

geno_matrix <- readRDS("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/genotypes/geno_matrix.rds")
coords <- rownames(geno_matrix)
sampleids <- colnames(geno_matrix)
coords <- data.frame(t(data.frame(strsplit(rownames(geno_matrix), ":"))), stringsAsFactors = F)

rownames(coords) <- c(1:dim(coords[1]))
colnames(coords) <- c("chr", "pos")

coords$chr <- as.numeric(coords$chr)
coords$pos <- as.numeric(coords$pos)

filename <- "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/genotypes/geno_all.gds"
plots.dir <- "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/plots/"

snpids <- rownames(geno_matrix)
geno_matrix <- as.matrix(geno_matrix)
snpgdsCreateGeno(filename, genmat = geno_matrix,
                 sample.id = sampleids, snp.id = snpids,
                 snp.chromosome = coords$chr,
                 snp.position = coords$pos,snpfirstdim=TRUE)


# Open the GDS file
genofile <- snpgdsOpen(filename)

set.seed(1000)

#It is suggested to use a pruned set of SNPs which are in 
#approximate linkage equilibrium with each other to avoid the 
#strong influence of SNP clusters in principal component analysis 
#and relatedness analysis. 

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(snpset)

# Run PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

# Calculate percent variance explained by the top PCs
# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Visualize
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the fourth eigenvector
                  EV5 = pca$eigenvect[,5],    # the fifth eigenvector
                  EV6 = pca$eigenvect[,6],    # the sixth eigenvector
                  EV7 = pca$eigenvect[,7],    # the seventh eigenvector
                  stringsAsFactors = FALSE)
head(tab)

# Read in the metadata for individuals 
meta_data_ind <- read.delim("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/meta_data/per_ind/metadata_B_cells.txt_per_ind.txt")
meta_data_ind$ind_cov <- paste0("X",meta_data_ind$ind_cov)
meta_data_ind <- meta_data_ind[meta_data_ind$ind_cov %in% pca$sample.id, ]
pop.cov <- c()
for(isample in tab$sample.id){
  pop.cov <- c(pop.cov,as.character(meta_data_ind$pop_cov[meta_data_ind$ind_cov == isample]))
}

colors <- ifelse(pop.cov == "WHITE", "blue", "red")

pdf(file = paste0(plots.dir, "genotype_PCs/genotype_PC1vsPC2.pdf"))
plot(tab[, 2], tab[ ,3], col = colors, pch = 20, xlab = "PC1", ylab = "PC2")
legend("topright", fill = c("red", "blue"), legend = c("Asian", "White"), cex = 0.8)
dev.off()

# Plot the scree plot 
pdf(file = paste0(plots.dir, "genotype_PCs/scree_plot.pdf"))
plot(pc.percent, pch = 20, xlim = c(1,20), xlab = "PC Number", ylab = "% Variance explained", main = "Scree plot")
dev.off()

pdf(file = paste0(plots.dir, "genotype_PCs/top7_PCs.pdf"))
plot(tab[2:8], col = colors, pch = 20)
legend("topright", fill = c("red", "blue"), legend = c("Asian", "White"), cex = 0.8)
dev.off()

eigen_vectors <- pca[["eigenvect"]]
rownames(eigen_vectors) <- pca$sample.id
eigen_vectors <- eigen_vectors[ ,1:7]
colnames(eigen_vectors) <- c("V1","V2","V3","V4","V5","V6","V7")

saveRDS(eigen_vectors, "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/cov_pr/genotype_PCs.rds")