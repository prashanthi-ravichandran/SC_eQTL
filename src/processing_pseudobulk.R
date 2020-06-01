rm(list = ls())
showfile.gds(closeall=TRUE)
# This script merges cells which belong to the same individual, experiment and day
# Further we format the genotype matrix
# We also compute expression PCs to be used in downstream analysis for our samples

library(rsvd)
library(tidyverse)
library(Seurat)
library(genoscapeRtools)
library(RColorBrewer)
require(data.table)
library(DESeq2)

cell_type <- "megakaryocytes"
dat.dir <- "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/"
plots.dir <- "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/plots/"

umi <- fread(paste0(dat.dir, "UMI_counts/expr/", cell_type, ".csv"), header = T)
rownames(umi) <- umi$index
umi$index <- NULL

metadata <- read.csv(paste0(dat.dir, "UMI_counts/metadata/", cell_type, ".csv"), stringsAsFactors = FALSE)
rownames(metadata) <- metadata$index
metadata$index <- NULL

metadata$ind_cov <- paste0("X", metadata$ind_cov)
uniq_individuals <- colnames(read.table(paste0(dat.dir, "genotypes/chr1.genotypes.matrix.eqtl.txt")))
#uniq_individuals <- gsub("X", "", uniq_individuals)

pop_cov <- c()
batch_cov <- c()
for(ind in uniq_individuals){
  pop_cov <- c(pop_cov, metadata$pop_cov[metadata$ind_cov == ind][1])
  batch_cov <- c(batch_cov, metadata$batch_cov[metadata$ind_cov == ind][1])
}

pseudobulk = sapply(uniq_individuals, function(ind_id){
  cell.ids <- rownames(metadata)[metadata$ind_cov == ind_id]
  sample_count_df = umi[rownames(umi) %in% cell.ids, ]
  sample_count_sums_df = colSums(sample_count_df)
  return(sample_count_sums_df)
})

gene_location_file_name = "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/gencode.v19.annotation.gene.txt"
genepos = read.delim(gene_location_file_name, stringsAsFactors = FALSE)
genepos$chr = gsub("chr", "", genepos$chr)
genepos = distinct(genepos, gene_name, .keep_all = TRUE)
common.genes = intersect(genepos$gene_name, rownames(pseudobulk))

pseudobulk <- pseudobulk[rownames(pseudobulk) %in% common.genes, ]
genepos <- genepos[match(rownames(pseudobulk), genepos$gene_name), ]
genepos_df <- cbind(genepos$gene_name, genepos$chr, genepos$start_pos, genepos$end_pos)
colnames(genepos_df) <- c("geneid", "chr", "left", "right")

write.table(genepos_df,
            file = paste0(dat.dir, "/genelocs/", cell_type,".txt"), quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")


# Normalize the pseudobulk data using rlog from deSEQ
sample_info <- data.frame(pop_cov, batch_cov)
rownames(sample_info) <- uniq_individuals
DESeq.ds <- DESeqDataSetFromMatrix(countData = pseudobulk, colData = sample_info, design = ~ batch_cov)
rld <- rlog(DESeq.ds, blind=TRUE)
#head(assay(rld), 3)
log_normalized  <- assay(rld)

scaled.pseudobulk <- scale(t(log_normalized))
genes.na <- colnames(scaled.pseudobulk)[colSums(is.na(scaled.pseudobulk)) > 0]
scaled.pseudobulk <- scaled.pseudobulk[ ,!colnames(scaled.pseudobulk) %in% genes.na]
expr.usv <- svd(scaled.pseudobulk)
expr_pcs <- expr.usv$u
pc.percent <- ((expr.usv$d^2)/sum((expr.usv$d^2)))*100

expr_save <- t(scaled.pseudobulk)
write.table(expr_save,
            file = paste0(dat.dir, "/expr/", cell_type,".txt"), quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

# To create the covariates file 
geno_matrix <- readRDS("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/genotypes/geno_matrix.rds")
coords <- rownames(geno_matrix)
sampleids <- colnames(geno_matrix)
chr  <- as.numeric(gsub(":.*", "", rownames(geno_matrix)))
pos <- as.numeric(gsub(".*:", "", rownames(geno_matrix)))

filename <- "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/genotypes/geno_all.gds"

snpids <- rownames(geno_matrix)
geno_matrix <- as.matrix(geno_matrix)
snpgdsCreateGeno(filename, genmat = geno_matrix,
                 sample.id = sampleids, snp.id = snpids,
                 snp.chromosome = chr,
                 snp.position = pos,snpfirstdim=TRUE)

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
tab <- data.frame(sample.id = uniq_individuals,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the fourth eigenvector
                  EV5 = pca$eigenvect[,5],    # the fifth eigenvector
                  EV6 = pca$eigenvect[,6],    # the sixth eigenvector
                  EV7 = pca$eigenvect[,7],    # the seventh eigenvector
                  stringsAsFactors = FALSE)
head(tab)

colors <- ifelse(pop_cov == "WHITE", "blue", "red")

pdf(file = paste0(plots.dir, "genotype_PCs/genotype_PC1vsPC2.pdf"))
plot(tab[, 2], tab[ ,3], col = colors, pch = 20, xlab = "PC1", ylab = "PC2")
legend("bottomright", fill = c("red", "blue"), legend = c("Asian", "White"), cex = 0.8)
dev.off()

# Plot the scree plot 
pdf(file = paste0(plots.dir, "genotype_PCs/scree_plot.pdf"))
plot(pc.percent, pch = 20, xlim = c(1,20), xlab = "PC Number", ylab = "% Variance explained", main = "Scree plot")
dev.off()

pdf(file = paste0(plots.dir, "genotype_PCs/top5_PCs.pdf"))
plot(tab[,2:6], col = colors, pch = 20)
legend("topright", fill = c("red", "blue"), legend = c("Asian", "White"), cex = 0.8)
dev.off()

eigen_vectors <- pca[["eigenvect"]]
rownames(eigen_vectors) <- uniq_individuals
eigen_vectors <- eigen_vectors[ ,1:7]
colnames(eigen_vectors) <- c("V1","V2","V3","V4","V5","V6","V7")

covariates <- cbind(eigen_vectors, expr_pcs[ ,1:15])
colnames(covariates) <- c("V1", "V2", "V3", "V4", "V5","V6","V7", "PC1", "PC2", "PC3", 
                          "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
                          "PC11", "PC12", "PC13", "PC14", "PC15") 
covariates <- t(covariates)

write.table(covariates,
            file = paste0(dat.dir, "/covariates/", cell_type,".txt"), quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")








