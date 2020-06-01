## ---------------------------
##
## Script name: call_eQTLs
##
## Purpose of script: This script uses Matrix eQTL to call eQTLs for each cell type individually from processed pseudobulk data
## 
##
## Author: Prashanthi Ravichandran
##
## Date Created: 2020-01-23
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

rm(list = ls())
library(MatrixEQTL)
library(dplyr)

base.dir <- "/work-zfs/abattle4/prashanthi/Single_cell_eQTL"

## Settings
inputArgs <- commandArgs(TRUE)

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
cell_type <- inputArgs[1]
chr <- inputArgs[2]
#cell_type <- "B_cells"
#chr <- 1
SNP_file_name = paste0(base.dir, "/data/genotypes/chr",as.character(chr),".genotypes.matrix.eqtl.txt")
snps_location_file_name = paste(base.dir, "/data/snpslocs/snpspos_chr",as.character(chr),".txt", sep="")
gene_location_file_name = paste(base.dir, "/data/genelocs/",cell_type, ".txt", sep="")

# Gene expression file name
expression_file_name = paste(base.dir, "/data/expr/",cell_type,".txt", sep="")


# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "/data/covariates/",cell_type,".txt", sep="")

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the space character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

cvrt$columnNames <- snps$columnNames
gene$columnNames <- snps$columnNames

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = F)

snpspos$chr = gsub("chr", "", snpspos$chr)
snpspos$chr = as.numeric(snpspos$chr)

# Check for minor allele frequency
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(sl:.*ice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)
length(maf[maf < 0.05])
## Look at the distribution of MAF
snps$RowReorder(maf>0.05)

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
#cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

## Plot the Q-Q plot of local p-values
pdf(file = paste0("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/plots/eQTL_calling_simple/",cell_type,"/chr",as.character(chr), ".pdf"))
plot(me, pch = 16, cex = 0.5, ylim = c(0, 20))
dev.off()

# Perform multiple testing correction 
summary_stats <- me$cis$eqtls

# Save the results from Matrix eQTL analysis
saveRDS(summary_stats, 
        paste0("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/results/eQTL_calling_simple/",cell_type,"/raw_summary_stats/chr",as.character(chr), ".rds"))

summary_stats <- group_by(summary_stats, gene)
unique_genes <- unique(as.character(summary_stats$gene))
summary_stats$gene <- as.character(summary_stats$gene)
summary_stats$snps <- as.character(summary_stats$snps)
pvalues_by_gene <- list()
for(igene in unique_genes){
  p_values <- summary_stats$pvalue[summary_stats$gene == igene]
  names(p_values) <- summary_stats$snps[summary_stats$gene == igene]
  pvalues_by_gene[[igene]] <- p_values
}

# Now get bonferroni corected pvalues
for(igene in unique_genes){
  pvalues_by_gene[[igene]] <- p.adjust(pvalues_by_gene[[igene]], method = "bonferroni")
}

# Now for each gene get the smallest p-value and corresponding SNP
smallest_p_value <- c()
for(igene in unique_genes){
  smallest_p_value <- c(smallest_p_value, pvalues_by_gene[[igene]][which.min(pvalues_by_gene[[igene]])])
}

# Summarize
final_summary_stats <- data.frame(unique_genes, names(smallest_p_value), smallest_p_value)
colnames(final_summary_stats) <- c("gene", "snp", "p_value")
final_summary_stats <- cbind(final_summary_stats, p.adjust(final_summary_stats$p_value, method = "BH"))
colnames(final_summary_stats) <- c("gene", "snp", "p_value", "adj_p_value")

saveRDS(final_summary_stats, 
        paste0("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/results/eQTL_calling_simple/",cell_type,"/final_summary_stats/chr",as.character(chr), ".rds"))



