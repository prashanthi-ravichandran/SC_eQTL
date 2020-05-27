# Single cell eQTL 

- First run **Process_per_cell_type.ipynb** which reads in the raw data from *CLUESImmVar_nonorm.V6.h5ad* and splits the samples based on the provided cell-types and retains only those samples that belong to individuals with Lupus. Then the following QC steps are performed
  - Remove cells with fewer than 200 genes expressed
  - Remove genes expressed in fewer than 3 cells
  - Remove cells with % mito reads greater than 0.05
  - Remove cells with more than 2500 genes expressed
  - We then subset to autosomal, protein-coding genes, and recalculate the QC metrics, the eQTL sepecific filters are then performed which include
  - Remove cells with fewer than 400 expressed genes
  - Remove genes expressed in fewer than 5% of all cells
 This script also computes the log1p normalized and scaled single cell data 
 
 - Run **processing_pseudobulk.R** to compute 
   - pseudobulk by summing counts from all cells belonging to a specific individual 
   - compute expression PCs 
   - compute genotype PCs
   - create covariates.txt file for downstream eQTL calling 
  
 - Run **call_eQTLs.R** for each combination of cell type and chromosome. This uses matrix-eQTL to call cis-eQTLs by testing all SNPs that are within a 1Mb distance from the gene TSS that have a MAF > 0.05. Further, for each gene, we select the gene-SNP pair that has the lowest p-value after Bonferroni correcting within each gene and then perform BH correction across genes to obtain significant SNPs at FDR of 0.05. 
