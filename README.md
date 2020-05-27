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
 
 - 
