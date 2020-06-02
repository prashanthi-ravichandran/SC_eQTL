rm(list = ls())

raw_res_dir <- "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/results/eQTL_calling_simple/"
significant_res_dir <- "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/results/eQTL_significant_all/"
cell_type <- "FCGR3A_monocytes"
for(i in c(1:22)){
  if(i == 1){
    all_summary_stats <- readRDS(paste0(raw_res_dir, cell_type, "/final_summary_stats/chr", i, ".rds"))
  }else{
    all_summary_stats <- rbind(all_summary_stats, readRDS(paste0(raw_res_dir, cell_type, "/final_summary_stats/chr", i, ".rds")))
  }
}

all_summary_stats <- all_summary_stats[all_summary_stats$adj_p_value < 0.05, ]
saveRDS(all_summary_stats, paste0(significant_res_dir, cell_type, ".rds"))

# Find MAF matched SNPs

