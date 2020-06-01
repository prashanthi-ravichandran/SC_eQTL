import numpy as np
import tensorflow as tf
import scipy.special as sp
import scipy.stats as st
import scqtl
import pandas as pd
import os
import json
import pickle
import argparse as ap

def recode(annotations, key):
  n = annotations.shape[0]
  cat = sorted(set(annotations[key]))
  onehot = np.zeros((n, len(cat)))
  onehot[np.arange(n), annotations[key].apply(cat.index)] = 1
  return onehot

def get_args():
	p = ap.ArgumentParser()
	p.add_argument("--cell-type", type=str, required=True)
	p.add_argument("--n-genes", type=int, default=200, help="Number of genes used in each subset")
	p.add_argument("--learning-rate", type=float, default=1e-4)
	p.add_argument("--train-epochs", type=int, default=1e5)
	return p.parse_args()	


args = get_args()

save_dir = "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/UMI_counts/expr/"
expr = pd.read_csv(save_dir + args.cell_type + ".csv")

pseudobulk = pd.read_csv("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/expr/" + args.cell_type  + ".txt", sep = "\t")
pseudobulk = pseudobulk[~(pseudobulk == 0).any(axis=1)]
expr = expr.loc[:, np.insert(pseudobulk.index.values, 0, "index", axis=0)]

meta_data = pd.read_csv("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/UMI_counts/metadata/" + args.cell_type + ".csv")
print(meta_data['index'].equals(expr['index']))
meta_data.head()


# Examine the meta data distributions
print(meta_data.batch.value_counts())
print(meta_data.batch_cov.value_counts())
print(meta_data.pop_cov.value_counts())
print(meta_data.well.value_counts())

umi = expr.iloc[:, 1:].values
umi

size_factor = pd.read_csv("/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data/size_factor/" + args.cell_type +".csv")
print(size_factor['index'].equals(expr['index']))
size_factor = size_factor['total_counts'].values.reshape((umi.shape[0],1))

batch     = recode(meta_data, "batch")
batch_cov = recode(meta_data, "batch_cov")
well      = recode(meta_data, "well")
pop_cov   = recode(meta_data, "pop_cov")
design    = np.concatenate((batch, batch_cov, well, pop_cov), axis = 1)
design -= design.mean(axis=0)
onehot = recode(meta_data, "ind_cov")

# Create variables to store the inferred log_mu and log_phi values
log_mu = np.empty([onehot.shape[1], umi.shape[1]])
log_phi = np.empty([onehot.shape[1], umi.shape[1]])
logodds = np.empty([onehot.shape[1], umi.shape[1]])
nb_llik = np.empty([onehot.shape[1]])
zinb_llik = np.empty([onehot.shape[1]])

print(design.shape)
print(umi.shape)
print(onehot.shape)
print(size_factor.shape)

for i in range(0, umi.shape[1], args.n_genes):
    if((i + args.n_genes) < umi.shape[1]):
        start_index = i
        end_index = i + args.n_genes
    else:
        start_index = i
        end_index = umi.shape[1]
    umi_sub = umi[:, start_index:end_index]
    print(umi_sub.shape)
    init = scqtl.tf.fit(
      umi=umi_sub.astype(np.float32),
      onehot=onehot.astype(np.float32),
      design=design.astype(np.float32),
      size_factor=size_factor.astype(np.float32),
      learning_rate=args.learning_rate,
      max_epochs=args.train_epochs,
      verbose=True)
    log_mu[:,start_index:end_index], log_phi[:,start_index:end_index], logodds[:,start_index:end_index], nb_llik[start_index:end_index], zinb_llik[start_index:end_index] = scqtl.tf.fit(
      umi=umi_sub.astype(np.float32),
      onehot=onehot.astype(np.float32),
      design=design.astype(np.float32),
      size_factor=size_factor.astype(np.float32),
      learning_rate=args.learning_rate,
      max_epochs=args.train_epochs,
      warm_start=init[:3],
      verbose=True)

res_dir = "/work-zfs/abattle4/prashanthi/Single_cell_eQTL/results/ZINB_summary_stats/"+ args.cell_type
np.savetxt((res_dir + "/log_mu.csv"), log_mu, delimiter=",")
np.savetxt((res_dir + "/log_phi.csv"), log_phi, delimiter=",")
np.savetxt((res_dir + "/logodds.csv"), logodds, delimiter=",")
np.savetxt((res_dir + "/nb_llik.csv"), nb_llik, delimiter=",")
np.savetxt((res_dir + "/zinb_llik.csv"), zinb_llik, delimiter=",")



