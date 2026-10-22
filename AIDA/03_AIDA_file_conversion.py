# sbatch --mem=120GB --time=3:15:00 --wrap="module load python/3.12.1; python 03_AIDA_file_conversion.py"
import pandas as pd
import scanpy as sc
from scipy import io


adata = sc.read_h5ad('/proj/yunligrp/users/chenweit/data/AIDA2023/local.h5ad')

donor_ids = adata.obs['donor_id']

# get unique donor_ids
donor_ids = donor_ids.unique()

donor_ids = pd.Series(donor_ids)

reference_ids = donor_ids.sample(n = 200, replace=False)

donor_ids = donor_ids[~donor_ids.isin(reference_ids)]

# save split data. First the reference

adata = adata[adata.obs['donor_id'].isin(reference_ids)]
out_dir = '/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/AIDA/datasets_for_comparison/reference_all'


with open(out_dir + '/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')

with open(out_dir + '/features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')

io.mmwrite(out_dir +'/matrix_raw', adata.raw.X.T)


adata.obs.to_csv(out_dir + '/metadata.csv')

# save split data. Then the test set

adata = sc.read_h5ad('/proj/yunligrp/users/chenweit/data/AIDA2023/local.h5ad')

adata = adata[adata.obs['donor_id'].isin(donor_ids)]
out_dir = '/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/AIDA/datasets_for_comparison/test_set_all'

with open(out_dir + '/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')

with open(out_dir + '/features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')

io.mmwrite(out_dir +'/matrix_raw', adata.raw.X.T)

adata.obs.to_csv(out_dir + '/metadata.csv')
