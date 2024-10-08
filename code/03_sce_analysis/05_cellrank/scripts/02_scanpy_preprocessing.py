#-------------------------------------------------------------------------------

# pre-process hsc adata for use with cellrank
# from cellrank and scanpy tutorials:
# https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/general/100_getting_started.html
# https://scanpy.readthedocs.io

# use ../../envs/cellrank.yml

#-------------------------------------------------------------------------------

import sys


import scanpy as sc
import numpy as np
# import random
# 
# random.seed(37)

#-------------------------------------------------------------------------------

adata = sc.read_h5ad(snakemake.input['adata_input'])

#-------------------------------------------------------------------------------
# Calculate neighbors according to scanpy tutorials based on PCA after BC
# use random_state=0 as mentioned in cellrank tutorial

# check out if the number 50 relates to PCs or neighbors!
sc.pp.neighbors(adata, 50, use_rep = 'PCA', random_state=0)

#-------------------------------------------------------------------------------
# calculate UMAP for visualisation in python
sc.tl.umap(adata)

#-------------------------------------------------------------------------------
# get diffusion map
sc.tl.diffmap(adata, n_comps=15)

#-------------------------------------------------------------------------------
# calculate diffusion pseudotime 
# define root cell type

adata.uns['iroot'] = np.flatnonzero(adata.obs['celltypes'] == 'HSCs')[0]
sc.tl.dpt(adata)
print(adata)

#-------------------------------------------------------------------------------
# save
adata.write_h5ad(snakemake.output['adata_output'])
