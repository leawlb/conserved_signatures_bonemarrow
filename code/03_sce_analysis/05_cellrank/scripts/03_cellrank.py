# using cellrank to obtain probabilities and branches from pseudotime
# https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/general/100_getting_started.html

# use ../../envs/cellrank.yml

#-------------------------------------------------------------------------------

import sys

import cellrank as cr
import scanpy as sc
import numpy as np

sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2

import warnings

warnings.simplefilter("ignore", category=UserWarning)

#-------------------------------------------------------------------------------

adata = sc.read_h5ad(snakemake.input['adata_input'])

n_states_total = snakemake.params['n_states_total']
n_states_terminal = snakemake.params['n_states_terminal']

#-------------------------------------------------------------------------------
# use Pseudotime kernel and pseudotime stored in coldata slot "dpt_pseudotime".
from cellrank.kernels import PseudotimeKernel

pk = PseudotimeKernel(adata, time_key="dpt_pseudotime")

#-------------------------------------------------------------------------------
# compute transition matrix between cells
pk.compute_transition_matrix()

#-------------------------------------------------------------------------------
# get cell states using an estimator
from cellrank.estimators import GPCCA

g = GPCCA(pk)
g.fit(n_states=n_states_total, cluster_key="celltypes")

#-------------------------------------------------------------------------------
# define terminal states 
g.predict_terminal_states(method="top_n", n_states=n_states_terminal)
#-------------------------------------------------------------------------------
# for each defines terminal state, get fate probabilities
# this will be added to adata automatically
g.compute_fate_probabilities(n_jobs=1) # trying to avoid weird bug

#-------------------------------------------------------------------------------
# add info from adata.obsm into coldata for nicer conversion to SCE later

print(adata.obsm['lineages_fwd'])
#adata.obs['Erythroid'] = adata.obsm['lineages_fwd']['Erythroid progs.']
#adata.obs['Lymphoid'] = adata.obsm['lineages_fwd']['Lymphoid progs._2']
#adata.obs['Neutro'] = adata.obsm['lineages_fwd']['Neutrophil progs.']

#-------------------------------------------------------------------------------
print(adata)
adata.write_h5ad(snakemake.output['adata_output'])

#-------------------------------------------------------------------------------
# export plots on loaded variables/objects

# https://www.geeksforgeeks.org/save-multiple-matplotlib-figures-in-single-pdf-file-using-python/
import matplotlib 
from matplotlib import pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages 

fig1 = plt.figure() 
sc.pl.embedding(adata, basis="X_umap", color=["celltypes"])

fig2 = plt.figure() 
sc.pl.embedding(adata, basis="X_umap", color=["dpt_pseudotime"])

#fig3 = plt.figure() 
#pk.plot_random_walks(
#    seed=0,
#    n_sims=20,
#    start_ixs={"celltypes": "HSCs"},
#    stop_ixs={"celltypes": ["Erythroid progs.", "Neutrophil progs."]},
#    basis="X_umap",
#    legend_loc="right",
#    dpi=150
#)

fig7 = plt.figure() 
g.plot_fate_probabilities(legend_loc="right", basis='X_umap', same_plot=False)

fig4 = plt.figure() 
g.plot_macrostates(which = 'all', basis = "X_umap")

fig5 = plt.figure() 
g.plot_macrostates(which = 'terminal', basis = "X_umap")

fig6 = plt.figure() #
cr.pl.circular_projection(adata, keys="celltypes", legend_loc="right")


#
def save_image(filename): 
    
    p = PdfPages(filename) 

    fig_nums = plt.get_fignums()   
    figs = [plt.figure(n) for n in fig_nums] 
      
    for fig in figs:  
        
      fig.savefig(p, format='pdf')  
      
    p.close()   
  
filename = snakemake.output['pdf_output']
# call the function
save_image(filename)   
