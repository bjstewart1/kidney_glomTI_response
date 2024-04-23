#perfusion analysis
import os
os.chdir('/lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response/code')
import anndata as ad
import scanpy as sc
import os
import pandas as pd
import numpy as np
sc.set_figure_params(figsize=(5,5), dpi = 150, fontsize = 10)
import sys
print(sys.executable)
import matplotlib as mpl
import matplotlib.pyplot as plt


from datetime import datetime
print(datetime.now(tz=None))

import warnings
warnings.filterwarnings('ignore')
#this will hide the awful pandas deprec warnings that are currently plaguing scanpy

#now print versions
print(sc.__version__)
print(ad.__version__)
print(pd.__version__)

os.chdir('/lustre/scratch126/cellgen/team297/bs16/current_projects/kidney_glomTI_response')

cells_use = pd.read_csv("data/da_results/perfusion/DE_testing_cells.csv", sep = '\t')

#set up some paths
cpdb_file_path = 'data/cellphoneDB/db/v4.1.0/cellphonedb.zip'
out_path = 'data/cellphoneDB/perfusion'

adata = sc.read_h5ad("./data/annotated/scRNAseq_perturbation_kidney.h5ad")
#adata = adata[cells_use['cells']]

#get raw count data from adata_ref
def get_counts_raw(adata):
    adata_counts = adata.uns['raw_adata'].copy()
    adata_counts = adata_counts[adata.obs_names]
    adata = adata.raw.to_adata()
    adata.layers = adata_counts.layers.copy()
    return(adata)

adata = get_counts_raw(adata)

adata.X = adata.layers['raw_counts']

#normalise for cellphoneDB
sc.pp.normalize_total(adata, target_sum = 1e4)

cts = np.array(adata.obs['cell_type'])

idx = adata.obs_names.isin(cells_use['cells'])
cts[idx]  = cts[idx]  + "_" + np.array(adata.obs['stimulation'][idx])

adata.obs['cell_type_stim'] = pd.Categorical(cts)

#save the adata
adata.write_h5ad("./data/cellphoneDB/perfusion/input_adata.h5ad")

#make metadata
metadata = pd.DataFrame({'Cell': adata.obs['barcode'], 'cell_type': adata.obs['cell_type_stim']}) #this has the stimulation data... 
metadata.to_csv('data/cellphoneDB/perfusion/metadata.tsv', sep = '\t')

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cpdb_file_path,
        meta_file_path = 'data/cellphoneDB/perfusion/metadata.tsv',
        counts_file_path =  "./data/cellphoneDB/perfusion/input_adata.h5ad",      
        counts_data = 'hgnc_symbol',
        threshold = 0.1,
        output_path = out_path,
        output_suffix = 'perfusion',
        separator = '|',
        result_precision = 3
        )
