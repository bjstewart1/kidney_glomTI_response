#glom niche analysis
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


#set up some paths
cpdb_file_path = 'data/cellphoneDB/db/v4.1.0/cellphonedb.zip'
out_path = 'data/cellphoneDB/glom'

adata = sc.read_h5ad("./data/annotated/scRNAseq_glom_TI_kidney.h5ad")
adata = adata[adata.obs['tissue'].isin(['Glomerulus_single_cell_suspension'])] #subset to GLOM


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

#save the adata
adata.write_h5ad("./data/cellphoneDB/glom/input_adata.h5ad")


#make metadata
metadata = pd.DataFrame({'Cell': adata.obs['barcode'], 'cell_type': adata.obs['broad_cell_type']})
metadata.to_csv('data/cellphoneDB/glom/metadata.tsv', sep = '\t')


#run cellphonedb
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cpdb_file_path,
        meta_file_path = 'data/cellphoneDB/glom/metadata.tsv',
        counts_file_path = "./data/cellphoneDB/glom/input_adata.h5ad",       
    # mandatory: normalized count matrix.
        counts_data = 'hgnc_symbol',             
    # defines the gene annotation in counts matrix.
        output_path = out_path,
        separator = '|',
        threshold = 0.1,                           
    # defines the min % of cells expressing a gene for this to be employed in the analysis.
        result_precision = 3,                    
    # Sets the rounding for the mean values in significan_means.
        debug = False,                           
    # Saves all intermediate tables emplyed during the analysis in pkl format.
        output_suffix = 'glom')
