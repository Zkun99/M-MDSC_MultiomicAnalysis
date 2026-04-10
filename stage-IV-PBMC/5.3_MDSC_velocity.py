# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 14:36:47 2024

@author: wllab_bioinformatics
"""

import os
import anndata
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

scv.set_figure_params()
plopath = "E:/NSCLC_PBMC_scRNA_seq/RNA_velocity"
prefix = "LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_"

pd.set_option('display.max_columns',70)
pd.set_option('display.width', 50)
pd.set_option('display.max_colwidth',50)

sc.settings.set_figure_params(dpi_save=300,facecolor="white",frameon=False,figsize=(10,10))
plt.rcParams['axes.grid'] = False
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


sample_obs = pd.read_csv("E:/NSCLC_PBMC_scRNA_seq/RNA_velocity/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_cellID_obs.csv")
umap = pd.read_csv("E:/NSCLC_PBMC_scRNA_seq/RNA_velocity/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_cell_embeddings.csv")
cell_clusters = pd.read_csv("E:/NSCLC_PBMC_scRNA_seq/RNA_velocity/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_cell_clusters.csv")
cell_celltype = pd.read_csv("E:/NSCLC_PBMC_scRNA_seq/RNA_velocity/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_cell_celltype.csv")

base_path = "G:/RNA_velocyto/NSCLC_late_12_sample/"

folders = ["post_0_F", "post_2_F", "post_3_F", "post_4_F", "post_5_F", "post_6_F", 
           "pre_0_F", "pre_2_F", "pre_3_F", "pre_4_F", "pre_5_F", "pre_6_F"]

adatas = []
for folder in folders:
    folder_path = os.path.join(base_path, folder, "velocyto")
    
    loom_files = [f for f in os.listdir(folder_path) if f.endswith('.loom')]
    
    if loom_files:
        loom_file = loom_files[0]  
        loom_path = os.path.join(folder_path, loom_file)
        
        adata = anndata.read_loom(loom_path)
        adata.var_names_make_unique()
        
        adata.obs = adata.obs.rename(index=lambda x: x.replace(':', '_').replace('x', ''))
        
        adata = adata[np.isin(adata.obs.index, sample_obs["meta_RNA_velo$barcodes_velocity"])]
        
        adatas.append(adata)

combined_adata = adatas[0].concatenate(*adatas[1:])
print(combined_adata)

adata = combined_adata
adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {"CellID":'barcodes_velocity'})
adata_index.barcodes_velocity.head()

rep=lambda x : x.split("-")[0]
adata_index["barcodes_velocity"]=adata_index["barcodes_velocity"].apply(rep)


umap = umap[np.isin(umap["barcodes_velocity"],adata_index["barcodes_velocity"])] 
umap = umap.drop_duplicates(subset=["barcodes_velocity"]) 
umap_ordered = adata_index.merge(umap, on = "barcodes_velocity")
umap_ordered = umap_ordered.iloc[:,1:] 
adata.obsm['X_umap'] = umap_ordered.values 


cell_clusters = cell_clusters[np.isin(cell_clusters["barcodes_velocity"],adata_index["barcodes_velocity"])]
cell_clusters = cell_clusters.drop_duplicates(subset=["barcodes_velocity"]) 
cell_clusters_ordered = adata_index.merge(cell_clusters, on = "barcodes_velocity")
cell_clusters_ordered = cell_clusters_ordered.iloc[:,1:]
adata.obs['clusters'] = cell_clusters_ordered.values

# Running RNA Velocity
scv.pp.filter_and_normalize(adata)#,min_shared_counts=30, n_top_genes=2000
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)


#adata.write('E:/NSCLC_PBMC_scRNA_seq/RNA_velocity/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_velo_done.h5ad', compression = 'gzip')
adata= scv.read('E:/NSCLC_PBMC_scRNA_seq/RNA_velocity/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_velo_done.h5ad')

##base on the zemin zhang pan-B 2024 setting
scv.pl.velocity_embedding_grid(adata, basis='X_umap',color = "clusters", 
                               palette = ident_colours,
                               arrow_size=3,
                               arrow_length=3,
                               density=0.4,
                               arrow_color="black",
                               legend_loc= "right margin",
                               legend_fontsize=5,
                               save=os.path.join(plopath,prefix +"MDSC_UMAP_stream_250213.pdf"), 
                               figsize=(5,5))

