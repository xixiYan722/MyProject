### RNA velocyto analysis of neural cells subset ###

# caculate the unsplice/spliced reads
velocyto run10x -m repeat_masker.gtf  cellranger/  Bos_taurus.ARS-UCD1.2.100.gtf 

# load python packages
import scvelo as scv
import sys
sys.executable ##examine the python version
import scanpy as sc 
import cellrank as cr 
import numpy as np 
import pandas as pd 
import anndata as ad

# read the data
adata = scv.read("neu.h5ad") 
adata
adata.obs.index

# modify barcodes ID to match with spliced/unspliced count data
barcodes = adata.obs.index.tolist()
barcodes = [bc[0:23] for bc in barcodes]
adata.obs.index = barcodes 
adata.obs.index

# read the loom files (store the unsplice/spliced reads)
ldata1 =scv.read('obj1.loom',cache=True)
ldata2 =scv.read('obj2.loom',cache=True)
ldata3 =scv.read('obj3.loom',cache=True)

# modify barcodes ID to match annadata objects to gene expression matrix data
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()] 
barcodes = ['100d_' + bc[0:len(bc)-1] + '-1' for bc in barcodes] 
ldata1.obs.index = barcodes 
ldata1.obs.index
barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()] 
barcodes = ['150d_' + bc[0:len(bc)-1] + '-1' for bc in barcodes] 
ldata2.obs.index = barcodes 
barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()] 
barcodes = ['200d_' + bc[0:len(bc)-1] + '-1' for bc in barcodes] 
ldata3.obs.index = barcodes 

# make variable names unique 
ldata1.var_names_make_unique() 
ldata2.var_names_make_unique() 
ldata3.var_names_make_unique()

# merge ldata
ldata = ldata1.concatenate([ldata2,ldata3])
ldata.obs.index

# modify barcodes ID
barcodes = ldata.obs.index.tolist()
barcodes = [bc[0:23] for bc in barcodes]
ldata.obs.index = barcodes 
ldata.obs.index

# merge matrices into the original adata object 
adata = scv.utils.merge(adata, ldata)

# change dtype from 'object' to 'categroy'
adata.obs['celltype']=adata.obs['celltype'].astype('category').values
adata.obs['celltype']

# examine the proportion of spliced/unspliced in each cell types
scv.pl.proportions(adata,groupby='celltype')

# preprocess
scv.pp.filter_and_normalize(adata) 
scv.pp.moments(adata)

# caculate RNA velocity using a steady-state model(random option)
scv.tl.velocity(adata, mode='stochastic') 
scv.tl.velocity_graph(adata) 

scv.pl.velocity_embedding(adata, basis='umap_harmony', frameon=False, save='embedding.pdf')

# visualization of RNA velocity
scv.pl.velocity_embedding_grid(adata, basis='umap_harmony', color='celltype', 
    save='embedding_grid.pdf', title='', scale=0.25,figsize =(6,7))
scv.pl.velocity_embedding_stream(adata, basis='umap_harmony', color='celltype', 
    save='embedding_stream.pdf', title='',figsize =(6,7))