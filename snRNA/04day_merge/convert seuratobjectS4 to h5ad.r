### seuratobject S4 to h5ad ###

library(Seurat)
library(SeuratData)
library(SeuratDisk)
InstallData("neu")
data("neu.final")
neu.final
SaveH5Seurat(neu.final, filename = "neu.h5Seurat")
Convert("neu.h5Seurat", dest = "h5ad")