### construct gene regulatory network using scMEGA ###


# read the snATAC-seq data 
 ATAC100 <- readRDS('../Save-ArchR-Project.rds')

# run UMAP and calculate the average frag of each cluster
p1<-plotEmbedding(ArchRProj=ATAC100,colorBy="cellColData",name="Clusters",
	embedding="UMAP",plotAs="points",size=0.5,labelAsFactors=FALSE,rastr=FALSE)+theme_cowplot()

df <-ATAC100@cellColData %>%as.data.frame()%>%group_by(Clusters)%>%summarise(mean_nFrags=mean(nFrags))
p1<-plotGroups(ArchRProj=ATAC100,groupBy="Clusters",colorBy="cellColData",name="TSSEnrichment",plotAs="violin",alpha=0.4,addBoxPlot=TRUE)
pdf(file="violin_ATAC100_clu_TSS.pdf", width = 8, height = 6)
p1
dev.off()

# read the original peak matrix, only check the original total number of cells
counts <- Read10X_h5("../filtered_peak_bc_matrix.h5")
dim(counts)
head(counts)[,1:5]

# extract the peak matrix from the ArchR object
getAvailableMatrices(ATAC100)
atac <- getMatrixFromProject(ArchRProj = ATAC100,useMatrix = "PeakMatrix")
peak_counts <- atac@assays@data$PeakMatrix
rownames(peak_counts) <- atac@elementMetadata$name
saveRDS(peak_counts, file = "../PeakMatrix.Rds")

peak_counts <- readRDS(glue::glue("../PeakMatrix.Rds"))
dim(peak_counts)

row <- atac@ rowRanges
row <-  as.data.frame(row)
row$use <-  paste(row$seqnames,":",sep="")
row$use <-  paste(row$use,row$start,sep="")
row$use <-  paste(row$use,"-",sep="")
row$use <-  paste(row$use,row$end,sep="")
head(row)[1:2,]
dim(row)
rownames(peak_counts) <- row$use
head(peak_counts)[1:2,1:3]

# save the peak matrix
saveRDS(peak_counts, file = "../peak_counts.Rds"


# read the metadata
metadata <- as.data.frame(ATAC100@cellColData)
dim(metadata)
metadata <- as.data.frame(ATAC100@cellColData) %%subset(., select = c("Sample","Clusters","Celltype"))


# peak matrix was extracted according to the selected cell names
peak_counts_sub <- peak_counts[,colnames(peak_counts) %in% rownames(metadata)]
dim(peak_counts_sub)
rownames(peak_counts_sub) <- paste("chr",rownames(peak_counts_sub),sep="")
head(peak_counts_sub)[,1:2]                       .                            .

# creat snATAC object
chrom_assay <- CreateChromatinAssay(counts = peak_counts_sub,sep = c(":", "-"),min.cells = 100)
obj.atac<-CreateSeuratObject(counts=chrom_assay,assay="ATAC",meta.data=metadata,names.field=1,names.delim="#")
obj.atac

## run UMAP
embedding <- ATAC100@embeddings$UMAP$df
colnames(embedding) <- paste0("UMAP_", 1:ncol(embedding))
embedding <- embedding[colnames(obj.atac), ]
obj.atac[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding),assay = "ATAC",key ="UMAP_")
obj.atac

# UMAP plot
cols.clusters <- ArchR::paletteDiscrete(metadata$Clusters)
p1 <- DimPlot(obj.atac, group.by = "Clusters", pt.size = 1, cols = cols.clusters,reduction = "umap")
pdf(file="UMAP_obj.atac_clu.pdf", width = 8, height = 6)
p1
dev.off()

# read the snRNA-seq data
obj.rna <- readRDS("/home/../horn.rna.rds")
obj.rna

# save the objects
saveRDS(obj.atac, "snATAC.rds")
saveRDS(obj.rna, "snRNA.rds")

# save the gene score matrix
atac <- getMatrixFromProject(ArchRProj = ATAC100,useMatrix = "GeneScoreMatrix")
gene_counts <- atac@assays@data$GeneScoreMatrix
rownames(gene_counts) <- atac@elementMetadata$name
saveRDS(gene_counts, file = "../GeneScoreMatrix.Rds")
gene.activity <- readRDS(glue::glue("../GeneScoreMatrix.Rds"))
gene.activity <- gene.activity[, colnames(obj.atac)]
dim(gene.activity)


# combine snRNA-seq and sATAC-seq
obj.coembed<-CoembedData(obj.rna,obj.atac, gene.activity, weight.reduction ="umap", dims = 1:2,verbose =TRUE)

# UMAP plot by tech
p2 <- DimPlot(obj.coembed, group.by = "tech", shuffle = TRUE, label = TRUE)
pdf(file="UMAP_co_tech.pdf", width = 8, height = 6)
p2
dev.off()

# reUMAP by harmony
obj.coembed<-RunHarmony(obj.coembed,group.by.vars =c( "tech"),reduction ="pca",
	max.iter.harmony =30,dims.use =1:30,project.dim =FALSE,plot_convergence =FALSE)
obj.coembed<-RunUMAP(obj.coembed,dims =1:30,reduction ='harmony',reduction.name ="umap_harmony",
	reduction.ke ='umapharmony_',verbose =FALSE,min.dist =0.4)
p2<-DimPlot(obj.coembed, group.by ="tech", reduction ="umap_harmony")
pdf(file="UMAP_co_har_tech.pdf", width = 8, height = 6)
p2
dev.off()

# add cluster
obj.coembed <- FindNeighbors(obj.coembed, reduction = "harmony", dims = 1:30)
obj.coembed <- FindClusters(obj.coembed, resolution = 1, verbose = FALSE)
cols <- ArchR::paletteDiscrete(obj.coembed@meta.data[, "RNA_snn_res.1"])
p <- DimPlot(obj.coembed, group.by = "RNA_snn_res.1", label = TRUE,
	reduction = "umap_harmony", shuffle = TRUE)+ scale_color_manual(values = cols) +xlab("UMAP1") + ylab("UMAP2")
pdf(file="UMAP_co_clu.pdf", width = 8, height = 6)
p
dev.off()

# statistical cell proportion
p2 <- CellPropPlot(obj.coembed,group.by = "tech",prop.in = "RNA_snn_res.1")
pdf(file="prop_co_tech.pdf", width = 8, height = 6)
p2
dev.off()

# top3 markers plot 
all.markers <- FindAllMarkers(obj.coembed,only.pos = TRUE,min.pct = 0.5, logfc.threshold = 0.5)
df<-all.markers%%group_by(cluster)%%slice_max(n =3, order_by =avg_log2FC)
p<-DotPlot(obj.coembed, features =unique(df$gene))+RotatedAxis()
pdf(file="marker_top3_clu.pdf", width = 18, height = 9)
p
dev.off()

# save the objects
saveRDS(obj.coembed, "./coembed.rds")


Idents(obj.coembed)<-"RNA_snn_res.1"
coembed.sub<-subset(obj.coembed, idents =c(5, 11, 12, 19, 29, 32, 35), invert =TRUE)
coembed.sub

options(repr.plot.height =6, repr.plot.width =6)
p<-DimPlot(coembed.sub, group.by ="RNA_snn_res.1", label =TRUE,reduction ="umap_harmony",
	shuffle =TRUE, cols =cols.clusters)+xlab("UMAP1")+ylab("UMAP2")
pdf(file="UAMP_select_clu.pdf", width = 8, height = 6)
p
dev.off()

# re-perform UMAP embedding and clustering at a lower resolution
coembed.sub<-RunUMAP(coembed.sub, dims =1:30, reduction ='harmony',reduction.name ="umap_harmony",
	reduction.key ='umap_harmony_',verbose =FALSE,min.dist =0.4)
coembed.sub <- FindNeighbors(coembed.sub, reduction = "harmony", dims = 1:30)
coembed.sub <- FindClusters(coembed.sub, resolution = 0.1, verbose = FALSE)
cols <- ArchR::paletteDiscrete(coembed.sub@meta.data[, "RNA_snn_res.0.1"])
p <- DimPlot(coembed.sub, group.by = "RNA_snn_res.0.1", label = TRUE,reduction = "umap_harmony",
	shuffle = TRUE, cols = cols) +xlab("UMAP1") + ylab("UMAP2")
pdf(file="UAMP_select_clu0.1.pdf", width = 8, height = 6)
p
dev.off()

# find and plot the markers
all.markers <- FindAllMarkers(coembed.sub,only.pos = TRUE,min.pct = 0.25, logfc.threshold = 0.5)
df<-all.markers%%group_by(cluster)%%slice_max(n =10, order_by =avg_log2FC)
p<-DotPlot(coembed.sub, features =unique(df$gene))+RotatedAxis()
pdf(file="marker_top10_subclu0.1.pdf", width = 18, height = 6)
p
dev.off()

# display the 
p<-DimPlot(coembed.sub, group.by ="RNA_snn_res.0.1", label =TRUE,reduction ="umap_harmony", 
	shuffle =TRUE, split.by ="tech",cols =cols)+xlab("UMAP1")+ylab("UMAP2")
pdf(file="UMAP_select_clu0.1_tech.pdf", width = 12, height = 6)
p
dev.off()

# save the object
saveRDS(coembed.sub, "./coembed.sub.rds")


# pair the cells
df.pair <- PairCells(object = coembed.sub, reduction = "harmony",pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")

# visualize the pair results
sel_cells<-c(df.pair$ATAC, df.pair$RNA)
coembed.sub2<-coembed.sub[, sel_cells]
coembed.sub2

options(repr.plot.height =5, repr.plot.width =10)
p <- DimPlot(coembed.sub2, reduction ="umap_harmony", group.by ="RNA_snn_res.0.1", split.by ="tech", cols =cols)
pdf(file="UMAP_pair_tech.pdf", width = 12, height = 6)
p
dev.off()
obj.pair<-CreatePairedObject(df.pair =df.pair, object =coembed.sub2,use.assay1 ="RNA", use.assay2 ="ATAC")
obj.pair

# marker genes plots
p <- FeaturePlot(obj.pair, features = c("COL1A1","COL3A1","PDGFRA","KRT15","KRT17","KRT14","PECAM1","KDR","GRIK2","PPP2R2B","ACTA2","PTPRC","PECAM1","KDR","PROX1","MMRN1","MYH3","TYR","KIT","GATA2"),reduction ="umap_harmony",min.cutoff = "q10", max.cutoff = "q90")
pdf(file="UMAP_cosub_RNAmarker.pdf", width =28, height = 26)
p
dev.off()

# identify the cell types using previous snRNA-seq results
p2 = DimPlot(obj.pair, group.by="celltype", label=T, label.size=5, reduction='umap_harmony')
pdf(file="newUMAP_cosub_celltype.pdf", width =8, height = 6)
p2
dev.off()

# add the trajectory
obj.pair<-AddTrajectory(object =obj.pair, trajectory =c(22,1,7,0,2,11),
	group.by ="RNA_snn_res.0.5", reduction ="umap_harmony",dims =1:2, use.all =TRUE)
obj <- obj.pair[, !is.na(obj.pair$Trajectory)]
p1 <- DimPlot(obj, reduction = "umap_harmony",group.by = "RNA_snn_res.0.5", cols = cols) 
	+xlab("UMAP 1") + ylab("UMAP 2") +ggtitle("Cluster")
p2 <- TrajectoryPlot(object = obj,reduction = "umap_harmony",continuousSet = "blueYellow",
	size = 1,addArrow = FALSE) +xlab("UMAP 1") + ylab("UMAP 2") +ggtitle("Trajectory")
pdf(file="traj_test_FIB.pdf", width =10, height = 6)
p1+p2
dev.off()

pfm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
obj <- AddMotifs(object = obj,genome = BSgenome.Btaurus.UCSC.bosTau9,pfm =pfm,assay ="ATAC")
obj<-RunChromVAR(object =obj,genome =BSgenome.Btaurus.UCSC.bosTau9,assay ="ATAC")

# select TFs
res<-SelectTFs(object =obj, return.heatmap =TRUE,cor.cutoff =0.1)
df.cor<-res$tfs
ht<-res$heatmap
pdf(file="traj_heatmap_FIB_TF2gene.pdf", width =10, height = 6)
ht
dev.off()

# select genes
res <- SelectGenes(object = obj,labelTop1 = 0,labelTop2 = 0)
df.p2g <- res$p2g
ht <- res$heatmap
pdf(file="traj_heatmap_FIB_peak2gene.pdf", width =10, height = 6)
ht
dev.off()

# construc the GRN
tf.gene.cor <- GetTFGeneCorrelation(object = obj,tf.use = df.cor$tfs,gene.use = unique(df.p2g$gene),
	tf.assay = "chromvar",gene.assay ="RNA",trajectory.name = "Trajectory")
ht <- GRNHeatmap(tf.gene.cor,tf.timepoint = df.cor$time_point)
pdf(file="GRN_heatmap_FIB.pdf", width =16, height = 6)
ht
dev.off()
motif.matching <- obj@assays$ATAC@motifs@data
colnames(motif.matching) <- obj@assays$ATAC@motifs@motif.names
motif.matching <-motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]
df.grn <- GetGRN(motif.matching = motif.matching,df.cor = tf.gene.cor,df.p2g = df.p2g)

# visualize the GRN
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs
df.grn2 <- df.grn %%subset(correlation  0.5) %%select(c(tf, gene, correlation)) %%rename(weights = correlation)
p <- GRNPlot(df.grn2,tfs.timepoint = tfs.timepoint,show.tf.labels = TRUE,seed = 42,
	plot.importance = FALSE,min.importance =remove.isolated = FALSE)
options(repr.plot.height = 20, repr.plot.width = 20)
pdf(file="GRN_FIB.pdf", width =16, height = 6)
p
dev.off()