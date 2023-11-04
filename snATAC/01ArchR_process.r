### process snATAC-seq using ArchR ###


# load packages and reference genome 
library('ArchR')
library('Seurat')
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Btaurus.UCSC.bosTau9)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene, OrgDb = org.Bt.eg.db)

# set the threads
addArchRThreads(threads = 8)  

# read the data
inputFiles <- list.files("fragments.tsv.gz", ##results of cellranger 
	pattern = ".gz",full.names = TRUE)

# creat arrowfiles
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,sampleNames = names(inputFiles),minTSS= 4,minFrags= 1000,addTileMat = TRUE,addGeneScoreMat = TRUE,geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation) 

# doublets inference
doubScores <- addDoubletScores(input = ArrowFiles, k = 10,knnMethod = "UMAP",LSIMethod = 1)   

# creat a ArchRproject
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "", copyArrows = TRUE,geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation) 

# plot QC metrics 
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(x = df[,1],y = df[,2],colorDensity = TRUE,continuousSet = "sambaNight",xlabel = "Log10 Unique Fragments",ylabel = "TSS Enrichment",xlim = c(log10(500), quantile(df[,1], probs = 0.99)),ylim = c(0, quantile(df[,2], probs = 0.99))) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
pdf("TSS-vs-Frags.pdf")   
p
dev.off()

p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample",colorBy = "cellColData", name = "TSSEnrichment",plotAs = "ridges")
p2 <- plotGroups(ArchRProj = proj,groupBy = "Sample",colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin", alpha = 0.4,addBoxPlot = TRUE)  
pdf("TSSEnrichment.pdf")
p1+p2
dev.off()

p3 <- plotGroups(ArchRProj = proj, groupBy = "Sample",colorBy = "cellColData", name = "log10(nFrags)",plotAs = "ridges")
p4 <- plotGroups(ArchRProj = proj,groupBy = "Sample",colorBy = "cellColData", name = "log10(nFrags)",plotAs = "violin", alpha = 0.4,addBoxPlot = TRUE)
pdf("UniqueUuclearFragments.pdf")
p3+p4
dev.off()

# save the ArchRproject
saveArchRProject(ArchRProj = proj, outputDirectory = "", load = FALSE)

# filter the doublets
proj <- filterDoublets(proj)

# dimensionality reduction by running the iterative LSI
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# clustering using Seurat’s FindClusters() function
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
table(projhorn2$Clusters)

# create a cluster confusion matrix across each sample
cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"),  border_color = "black")
pdf("pheatmap_fivedata")
p
dev.off()

# UMAP embedding 
proj <- addUMAP(ArchRProj =proj, reducedDims = "IterativeLSI")

# visualize the UMAP result
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", highlightCells=result,embedding = "UMAP",plotAs ="points")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj =proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",ArchRProj =proj, addDOC = FALSE, width = 5, height = 5)

# identify marker genes
markersGS <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix",groupBy = "Clusters",  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <=0.05& Log2FC >=1") 
markerList$C6

# display the interesting marker genes
markerGenes<-c(" "," ")
heatmapGS <- plotMarkerHeatmap(seMarker = markersGS,  labelMarkers = markerGenes,transpose = TRUE)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj =proj, addDOC = FALSE)

# visualize the marker genes on UMAP
p <- plotEmbedding( ArchRProj = proj,  colorBy = "GeneScoreMatrix", name = markerGenes,  embedding = "UMAP", quantCut = c(0.01, 0.95), imputeWeights = NULL)
proj <- addImputeWeights(proj,reducedDims = "IterativeLSI")
p <- plotEmbedding( ArchRProj = proj, colorBy = "GeneScoreMatrix",   name = markerGenes, embedding = "UMAP",imputeWeights = getImputeWeights(proj))
plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",  ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# visualizing Genome Browser Tracks
p <- plotBrowserTrack(ArchRProj = proj, groupBy = "Clusters", geneSymbol = markerGenes, upstream = 50000, downstream = 50000)
grid::grid.newpage()
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# save the ArchRProject
proj <- saveArchRProject(ArchRProj = proj)
proj <- loadArchRProject(path = "projaggrAl")

# integrated snRNA-seq
seRNA <- readRDS(".rds")
seRNA <- as.SingleCellExperiment(seRNA) ##conver single cell data format
colnames(colData(seRNA))
table(colData(seRNA)$celltype) 

# unconstrained integration
proj <- addGeneIntegrationMatrix(ArchRProj = proj, useMatrix = "GeneScoreMatrix",matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA,addToArrow = FALSE,groupRNA = "celltype",nameCell = "predictedCell_Un", nameGroup = "predictedGroup_Un",nameScore = "predictedScore_Un")

# visualize the integration results
pal <- paletteDiscrete(values = colData(seRNA)$celltype)
p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup_Un", pal = pal)
plotPDF(p1, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# use snRNA-seq information to identify snATAC-seq
cM <- confusionMatrix(projhorn5$Clusters, projhorn5$predictedGroup_Un)
labelOld <- rownames(cM)
labelOld
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
labelNew2
proj$Clusters2 <- mapLabels(proj$Clusters, newLabels = labelNew, oldLabels = labelOld)
p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters2")
p1
plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = proj addDOC = FALSE, width = 5, height = 5)

# generation of pseudo-bulk replicates
proj <- addGroupCoverages(ArchRProj = proj,groupBy = "predictedGroup_Un/Clusters2/Clusters/Sample")

# peak-calling
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "", pathToMacs2 = pathToMacs2,genomeSize = 2.7e+09)
getPeakSet(proj)
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)

table(proj$Clusters2)
markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "",bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markerlist <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

# plot markerpeak heatmap
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj addDOC = FALSE)

# markerpeak MA and Volcano Plots
pma <- markerPlot(seMarker = markersPeaks, name = "Neuron", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma
pv <- markerPlot(seMarker = markersPeaks, name = "Neuron", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pv
plotPDF(pma, pv, name = "Neuron-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

# markerpeak in Browser Tracks
p <- plotBrowserTrack(ArchRProj = proj, groupBy = "", geneSymbol = c("KLK9"),features = getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Keratinocytes"],upstream = 50000,downstream = 50000)
grid::grid.newpage()
grid::grid.draw(p$KLK9)
plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
