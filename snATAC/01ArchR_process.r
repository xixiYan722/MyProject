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

# 
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "输出文件的目录",copyArrows = TRUE,geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation)     
