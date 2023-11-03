## 需要提前准备文件以及软件安装   
1. 物种重复序列注释文件repeat_masker.gtf，可以直接在ucsc上下载（如果没有的话，可以用RepeatMasker自行生成，此文件非必需）；
2. 物种基因组注释文件gtf；
3. 每个scRNA样本cellranger后的结果文件，主要是需要bam文件；
4. 通过seurat包分析注释完成的rds文件；
5. 软件：velocyto，scVelo；
6. R包：Seurat，SeuratData，SeuratDisk