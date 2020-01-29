#!/bin/env Rscript
#SBATCH --mem=40G 
#SBATCH --time=24:00:00 
#SBATCH --partition=BioCompute,Lewis,hpc5
#SBATCH --cpus-per-task=12
#SBATCH --account=warrenlab
#SBATCH -o job_files.dir/30_find_markers_%j_o.txt

clusters_RDS_filename <- 'objects/clusters_pre_mRNA.rds'

library(Seurat)

# multiprocessor magic (use 12 cores, 40GB RAM)
library(future)
plan("multiprocess", workers=12)
options(future.globals.maxSize = 40000 * 1024^2) # default is like 512MiB (512 * 1024^2)

library(cowplot)

clusters <- readRDS(clusters_RDS_filename)

all_cell_markers <- FindAllMarkers(clusters, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(all_cell_markers, file='markers/FindAllMarkers.tsv', sep="\t", quote=FALSE)
