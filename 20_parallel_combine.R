#!/bin/env Rscript
#SBATCH --mem=200G 
#SBATCH --time=24:00:00 
#SBATCH --partition=BioCompute,Lewis,hpc5
#SBATCH --account=warrenlab
#SBATCH --cpus-per-task=12
#SBATCH -o job_files.dir/parallel_combined_%j_o.txt
library(cowplot)
library(Seurat)

# Parallel magic
library(future)
plan(multiprocess)
options(future.globals.maxSize = 200000 * 1024^2) # (this is 200G, default is like 512MiB (512 * 1024^2))

pdf_height <- 6     # inches
obj_dir    <- 'objects'
plot_dir   <- 'plots'

paste_dir <- function(x) {
    return(paste(x,collapse='/'))
}

simple_print <- function(x) {
    print(paste0(x))
}

samp_names <- list.files('count')

#=====================================================
# PLOTS
#-------------------------------------------------------

plot_pdf <- function(name, plot) {
    filename <- paste0(name,'.pdf')
    pdf(filename,width=pdf_height,height=pdf_height)
        print(plot)
    dev.off()
    return(filename)
}

create_plots <- function(name, plot) {
    pdf_filename <- plot_pdf(name, plot)
    png_filename <- paste0(name,'.png')
    system(paste0('convert -density 600 ',pdf_filename, ' ', png_filename))
}

#-------------------------------------------------------
# PLOTS
#=====================================================

ref <- 'pre_mRNA'


integrated_data_RDS_filename <- paste_dir(c(obj_dir, paste0('integrated_all_samples_',ref,'.rds')))
integrated_data <- c()
if (file.exists(integrated_data_RDS_filename)) {
    simple_print(c('loading ',integrated_data_RDS_filename))
    integrated_data <- readRDS(integrated_data_RDS_filename)
} else {
    vfeatures <- list()
    for (samp_name in samp_names) {
        vfeature_RDS_filename=paste_dir(c(obj_dir,c(paste0(samp_name,'_vfeature.rds'))))
        simple_print(c('loading ',vfeature_RDS_filename))
        vfeatures[[samp_name]] = readRDS(vfeature_RDS_filename)
    }

    anchors         <- FindIntegrationAnchors(object.list=vfeatures, dims=1:20)
    integrated_data <- IntegrateData(anchorset=anchors, dims=1:20)
    simple_print(c('saving ',integrated_data_RDS_filename))
    saveRDS(integrated_data,file=integrated_data_RDS_filename)
}

clusters_RDS_filename <- paste_dir(c(obj_dir,paste0('clusters_', ref, '.rds')))
clusters <- c()
if (file.exists(clusters_RDS_filename)) {
    simple_print(c('loading ',clusters_RDS_filename))
    clusters <- readRDS(clusters_RDS_filename)
} else {
    
    scaled_data <- ScaleData(integrated_data, verbose=FALSE)
    pcas        <- RunPCA(scaled_data, npcs=30, verbose=FALSE)
    ump         <- RunUMAP(pcas, reduction="pca", dims=1:20)
    neighbors   <- FindNeighbors(ump, reduction="pca", dims=1:20)
    clusters    <- FindClusters(neighbors, resolution=0.5)
    simple_print(c('saving ',clusters_RDS_filename))
    saveRDS(clusters,clusters_RDS_filename)
}


for (label_bool in c(TRUE,FALSE)) {
    labelling <- ''
    if (! label_bool) {
        labelling <- 'unlabeled_'
    }
    subgroup_comparison <- paste_dir(c(plot_dir,paste0(labelling,'Subgroup_comparison_type_',ref)))

    create_plots(subgroup_comparison,DimPlot(clusters,reduction="umap",label=label_bool))
    
    subgroup_comparison_split <- paste_dir(c(plot_dir,paste0(labelling,'Subgroup_comparison_type_',ref,'split_by_sample')))

    create_plots(subgroup_comparison_split,DimPlot(clusters, reduction="umap",split.by="sample_type",label=label_bool)) # was split.by="sample_sub_group" before
    
    p1 <- DimPlot(clusters, reduction="umap", group.by="sample_type", label=label_bool)
    p2 <- DimPlot(clusters, reduction="umap", label=label_bool)
    
    grid_compare <- paste_dir(c(plot_dir,paste0(labelling,'Grid_comparison_',ref)))

    create_plots(grid_compare,plot_grid(p1,p2))
}
