# single-cell-ircf
Single cell scripts from MU's IRCF

# Why another repo for single cell scripts?
Since these are so different from WarrenLab/single-cell, I've created a separate repository for them. We can work towards a unified approach, if desired, or simply provide multiple options.

# Note on previously installed Seurat on Lewis

This approach uses Seurat as already installed on Lewis. You can access the same version via these two module load commands:

    module load biocompute/biocompute-modules
    module load seurat/seurat-3.1.0_R-3.6.1

# Pipeline synopsis 

    # Load the required modules and run 10_get_vfeature.py on the "count" directory.
    # (Note: currently need to configure a regex in 10_get_vfeature.py to correctly identify sample types)
    run_10_get_vfeature.sh

    # Call clusters on combined data and create plots
    seurat_run 20_parallel_combine.R

    # Find and output markers 
    seurat_run 30_find_markers.R

    # Run 40_top_genes_R_template, which generates the heatmap figures
    run_40_generate_heatmap.sh

# Utilities
    seurat_run # load seurat module and then run the specified R script

    clean_up_here # Delete pipeline generated directories (useful after having made changes to a pipeline script)

    sc_count_normal_3.1.0 # Run cellranger count on fastq data

    tenX_cellranger_3.1.0_mkfastq # Convert BCL files to FASTQ files (via cellranger mkfastq)


# TODO: Use a config file (thus making it more accessible to general users)

    1. Use a config file (instead of detecting sample types by updating the pattern matching line in 10_get_vfeature.py.)
        a. List samples with corresponding sample type
        b. List SLURM account to use for data processing (currently hard-coded to "warrenlab" for this repository)
        c. Generalize mitochondrial and ribosomal filtering (or create mito-ribo free references). One issue with this is that mitochondrial and ribosomal genes are not consistently named. In chicken, for example, I manually had to prefix mitochondrial genes with "MT-" (in the "features" file inside "bc_matrix" directories).
