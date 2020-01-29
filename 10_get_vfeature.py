#!/bin/env python
import argparse
import os

debugging = "TRUE" 

if (debugging == "TRUE"):
    os.system("mkdir -p debug_obj")

def main(args):

    sample_dirs = os.listdir(args.count_dir) 

    script_dir = 'scripts'
    object_dir = 'objects'
    plot_dir = 'plots'
    job_files_dir = 'job_files.dir'
    marker_dir = 'markers'

    for i in (script_dir, object_dir, plot_dir, job_files_dir, marker_dir):
        os.system('mkdir -p ' + i)

    for sample_dir in sample_dirs:

        # Incorporate cores and mem requests
        script = """#!/bin/env Rscript
#SBATCH --mem={mem}G 
#SBATCH --time=2-00:00:00 
#SBATCH --partition={partition_string}
#SBATCH --cpus-per-task={cores}
#SBATCH --account=warrenlab
#SBATCH --job-name=create_vfeatures_for{sample_dir}
#SBATCH -o job_files.dir/create_vfeatures_for{sample_dir}_%j_o.txt
debugging <- {debugging}
library(Seurat)
library(future)
plan(multiprocess, workers={cores})

min_cells    <- 3
min_nFeature <- 200
max_nFeature <- 2500
max_mito     <- 5    
max_ribo     <- 30

options(future.globals.maxSize = {mem}000 * 1024^2) # default is like 512MiB (512 * 1024^2)
sample_dir <- '{sample_dir}'
""".format(mem=args.mem,
           partition_string=args.partition_string,
           debugging=debugging,
           cores=args.cores,
           sample_dir=sample_dir)

        script += """
paste_dir <- function(x) {
    return(paste(x,collapse='/'))
}
"""

        script += """
full_dir_name = paste_dir(c('count',sample_dir,'{bc_matrix_dir}'))
""".format(bc_matrix_dir=args.bc_matrix_dir)

        script += """

simple_paste <- function(x) {
    return(paste(x,collapse=''))
}

debug_print <- function(x) {
    if (debugging) {
        print(simple_paste(x))
    }
}

# save an RDS file (and mention it if debugging is on)
save_file <- function(x,y) {
    filename <- simple_paste(y)
    debug_print(c('saving ',filename))
    saveRDS(x,filename)
    return()
}

debug_save <- function(x,y) {
    if (debugging) {
        filename <- paste_dir(c("debug_obj",simple_paste(y)))
        debug_print(c('saving ',filename))
        saveRDS(x,filename)
    }
    return()
}

# WARNING: Getting sample_type is kludgey. Use json file in future version.
# Extract "Surface" or "Cave" from sample_dir
sample_type <- gsub('\\\\d\\\\d(C|M)','\\\\1',sample_dir,perl=TRUE)
abbreviated_sample_name <- sample_dir
debug_print(c('Extracted sample_type from sample name: ', sample_type))

debug_print(c('Reading in ', full_dir_name))
my_counts <- Read10X(data.dir = full_dir_name)
debug_save(my_counts,simple_paste(c(sample_dir,"__my_counts.rds")))
debug_print(c('Filter for genes present in at least three cells (', sample_dir,')'))
samp_obj <- CreateSeuratObject(counts = my_counts, min.cells = 3)
samp_obj$sample_type <- sample_type
samp_obj$abbreviated_sample_name <- abbreviated_sample_name
debug_save(samp_obj,simple_paste(c(sample_dir,"samp_obj.rds")))

debug_print(c('Filtering for at least 200 genes but fewer than 2500 and for less than 5% mitochondrial genes'))
samp_obj[["percent.mt"]] <- PercentageFeatureSet(samp_obj, pattern = "^MT-")
samp_obj[["percent.ribo"]] <- PercentageFeatureSet(samp_obj, pattern = "^RP[LS]")

vfeature <- SCTransform(samp_obj, vars.to.regress = c('percent.mt','percent.ribo'), verbose=FALSE)

debug_print(c(sample_dir, ': Running PCA'))
pca_obj <- RunPCA(vfeature)
pca_obj$sample_type <- sample_type 
"""

        script +="""
object_file_name = paste_dir(c("{object_dir}",simple_paste(c(abbreviated_sample_name,'_vfeature.rds'))))
save_file(pca_obj, object_file_name)
""".format(object_dir=object_dir)
        script_name = os.path.join(script_dir, sample_dir + '.get_vfeature.R')
        with open(script_name,"w") as out_file:
            out_file.write(script)
        command = args.load_modules + '; sbatch ' + script_name
        os.system(command)

# command line interface (making this a modulino)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download a directory tree from an FTP server')
    parser.add_argument('--cores',default=12,help='number of cores to request')
    parser.add_argument('--mem',default=80,help='memory to request (in GB)')
    parser.add_argument('--load_modules',
                         default="module load biocompute/biocompute-modules; module load seurat/seurat-3.1.0_R-3.6.1",
                         help='module loading commands')
    parser.add_argument('--partition_string',default='BioCompute,Lewis,hpc5',help='partition(s) to request')
    parser.add_argument('--bc_matrix_dir',default='outs/filtered_feature_bc_matrix',help='directory containing barcode matrix files (output from cell ranger)')
    parser.add_argument('count_dir',help='Directory containing the results of cellranger count run on each sample directory')

    args = parser.parse_args()
    main(args)
