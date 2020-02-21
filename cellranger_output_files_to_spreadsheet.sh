#!/bin/bash
#SBATCH --account=warrenlab
#SBATCH --partition=BioCompute,hpc5,Lewis
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=SAMPLE_ID_combine_raw_cellranger_data
#SBATCH -o SAMPLE_ID_combine_raw_cellranger_data_%j_o.txt

# Print line number along with contents of barcodes.tsv.gz and genes.tsv.gz 
pigz -d -c barcodes.tsv.gz | awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1}' | sort -t, -k 1b,1 > numbered_barcodes.csv
pigz -d -c features.tsv.gz | awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1,$2,$3}' | sort -t, -k 1b,1 > numbered_features.csv

# Skip the header lines and sort matrix.mtx.gz
pigz -d -c matrix.mtx.gz | tail -n +4 | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 1b,1 > feature_sorted_matrix.csv
pigz -d -c matrix.mtx.gz | tail -n +4 | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 2b,2 > barcode_sorted_matrix.csv

# Use join to replace line number with barcodes and genes
join -t, -1 1 -2 1 numbered_features.csv feature_sorted_matrix.csv | cut -d, -f 2,3,4,5,6 | sort -t, -k 4b,4 | join -t, -1 1 -2 4 numbered_barcodes.csv - | cut -d, -f 2,3,4,5,6 > final_matrix.csv

# Remove temp files
rm -f barcode_sorted_matrix.csv feature_sorted_matrix.csv numbered_barcodes.csv numbered_features.csv
python3 flat_to_table.py final_matrix.csv
