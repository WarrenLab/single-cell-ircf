#!/bin/env python3
#SBATCH --account=warrenlab
#SBATCH --partition=BioCompute,hpc5,Lewis
#SBATCH --mem=40G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=SAMPLE_ID_extract_spreadsheet
#SBATCH -o SAMPLE_ID_extract_spreadsheet_%j_o.txt

import os
import argparse
import json

DEBUG = False

def value_or_default_for(my_dictionary,my_key,default):
    if my_key in my_dictionary:
        return my_dictionary[my_key]
    else:
        return default 

def dictionary_slice(my_dictionary,desired_keys,default):
    values = [value_or_default_for(my_dictionary,k,default) for k in desired_keys]
    return values

def main(args):
    DEBUG = args.debug
    
    counts = {}
    gene_list = []
    barcode_list = []
    gene_id_for = {}
    gene_symbol_for = {}

    with open(args.out_tsv, "w") as fh_out:
        with open(args.flat_csv) as fh_in:
            for line in fh_in:
                barcode,gene_id,gene_symbol,skip_me,gene_count = line.rstrip().split(",")   
                gene_symbol_id_combo = gene_symbol + "(" + gene_id + ")"

                if DEBUG:
                    print("line: " + line)
                    print(f"Barcode/gene_symbol/count: {barcode}/{gene_symbol}/{gene_count}")

                if gene_symbol_id_combo not in counts:
                    counts[gene_symbol_id_combo] = {}
                    gene_list.append(gene_symbol_id_combo)

                if barcode in counts[gene_symbol_id_combo]:
                    print(f"ERROR: Duplicate combination of barcode {barcode}, gene symbol {gene_symbol}, and gene id {gene_id}") 
                    exit(1)
                else:
                    counts[gene_symbol_id_combo][barcode] = gene_count
                    gene_id_for[gene_symbol_id_combo] = gene_id
                    gene_symbol_for[gene_symbol_id_combo] = gene_symbol
                    if barcode not in barcode_list:
                        barcode_list.append(barcode)

        gene_list.sort()
        barcode_list.sort()

        # Print header (i.e. first column is "gene_symbol", remaining columns are barcodes)
        fh_out.write("\t".join(["gene_symbol_id_combo","gene_symbol","gene_id","total_counts","\t".join(barcode_list)]) + "\n")

        for gene_symbol_id_combo in gene_list:

            if DEBUG:
                print(counts[gene_symbol_id_combo])

            # Get gene counts for all barcodes for this gene
            values = dictionary_slice(
                         counts[gene_symbol_id_combo],   # dictionary containg all the counts for this gene_symbol_id_combo
                         barcode_list,                   # list of all barcodes (previously sorted)
                         "0"                             # Default of "0" for missing values
                     )
            total_gene_counts = sum(values)
            # Print out gene identifcation info and gene counts for each barcode
            fh_out.write(
                "\t".join([
                        gene_symbol_id_combo,
                        gene_symbol_for[gene_symbol_id_combo],
                        gene_id_for[gene_symbol_id_combo],
                        total_gene_counts,
                        "\t".join(values)
                    ])
                + "\n"
            )

# command line interface (making this a modulino)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                description='Produce a large "dense" tab-separated-value file from a "flat" csv file which was generated from barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz'
             )
    parser.add_argument(
        '--debug',
        dest='debug',
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--out_tsv',
        type=str,
        default='final_spreadsheet.tsv',
        help='Output file containing all of the barcodes, genes, and their counts, with barcodes as columns, genes as rows, and gene counts as cell values',
    )
    parser.add_argument(
        'flat_csv',
        type=str,
        default='final_matrix.csv',
        help='File containing all of the barcodes, genes, and their counts, with one line per barcode/gene combination',
    )

    args = parser.parse_args()

    main(args)


