#!/bin/env perl
#SBATCH -J flat_to_table
#SBATCH -o flat_to_table_%j_o.txt
#SBATCH --partition BioCompute,Lewis,hpc5
#SBATCH --mem 240G
#SBATCH --ntasks 1
#SBATCH --nodes 1
use strict;
use warnings;
use v5.10;
use Memoize;

memoize 'gene_symbol_and_id';

use List::Util 'sum0';

my %counts_for;
my %barcodes_seen;
my %gene_symbol_and_ids_seen;
my %sum_for_barcode;
my %sum_for_gene;

while(<>) {
    chomp;
    my ($barcode, $gene_id, $gene_symbol, undef, $gene_count) = split /,/;

    my $gene_symbol_and_id = gene_symbol_and_id($gene_symbol, $gene_id);

    $counts_for{$gene_symbol_and_id}{$barcode} = $gene_count;

    $sum_for_gene{$gene_symbol_and_id} += $gene_count;
    $sum_for_barcode{$barcode}         += $gene_count;

    # Record any new barcodes or gene_symbol_and_ids
    $barcodes_seen{$barcode}                       ||= 1;
    $gene_symbol_and_ids_seen{$gene_symbol_and_id} ||=1;
}

my @gene_list = sort keys %gene_symbol_and_ids_seen;
my @barcode_list = sort keys %barcodes_seen;

say (join("\t","gene_symbol_and_id","total_counts", @barcode_list));

for my $gene_symbol_and_id ( @gene_list ) {

    print "$gene_symbol_and_id";

    print "\t", $sum_for_gene{$gene_symbol_and_id} // "0";

    for my $barcode (@barcode_list) {
        print "\t", $counts_for{$gene_symbol_and_id}{$barcode} // "0";
    }

    # Newline
    print "\n";
}

my $grand_total = sum0(values %sum_for_gene);

say (join("\t","ALL_GENES",$grand_total, @sum_for_barcode{ @barcode_list })); 

sub gene_symbol_and_id {
    my $symbol = shift;
    my $id = shift;
    return "$symbol ($id)";
}
