library(Seurat)
clusters <- readRDS("seurat.rds")

grand_total <- 0
for (cluster_idx in sort(unique(clusters$seurat_clusters))) {
    cat(cluster_idx)
    cluster_name_total = sum(clusters$seurat_clusters == cluster_idx)
    cat("\t")
    cat(cluster_name_total)
    grand_total <- grand_total + cluster_name_total
    cat("\n")
}
cat("TOTAL\t")
cat(grand_total)
cat("\n")


big_count_table <- as.data.frame(clusters$RNA@counts)
raw_all_feature_sum <- rowSums(big_count_table)

for (cluster_idx in sort(unique(clusters$seurat_clusters))) {
    selected_barcodes <- names(clusters$seurat_clusters[clusters$seurat_clusters == cluster_idx])
    cluster_specific <- big_count_table[selected_barcodes]

    cluster_specific$feature_raw_sum     <- rowSums(cluster_specific[selected_barcodes])
    cluster_specific$feature_raw_percent <- round(100 * (cluster_specific$feature_raw_sum / raw_all_feature_sum), 2)

    cells_in_cluster                     <- length(selected_barcodes)
    cluster_specific$feature_mean_count  <- round(cluster_specific$feature_raw_sum / cells_in_cluster, 2)

    write.table(cluster_specific,
                 file=paste0("cluster_",
                    cluster_idx,
                    ".counts.txt"
                 ),
                 sep="\t",
                 quote=FALSE
    )

    write.table(cluster_specific[c("feature_raw_sum","feature_raw_percent","feature_mean_count")],
                row.names=rownames(cluster_specific),
                col.names=TRUE,
                file=paste0("feature_stats_for_",
                   cluster_idx,
                   ".txt"
                ),
                sep="\t",
                quote=FALSE
    )

}

write.table(big_count_table,
             file="all_feature_counts.txt",
             sep="\t",
             quote=FALSE
)
