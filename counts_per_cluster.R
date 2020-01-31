library(Seurat)
clusters <- readRDS("objects/clusters_pre_mRNA.rds")

cat("sample")
for (sample_id in sort(unique(clusters$abbreviated_sample_name))) {
    cat("\t")
    cat(sample_id)
}
cat("\n")

for (cluster_idx in sort(unique(clusters$seurat_clusters))) {
    cat(cluster_idx)
    for (sample_name in sort(unique(clusters$abbreviated_sample_name))) {
        cluster_name_total = sum(clusters$seurat_clusters == cluster_idx & clusters$abbreviated_sample_name == sample_name)
        cat("\t")
        cat(cluster_name_total)
    }
    cat("\n")
}
