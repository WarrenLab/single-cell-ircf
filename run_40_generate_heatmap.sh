module load biocompute/biocompute-modules
module load seurat/seurat-3.1.0_R-3.6.1
#MARKERS_FILE=${1}
MARKERS_FILE=markers/FindAllMarkers.tsv

# Find the last cluster ID
LAST_CLUSTER_ID=`tail -n1 ${MARKERS_FILE}| awk '{print $7}'`

# # Create lists of markers for each cluster
# for i in `seq 0 ${LAST_CLUSTER_ID}`; do awk -v var=$i 'BEGIN{OFS="\t"} $7 == var {print $8,$2,$3,$4,$5,$6}' ${MARKERS_FILE} > markers_for_${i}.tsv; done

# Create incrementally larger sets of "Top markers"
TOP_COUNTS=("3" "5" "10")

for TOP_COUNT in "${TOP_COUNTS[@]}"; do

    TOP_MARKER_LIST_FILE="markers/list_of_top_${TOP_COUNT}_markers.txt"

    # Create (or truncate) marker list file
    truncate -s 0 ${TOP_MARKER_LIST_FILE}
    for i in `seq 0 ${LAST_CLUSTER_ID}`; do echo -e "\n#cluster${i}" >> ${TOP_MARKER_LIST_FILE}; awk -v clusterID=$i '$7 == clusterID {print $8}' ${MARKERS_FILE} | head -n ${TOP_COUNT} >> ${TOP_MARKER_LIST_FILE}; done

    # Create code file
    TOP_MARKER_R_LIST_FILE="markers/Generate_heatmap_for_top_${TOP_COUNT}_markers.R"

    cp 40_top_genes_R_template $TOP_MARKER_R_LIST_FILE

    echo "top_${TOP_COUNT} <- c(" >> $TOP_MARKER_R_LIST_FILE

    awk '/^#/ {print "\n\t" $1} /^\w/{print "\t\"" $1 "\","}' $TOP_MARKER_LIST_FILE >> $TOP_MARKER_R_LIST_FILE 

    # Replace last comma in the file with close parenthesis
    sed -i '$s/,$/)/' ${TOP_MARKER_R_LIST_FILE}  # The first dollar sign makes this match only the last line in the file

    cat <<END_OF_SCRIPT >> ${TOP_MARKER_R_LIST_FILE}

# Draw a heatmap of all cells for these marker genes
# Gene names are a bit too small to see in the tutorial, but you can blow up the heatmap in R
pdf('plots/Heatmap_top_${TOP_COUNT}_genes_per_cluster.pdf')
    print(DoHeatmap(clusters,features=top_${TOP_COUNT}))
dev.off()
END_OF_SCRIPT

    sbatch ${TOP_MARKER_R_LIST_FILE}
done

