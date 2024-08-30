#####################################################
## Metabolomic analysis of OMM12 timecourse.       ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Heatmap clustering
######################################################################################################################

int_scaled <- read.csv("results/ms_scaled_sorted.csv")
rownames(int_scaled) <- int_scaled$X
int_scaled <- int_scaled[,-1]

#create hierarchical clustering for metabolites
my_hclust_metabolites <- hclust(dist(int_scaled), method = "complete")
dend <- as.dendrogram(my_hclust_metabolites)
clust_col <- cutree(tree = dend, k = 4)
clust_col2 <- data.frame(metabolite = names((clust_col)), cluster = unname((clust_col)))
clust_col2$cluster <- factor(clust_col2$cluster, levels = 1:4)
rownames(clust_col2) <- clust_col2$metabolite
clust_col2 <- clust_col2[, 2, drop = FALSE]


# Put in same order as scaled matrix
int_scaled <- int_scaled[match(rownames(clust_col2), rownames(int_scaled)),]

# Define colours for clusters
annotation_colors_list <- list(
  day = setNames(colorRampPalette(c("black", "white"))(length(levels(ann_df$day))), levels(ann_df$day)),
  hour = setNames(brewer.pal(n = length(levels(ann_df$hour)), name = "PiYG"), levels(ann_df$hour)),
  cluster = setNames(c("violetred3", "tomato3", "tomato4", "violet"), c("1", "2", "3", "4"))
)


# Create the heatmap with the sorted columns
heatmap <- pheatmap(ms_scaled_sorted, show_rownames = F, show_colnames = F, cluster_cols = F, annotation_col = ann_df_sorted, annotation_row = clust_col2, annotation_colors = annotation_colors_list, cutree_rows = 4)
save_pheatmap_pdf(heatmap, filename = "plots/heatmap_metabolites_sorted_cluster_coloured.pdf", width = 7, height = 10)


# Get clusters from heatmap
######################################################################################################################
int_scaled_clust <- cbind(int_scaled, cluster = cutree(heatmap$tree_row, k = 4))

write.csv(int_scaled_clust, "results/_int_scaled_cluster_heatmap.csv")

# Split the data by cluster
int_scaled_clust <- as.data.frame(int_scaled_clust)
int_scaled_clust$metabolite <- rownames(int_scaled_clust)
int_scaled_clust <- int_scaled_clust %>% 
  relocate(metabolite)

# Split the data by cluster and keep the row names
cluster_groups <- split(int_scaled_clust, int_scaled_clust$cluster, drop = TRUE)


# Create a new Excel workbook
wb <- createWorkbook()

# Loop through each cluster group
used_names <- character()

# Loop through each cluster group
for (i in seq_along(names(cluster_groups))) {
  # Extract the cluster data
  cluster_num <- names(cluster_groups)[i]
  cluster_data <- cluster_groups[[cluster_num]]
  
  # Create a new sheet name
  sheet_name <- paste("Cluster", cluster_num, sep = "_")
  
  # Check if the sheet name is already used, if yes, modify it
  while (sheet_name %in% used_names) {
    cluster_num <- cluster_num + 1
    sheet_name <- paste("Cluster", cluster_num, sep = "_")
  }
  
  # Add the sheet name to the vector of used names
  used_names <- c(used_names, sheet_name)
  
  # Add a new worksheet and write data
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = cluster_data, startRow = 1, startCol = 1)
}

# Save the workbook to a file
saveWorkbook(wb, "results/cluster_metabolites.xlsx", overwrite = TRUE)
