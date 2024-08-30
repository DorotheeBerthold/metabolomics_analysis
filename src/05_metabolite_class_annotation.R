#####################################################
## Metabolomic analysis of OMM12 timecourse.       ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

library(gridExtra)
library(grid)
# Clean up metabolite class annotations
######################################################################################################################

metabolite_classes <- read.csv("tables/metabolite_classes.csv")
int_clust <- read.csv("results/_int_scaled_cluster_heatmap.csv")
colnames(int_clust)[1] <- "Metabolite"
int_clust <- int_clust[,c(1,32)]
int_clust$Metabolite <- gsub("D-|L-| 1-PHOSPHATE| 6-PHOSPHATE", "", int_clust$Metabolite)
int_clust$Metabolite <- trimws(int_clust$Metabolite)

sum_classes_before_ann <- int_clust %>% 
  group_by(cluster) %>% 
  summarise(count = n())




metabolite_classes$BIOCHEMICAL <- toupper(metabolite_classes$BIOCHEMICAL)
metabolite_classes <- metabolite_classes[,c(1:3)]
colnames(metabolite_classes)[1:3] <- c("Metabolite", "Superclass", "Class")

metabolite_cluster_annotated <- inner_join(int_clust, metabolite_classes, by = "Metabolite")

not_annotated <- anti_join(int_clust, metabolite_classes, by = "Metabolite")


#set distinct colours for the superclasses
unique_superclasses <- unique(metabolite_cluster_annotated$Superclass)

# Define a palette of colors with enough colors for all unique superclasses
palette <- viridis(length(unique_superclasses))

# Create a named vector to map each superclass to a color
color_mapping <- setNames(palette, unique_superclasses)


plots <- list()

for (i in 1:4) {
  int <- metabolite_cluster_annotated %>%
    filter(cluster == i) %>%
    group_by(Superclass) %>%
    summarise(count = n())
  
  p <- ggplot(int, aes(x = "", y = count, fill = Superclass)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = color_mapping) +
    theme_void() +
    labs(title = paste("Cluster", i))
  
  plots[[i]] <- p
}

grid.arrange(grobs = plots, ncol = 2)

ggsave(paste0("plots/superclass_distribution_heatmap.pdf"), 
       arrangeGrob(grobs = plots, ncol = 3, 
                   top = textGrob("Metabolite Superclass Distribution by Cluster", 
                                  gp = gpar(fontsize = 16, fontface = "bold"))), width = 12, height = 8)

# Check how many metabolites in each cluster
sum_classes <- metabolite_cluster_annotated %>% 
  group_by(cluster) %>% 
  summarise(count = n())
