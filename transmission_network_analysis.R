getwd()
setwd("D:/Dropbox/bovineWGSA/BTBGhAnalysis/transmission_network")

###1. Required R packages
install.packages(c("tidyverse","igraph","pheatmap","ggraph","tidygraph","RColorBrewer","ggrepel"))
install.packages("ape")

install.packages("BiocManager")
BiocManager::install("ggtree")


###2. Load SNP distance matrix
#Generate SNP distances beforehand using pairsnp or snp-dists
#pairsnp -sc filtered_alignment_122cn.fas > filtered_alignment_122cn_snp_matrix_sparse.csv
#snp-dists aligned_pseudogenomes_122cn_masked.fas > aligned_pseudogenomes_122cn_masked_distances.tab


library(tidyverse)

dist_matrix <- read.table("filtered_alignment_122cn_snp_matrix.tsv",
                          header=TRUE,
                          sep = "\t",
                          row.names=1,
                          check.names=FALSE)

dist_matrix <- as.matrix(dist_matrix)


###3. Define SNP transmission threshold
#NB. repeat this step for cluster sizes of 5 and 10
threshold <- 12


###4. Build transmission network
#Convert SNP distances into a graph.
library(igraph)

adj_matrix <- dist_matrix <= threshold
diag(adj_matrix) <- 0

g <- graph_from_adjacency_matrix(adj_matrix,
                                 mode="undirected")

###5a. Identify clusters
clusters <- components(g)

cluster_df <- data.frame(
  sample = names(clusters$membership),
  old_cluster = clusters$membership
)

### Renumber clusters so largest cluster = 1
cluster_order <- cluster_df %>%
  count(old_cluster, name = "cluster_size") %>%
  arrange(desc(cluster_size), old_cluster) %>%
  mutate(cluster = row_number())

cluster_df <- cluster_df %>%
  left_join(cluster_order, by = "old_cluster") %>%
  select(sample, cluster, cluster_size, old_cluster)

print(cluster_df)

###5b Assign reordered cluster membership to graph vertices and calculate cluster size
V(g)$cluster <- cluster_df$cluster[match(V(g)$name, cluster_df$sample)]
V(g)$cluster_size <- cluster_df$cluster_size[match(V(g)$name, cluster_df$sample)]

###5d. Remove singleton clusters
g_clustered <- induced_subgraph(g, vids = V(g)[cluster_size > 1])

###5e. Export cluster table
write.csv(cluster_df,
          "transmission_clusters_threshold_12SNPs.csv",
          row.names = FALSE)


###6a. Attach metadata. Remember to recode "NA" to "MD" for region.
meta_all <- read.csv("metadata.csv", stringsAsFactors = FALSE)
#Use filtered metadata only for network plotting
meta_net <- meta_all[match(V(g_clustered)$name, meta_all$Genome_ID),]

V(g_clustered)$district <- meta_net$District
V(g_clustered)$ID <- meta_net$Genome_ID
V(g_clustered)$gender <- meta_net$Gender
V(g_clustered)$region <- meta_net$Region
V(g_clustered)$breed <- meta_net$Cattle_Breed
V(g_clustered)$sb <- meta_net$SB_number
V(g_clustered)$host <- meta_net$Host


###6b. Convert to tidygraph

library(tidygraph)
library(RColorBrewer)

tg <- as_tbl_graph(g_clustered)


###6c. Define color palettes
## District colors (~20 categories)
districts <- sort(unique(meta_net$District))

district_colors <- setNames(
  colorRampPalette(brewer.pal(12,"Set3"))(length(districts)),
  districts
)

##Region shapes 
region_shapes <- c(North_East = 21, Northern = 22, Savannah = 23, Upper_West = 24, MD = 25)

## Create edge attribute (same vs different breed)
tg <- tg %>%
  activate(edges) %>%
  mutate(
    same_breed = case_when(
      is.na(.N()$breed[from]) | is.na(.N()$breed[to]) ~ "NA",
      .N()$breed[from] == .N()$breed[to] ~ "Same breed",
      TRUE ~ "Different breed"
    )
  )

###6d. Build the transmission network (ggraph)

library(ggraph)

set.seed(123)

ggraph(tg, layout = "fr") +
  
  # edges with breed logic
  geom_edge_link(aes(color = same_breed),
                 width = 1.2,
                 alpha = 0.7) +
  
   # nodes
  geom_node_point(aes(fill = district,
                      shape = region,
                      color = gender),
                  size = 9,
                  stroke = 1.2) +
  
  # labels (SB number with ID). You can remove the SB numbers manually after plotting and rearranging to produce a cleaner image.
  geom_node_text(aes(label = paste0(sb, "\n", ID)),
                 repel = TRUE,
                 size = 3) +
  
  # edge colors
  scale_edge_color_manual(
    values = c(
      "Same breed" = "#1b9e77",
      "Different breed" = "#d95f02",
      "NA" = "grey70"
    )
  ) +
  
  # scales (node aesthetics)
  scale_fill_manual(values = district_colors) +
  scale_shape_manual(values = region_shapes) +
  scale_color_manual(values = c(Male="blue", Female="red")) +
  
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  
  # theme
  theme_graph() +
  labs(title = "Putative Transmission Network (M. bovis)",
       edge_color = "Cattle Breed",
       fill = "District",
       shape = "Region",
       color = "Gender")



###6e. Export publication-quality figure
ggsave("transmission_network_threshold_12SNPs.svg",
       width = 14,
       height = 14,
       dpi = 400)


###7. Heatmap of SNP distances
library(pheatmap)
###7a. Heatmap of SNP distances for all cases
pheatmap(dist_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE)

###7b. Plot Heatmap for only clusters
#Identify clustered isolates. First, remove singletons (isolates not connected to others).
cluster_sizes <- table(cluster_df$cluster)

clustered <- cluster_df$sample[cluster_df$cluster %in% names(cluster_sizes[cluster_sizes > 1])]

#Now clustered contains only samples that belong to clusters with ≥2 isolates.

#Subset the SNP distance matrix. Filter the matrix to keep only clustered isolates.
dist_clustered <- dist_matrix[clustered, clustered]

#Plot heatmap for clusters
pheatmap(dist_clustered,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "SNP Distances within Transmission Clusters")

#color clusters on the heatmap
annotation <- data.frame(
  Cluster = factor(cluster_df$cluster[match(clustered, cluster_df$sample)])
)

rownames(annotation) <- clustered

# Create heatmap object
p <- pheatmap(
  dist_clustered,
  annotation_col = annotation,
  annotation_row = annotation
)

# Open SVG device
svg("transmission_cluster_SNP_distances_heatmap_threshold_12SNPs.svg",
    width = 12,
    height = 12)

# Draw plot
grid::grid.newpage()
grid::grid.draw(p$gtable)

# Close device
dev.off()


###8. Calculate Recent Transmission Index (RTI). Formulae: RTI = (clustered cases − number of clusters) / total cases
total_cases <- nrow(cluster_df)

clustered_cases <- cluster_df %>%
  filter(cluster_size > 1) %>%
  nrow()

num_clusters <- cluster_df %>%
  filter(cluster_size > 1) %>%
  distinct(cluster) %>%
  nrow()

RTI <- (clustered_cases - num_clusters) / total_cases

RTI



###9. Calculate Recent Transmission Index (RTI) per SB number using all isolates
#NB. This session requires tidyverse package
# Merge transmission cluster results with metadata
# cluster_df must contain all isolates, including singletons
# meta must contain all 122 isolates

# A. Denominator: all isolates per SB number from full metadata
sb_denominator <- meta_all %>%
  filter(!is.na(SB_number), SB_number != "") %>%
  count(SB_number, name = "total_cases")

# B. Add cluster information to all isolates
cluster_meta_all <- meta_all %>%
  filter(!is.na(SB_number), SB_number != "") %>%
  left_join(cluster_df, by = c("Genome_ID" = "sample")) %>%
  mutate(
    cluster_size = ifelse(is.na(cluster_size), 1, cluster_size),
    clustered = cluster_size > 1
  )

# C. Numerator components per SB number
sb_cluster_summary <- cluster_meta_all %>%
  group_by(SB_number) %>%
  summarise(
    clustered_cases = sum(clustered, na.rm = TRUE),
    num_clusters = n_distinct(cluster[clustered]),
    .groups = "drop"
  )

# D. Calculate SB-specific RTI
sb_rti_all <- sb_denominator %>%
  left_join(sb_cluster_summary, by = "SB_number") %>%
  mutate(
    clustered_cases = replace_na(clustered_cases, 0),
    num_clusters = replace_na(num_clusters, 0),
    RTI = (clustered_cases - num_clusters) / total_cases,
    proportion_clustered = clustered_cases / total_cases
  ) %>%
  arrange(desc(RTI), desc(total_cases))

sb_rti_all

# E. export the results
write.csv(sb_rti_all,
          "RTI_by_SB_number_all_122_cases_threshold_12SNPs.csv",
          row.names = FALSE)


#Filter for major SB numbers only
#keep SB numbers represented by at least 20 isolates:
sb_rti_major <- sb_rti_all %>%
  filter(total_cases >= 20)

write.csv(sb_rti_major,
          "RTI_by_major_SB_number_all_122_cases_threshold_12SNPs.csv",
          row.names = FALSE)

#Plot RTI by SB number
library(ggplot2)

ggplot(sb_rti_major,
       aes(x = reorder(SB_number, RTI),
           y = RTI)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "SB number",
    y = "Recent Transmission Index",
    title = "Recent Transmission Index by major SB number"
  )


#########################################################################################
#########################################################################################
###################### REPEAT ANALYSIS FOR THRESHOLD OF 10 SNPS #########################
#########################################################################################
#########################################################################################

print("Clear the Global Environment and continue from here")
print("Clear the Global Environment and continue from here")
print("Clear the Global Environment and continue from here")
print("Clear the Global Environment and continue from here")

###2. Load SNP distance matrix
library(tidyverse)

dist_matrix <- read.table("filtered_alignment_122cn_snp_matrix.tsv",
                          header=TRUE,
                          sep = "\t",
                          row.names=1,
                          check.names=FALSE)

dist_matrix <- as.matrix(dist_matrix)

###3. Define SNP transmission threshold
threshold <- 10

###4. Build transmission network
#Convert SNP distances into a graph.
library(igraph)

adj_matrix <- dist_matrix <= threshold
diag(adj_matrix) <- 0

g <- graph_from_adjacency_matrix(adj_matrix,
                                 mode="undirected")
###5a. Identify clusters
clusters <- components(g)

cluster_df <- data.frame(
  sample = names(clusters$membership),
  old_cluster = clusters$membership
)

### Renumber clusters so largest cluster = 1
cluster_order <- cluster_df %>%
  count(old_cluster, name = "cluster_size") %>%
  arrange(desc(cluster_size), old_cluster) %>%
  mutate(cluster = row_number())

cluster_df <- cluster_df %>%
  left_join(cluster_order, by = "old_cluster") %>%
  select(sample, cluster, cluster_size, old_cluster)

print(cluster_df)

###5b Assign reordered cluster membership to graph vertices and calculate cluster size
V(g)$cluster <- cluster_df$cluster[match(V(g)$name, cluster_df$sample)]
V(g)$cluster_size <- cluster_df$cluster_size[match(V(g)$name, cluster_df$sample)]

###5d. Remove singleton clusters
g_clustered <- induced_subgraph(g, vids = V(g)[cluster_size > 1])

###5e. Export cluster table
write.csv(cluster_df,
          "transmission_clusters_threshold_10SNPs.csv",
          row.names = FALSE)


###6a. Attach metadata. Remember to recode "NA" to "MD" for region.
meta_all <- read.csv("metadata.csv", stringsAsFactors = FALSE)
#Use filtered metadata only for network plotting
meta_net <- meta_all[match(V(g_clustered)$name, meta_all$Genome_ID),]

V(g_clustered)$district <- meta_net$District
V(g_clustered)$ID <- meta_net$Genome_ID
V(g_clustered)$gender <- meta_net$Gender
V(g_clustered)$region <- meta_net$Region
V(g_clustered)$breed <- meta_net$Cattle_Breed
V(g_clustered)$sb <- meta_net$SB_number
V(g_clustered)$host <- meta_net$Host


###6b. Convert to tidygraph

library(tidygraph)
library(RColorBrewer)

tg <- as_tbl_graph(g_clustered)


###6c. Define color palettes
## District colors (~20 categories)
districts <- sort(unique(meta_net$District))

district_colors <- setNames(
  colorRampPalette(brewer.pal(12,"Set3"))(length(districts)),
  districts
)

##Region shapes 
region_shapes <- c(North_East = 21, Northern = 22, Savannah = 23, Upper_West = 24, MD = 25)

## Create edge attribute (same vs different breed)
tg <- tg %>%
  activate(edges) %>%
  mutate(
    same_breed = case_when(
      is.na(.N()$breed[from]) | is.na(.N()$breed[to]) ~ "NA",
      .N()$breed[from] == .N()$breed[to] ~ "Same breed",
      TRUE ~ "Different breed"
    )
  )

###6d. Build the transmission network (ggraph)

library(ggraph)

set.seed(123)

ggraph(tg, layout = "fr") +
  
  # edges with breed logic
  geom_edge_link(aes(color = same_breed),
                 width = 1.2,
                 alpha = 0.7) +
  
  # nodes
  geom_node_point(aes(fill = district,
                      shape = region,
                      color = gender),
                  size = 9,
                  stroke = 1.2) +
  
  # labels (SB number with ID). You can remove the SB numbers manually after plotting and rearranging to produce a cleaner image.
  geom_node_text(aes(label = paste0(sb, "\n", ID)),
                 repel = TRUE,
                 size = 3) +
  
  # edge colors
  scale_edge_color_manual(
    values = c(
      "Same breed" = "#1b9e77",
      "Different breed" = "#d95f02",
      "NA" = "grey70"
    )
  ) +
  
  # scales (node aesthetics)
  scale_fill_manual(values = district_colors) +
  scale_shape_manual(values = region_shapes) +
  scale_color_manual(values = c(Male="blue", Female="red")) +
  
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  
  # theme
  theme_graph() +
  labs(title = "Putative Transmission Network (M. bovis)",
       edge_color = "Cattle Breed",
       fill = "District",
       shape = "Region",
       color = "Gender")



###6e. Export publication-quality figure
ggsave("transmission_network_threshold_10SNPs.svg",
       width = 14,
       height = 14,
       dpi = 400)


###7. Heatmap of SNP distances
library(pheatmap)
###7a. Heatmap of SNP distances for all cases
pheatmap(dist_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE)

###7b. Plot Heatmap for only clusters
#Identify clustered isolates. First, remove singletons (isolates not connected to others).
cluster_sizes <- table(cluster_df$cluster)

clustered <- cluster_df$sample[cluster_df$cluster %in% names(cluster_sizes[cluster_sizes > 1])]

#Now clustered contains only samples that belong to clusters with ≥2 isolates.

#Subset the SNP distance matrix. Filter the matrix to keep only clustered isolates.
dist_clustered <- dist_matrix[clustered, clustered]

#Plot heatmap for clusters
pheatmap(dist_clustered,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "SNP Distances within Transmission Clusters")

#color clusters on the heatmap
annotation <- data.frame(
  Cluster = factor(cluster_df$cluster[match(clustered, cluster_df$sample)])
)

rownames(annotation) <- clustered

# Create heatmap object
p <- pheatmap(
  dist_clustered,
  annotation_col = annotation,
  annotation_row = annotation
)

# Open SVG device
svg("transmission_cluster_SNP_distances_heatmap_threshold_10SNPs.svg",
    width = 12,
    height = 12)

# Draw plot
grid::grid.newpage()
grid::grid.draw(p$gtable)

# Close device
dev.off()


###8. Calculate Recent Transmission Index (RTI). Formulae: RTI = (clustered cases − number of clusters) / total cases
total_cases <- nrow(cluster_df)

clustered_cases <- cluster_df %>%
  filter(cluster_size > 1) %>%
  nrow()

num_clusters <- cluster_df %>%
  filter(cluster_size > 1) %>%
  distinct(cluster) %>%
  nrow()

RTI <- (clustered_cases - num_clusters) / total_cases

RTI



###9. Calculate Recent Transmission Index (RTI) per SB number using all isolates
#NB. This session requires tidyverse package
# Merge transmission cluster results with metadata
# cluster_df must contain all isolates, including singletons
# meta must contain all 122 isolates

# A. Denominator: all isolates per SB number from full metadata
sb_denominator <- meta_all %>%
  filter(!is.na(SB_number), SB_number != "") %>%
  count(SB_number, name = "total_cases")

# B. Add cluster information to all isolates
cluster_meta_all <- meta_all %>%
  filter(!is.na(SB_number), SB_number != "") %>%
  left_join(cluster_df, by = c("Genome_ID" = "sample")) %>%
  mutate(
    cluster_size = ifelse(is.na(cluster_size), 1, cluster_size),
    clustered = cluster_size > 1
  )

# C. Numerator components per SB number
sb_cluster_summary <- cluster_meta_all %>%
  group_by(SB_number) %>%
  summarise(
    clustered_cases = sum(clustered, na.rm = TRUE),
    num_clusters = n_distinct(cluster[clustered]),
    .groups = "drop"
  )

# D. Calculate SB-specific RTI
sb_rti_all <- sb_denominator %>%
  left_join(sb_cluster_summary, by = "SB_number") %>%
  mutate(
    clustered_cases = replace_na(clustered_cases, 0),
    num_clusters = replace_na(num_clusters, 0),
    RTI = (clustered_cases - num_clusters) / total_cases,
    proportion_clustered = clustered_cases / total_cases
  ) %>%
  arrange(desc(RTI), desc(total_cases))

sb_rti_all

# E. export the results
write.csv(sb_rti_all,
          "RTI_by_SB_number_all_122_cases_threshold_10SNPs.csv",
          row.names = FALSE)


#Filter for major SB numbers only
#keep SB numbers represented by at least 20 isolates:
sb_rti_major <- sb_rti_all %>%
  filter(total_cases >= 20)

write.csv(sb_rti_major,
          "RTI_by_major_SB_number_all_122_cases_threshold_10SNPs.csv",
          row.names = FALSE)


#Plot RTI by SB number
library(ggplot2)

ggplot(sb_rti_major,
       aes(x = reorder(SB_number, RTI),
           y = RTI)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "SB number",
    y = "Recent Transmission Index",
    title = "Recent Transmission Index by major SB number"
  )



#########################################################################################
#########################################################################################
###################### REPEAT ANALYSIS FOR THRESHOLD OF 5 SNPS #########################
#########################################################################################
#########################################################################################

print("Clear the Global Environment and continue from here")
print("Clear the Global Environment and continue from here")
print("Clear the Global Environment and continue from here")
print("Clear the Global Environment and continue from here")

###2. Load SNP distance matrix
library(tidyverse)

dist_matrix <- read.table("filtered_alignment_122cn_snp_matrix.tsv",
                          header=TRUE,
                          sep = "\t",
                          row.names=1,
                          check.names=FALSE)

dist_matrix <- as.matrix(dist_matrix)

###3. Define SNP transmission threshold
threshold <- 5

###4. Build transmission network
#Convert SNP distances into a graph.
library(igraph)

adj_matrix <- dist_matrix <= threshold
diag(adj_matrix) <- 0

g <- graph_from_adjacency_matrix(adj_matrix,
                                 mode="undirected")
###5a. Identify clusters
clusters <- components(g)

cluster_df <- data.frame(
  sample = names(clusters$membership),
  old_cluster = clusters$membership
)

### Renumber clusters so largest cluster = 1
cluster_order <- cluster_df %>%
  count(old_cluster, name = "cluster_size") %>%
  arrange(desc(cluster_size), old_cluster) %>%
  mutate(cluster = row_number())

cluster_df <- cluster_df %>%
  left_join(cluster_order, by = "old_cluster") %>%
  select(sample, cluster, cluster_size, old_cluster)

print(cluster_df)

###5b Assign reordered cluster membership to graph vertices and calculate cluster size
V(g)$cluster <- cluster_df$cluster[match(V(g)$name, cluster_df$sample)]
V(g)$cluster_size <- cluster_df$cluster_size[match(V(g)$name, cluster_df$sample)]

###5d. Remove singleton clusters
g_clustered <- induced_subgraph(g, vids = V(g)[cluster_size > 1])

###5e. Export cluster table
write.csv(cluster_df,
          "transmission_clusters_threshold_5SNPs.csv",
          row.names = FALSE)


###6a. Attach metadata. Remember to recode "NA" to "MD" for region.
meta_all <- read.csv("metadata.csv", stringsAsFactors = FALSE)
#Use filtered metadata only for network plotting
meta_net <- meta_all[match(V(g_clustered)$name, meta_all$Genome_ID),]

V(g_clustered)$district <- meta_net$District
V(g_clustered)$ID <- meta_net$Genome_ID
V(g_clustered)$gender <- meta_net$Gender
V(g_clustered)$region <- meta_net$Region
V(g_clustered)$breed <- meta_net$Cattle_Breed
V(g_clustered)$sb <- meta_net$SB_number
V(g_clustered)$host <- meta_net$Host


###6b. Convert to tidygraph

library(tidygraph)
library(RColorBrewer)

tg <- as_tbl_graph(g_clustered)


###6c. Define color palettes
## District colors (~20 categories)
districts <- sort(unique(meta_net$District))

district_colors <- setNames(
  colorRampPalette(brewer.pal(12,"Set3"))(length(districts)),
  districts
)

##Region shapes 
region_shapes <- c(North_East = 21, Northern = 22, Savannah = 23, Upper_West = 24, MD = 25)

## Create edge attribute (same vs different breed)
tg <- tg %>%
  activate(edges) %>%
  mutate(
    same_breed = case_when(
      is.na(.N()$breed[from]) | is.na(.N()$breed[to]) ~ "NA",
      .N()$breed[from] == .N()$breed[to] ~ "Same breed",
      TRUE ~ "Different breed"
    )
  )

###6d. Build the transmission network (ggraph)

library(ggraph)

set.seed(123)

ggraph(tg, layout = "fr") +
  
  # edges with breed logic
  geom_edge_link(aes(color = same_breed),
                 width = 1.2,
                 alpha = 0.7) +
  
  # nodes
  geom_node_point(aes(fill = district,
                      shape = region,
                      color = gender),
                  size = 9,
                  stroke = 1.2) +
  
  # labels (SB number with ID). You can remove the SB numbers manually after plotting and rearranging to produce a cleaner image.
  geom_node_text(aes(label = paste0(sb, "\n", ID)),
                 repel = TRUE,
                 size = 3) +
  
  # edge colors
  scale_edge_color_manual(
    values = c(
      "Same breed" = "#1b9e77",
      "Different breed" = "#d95f02",
      "NA" = "grey70"
    )
  ) +
  
  # scales (node aesthetics)
  scale_fill_manual(values = district_colors) +
  scale_shape_manual(values = region_shapes) +
  scale_color_manual(values = c(Male="blue", Female="red")) +
  
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  
  # theme
  theme_graph() +
  labs(title = "Putative Transmission Network (M. bovis)",
       edge_color = "Cattle Breed",
       fill = "District",
       shape = "Region",
       color = "Gender")



###6e. Export publication-quality figure
ggsave("transmission_network_threshold_5SNPs.svg",
       width = 14,
       height = 14,
       dpi = 400)


###7. Heatmap of SNP distances
library(pheatmap)
###7a. Heatmap of SNP distances for all cases
pheatmap(dist_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE)

###7b. Plot Heatmap for only clusters
#Identify clustered isolates. First, remove singletons (isolates not connected to others).
cluster_sizes <- table(cluster_df$cluster)

clustered <- cluster_df$sample[cluster_df$cluster %in% names(cluster_sizes[cluster_sizes > 1])]

#Now clustered contains only samples that belong to clusters with ≥2 isolates.

#Subset the SNP distance matrix. Filter the matrix to keep only clustered isolates.
dist_clustered <- dist_matrix[clustered, clustered]

#Plot heatmap for clusters
pheatmap(dist_clustered,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "SNP Distances within Transmission Clusters")

#color clusters on the heatmap
annotation <- data.frame(
  Cluster = factor(cluster_df$cluster[match(clustered, cluster_df$sample)])
)

rownames(annotation) <- clustered

# Create heatmap object
p <- pheatmap(
  dist_clustered,
  annotation_col = annotation,
  annotation_row = annotation
)

# Open SVG device
svg("transmission_cluster_SNP_distances_heatmap_threshold_5SNPs.svg",
    width = 12,
    height = 12)

# Draw plot
grid::grid.newpage()
grid::grid.draw(p$gtable)

# Close device
dev.off()


###8. Calculate Recent Transmission Index (RTI). Formulae: RTI = (clustered cases − number of clusters) / total cases
total_cases <- nrow(cluster_df)

clustered_cases <- cluster_df %>%
  filter(cluster_size > 1) %>%
  nrow()

num_clusters <- cluster_df %>%
  filter(cluster_size > 1) %>%
  distinct(cluster) %>%
  nrow()

RTI <- (clustered_cases - num_clusters) / total_cases

RTI



###9. Calculate Recent Transmission Index (RTI) per SB number using all isolates
#NB. This session requires tidyverse package
# Merge transmission cluster results with metadata
# cluster_df must contain all isolates, including singletons
# meta must contain all 122 isolates

# A. Denominator: all isolates per SB number from full metadata
sb_denominator <- meta_all %>%
  filter(!is.na(SB_number), SB_number != "") %>%
  count(SB_number, name = "total_cases")

# B. Add cluster information to all isolates
cluster_meta_all <- meta_all %>%
  filter(!is.na(SB_number), SB_number != "") %>%
  left_join(cluster_df, by = c("Genome_ID" = "sample")) %>%
  mutate(
    cluster_size = ifelse(is.na(cluster_size), 1, cluster_size),
    clustered = cluster_size > 1
  )

# C. Numerator components per SB number
sb_cluster_summary <- cluster_meta_all %>%
  group_by(SB_number) %>%
  summarise(
    clustered_cases = sum(clustered, na.rm = TRUE),
    num_clusters = n_distinct(cluster[clustered]),
    .groups = "drop"
  )

# D. Calculate SB-specific RTI
sb_rti_all <- sb_denominator %>%
  left_join(sb_cluster_summary, by = "SB_number") %>%
  mutate(
    clustered_cases = replace_na(clustered_cases, 0),
    num_clusters = replace_na(num_clusters, 0),
    RTI = (clustered_cases - num_clusters) / total_cases,
    proportion_clustered = clustered_cases / total_cases
  ) %>%
  arrange(desc(RTI), desc(total_cases))

sb_rti_all

# E. export the results
write.csv(sb_rti_all,
          "RTI_by_SB_number_all_122_cases_threshold_5SNPs.csv",
          row.names = FALSE)


#Filter for major SB numbers only
#keep SB numbers represented by at least 20 isolates:
sb_rti_major <- sb_rti_all %>%
  filter(total_cases >= 20)

write.csv(sb_rti_major,
          "RTI_by_major_SB_number_all_122_cases_threshold_5SNPs.csv",
          row.names = FALSE)


#Plot RTI by SB number
library(ggplot2)

ggplot(sb_rti_major,
       aes(x = reorder(SB_number, RTI),
           y = RTI)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "SB number",
    y = "Recent Transmission Index",
    title = "Recent Transmission Index by major SB number"
  )





