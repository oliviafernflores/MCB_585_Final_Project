# Load libraries
library(tidyverse)

# Set seed for reproducibility
set.seed(123)

# Read in raw, upregulated, and downregulated data
df <- read.csv("/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/xenograft_data_normalized.csv")
up <- read.csv('/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/xenograft_upregulated.csv')
down <- read.csv('/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/xenograft_downregulated.csv')

# Get a list of genes that are all upregulated or downregulated
valid_genes <- c(up$Probe, down$Probe)
valid_genes <- ifelse(grepl("^\\d", valid_genes), paste0("X", valid_genes), valid_genes)


# Create a new dataframe, df_numeric, with only numeric data and only the differentially expressed (either up or down) genes
df_numeric<-df %>%mutate(across(everything(), ~as.numeric(.)))

df_numeric <- df_numeric %>% select(any_of(valid_genes))

# df_numeric <- df_numeric[1:(nrow(df_numeric)-1), ]

# Perform PCA using prcomp, https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
pca_result <- prcomp(df_numeric, center = TRUE, scale. = TRUE)
top_10_pcs <- pca_result$x[, 1:10]

# Get information about short or long time to leukemia
time_status <- df[1:(nrow(df)-1), "Time_Status"]

#Plot PC1 vs PC2, color by Time_Status
time_status <- as.factor(time_status)
plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = time_status, pch = 19, 
     xlab = "PC1", ylab = "PC2", 
     main = "PC1 vs PC2")
legend("bottom", legend = levels(time_status), 
       col = 1:length(levels(time_status)), pch = 19)


# Perform K-means with 2 centers, https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/kmeans
kmeans_result<- kmeans(df_numeric, centers = 2)

#Plot PC1 vs PC2, color by Time_Status, shape is based on cluster assignment
plot(pca_result$x[, 1], pca_result$x[, 2], 
     col = as.integer(factor(time_status)),
     pch =as.integer(kmeans_result$cluster),
     xlab = "PC1", ylab = "PC2", 
     main = "PC1 vs. PC2, K-means with 2 Clusters")
legend("top", legend = c("C1, Long_TTL", "C2, Long_TTL", "C1, Short_TTL", "C2, Short_TTL"), 
       col = c(1, 1, 2, 2), pch = c(1, 2, 1, 2), title = "K-means Clusters")


# Perform hierarchical clustering, followed https://www.datacamp.com/tutorial/hierarchical-clustering-R as an example implementation in R
dist_matrix <- dist(df_numeric, method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "complete")
df_numeric_matrix <- as.matrix(df_numeric)

# Plot a dendrogram of the hierarchical clustering results
plot(hclust_result, main = "Hierarchical Clustering Dendrogram",
     xlab = "Samples", ylab = "Height", col = "darkblue", cex = 0.6)

# Create the scale of colors for the heatmap
my_palette <- colorRampPalette(c("red", "black", "green"))(256)
# Get the labels for each sample (short or long time to leukemia)
custom_labels <- df$Time_Status[1:12]
# Plot the heatmap with the sample labels and clustering results as a dendrogram
heatmap(df_numeric_matrix,
        main = "Xenograft Hierarchical Clustering",
        scale = "row",
        col = my_palette,
        Rowv = as.dendrogram(hclust_result),
        Colv = NA,
        margins = c(6, 6), 
        labRow = custom_labels) 