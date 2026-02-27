# ------------------------------
# Load Required Libraries
# ------------------------------
library(readxl)
library(dplyr)
library(pheatmap)
library(writexl)

# ------------------------------
# Step 1: Load Expression Data
# ------------------------------
expr_data <- read_excel("/Users/tanyatripathi/Desktop/MS_data/FINAL/CPM_normalized_data.xlsx")

# Set gene names as rownames
expr_data <- expr_data %>% column_to_rownames("hgnc_symbol")

# ------------------------------
# Step 2: Select Day 0 Samples (Healthy and MS)
# ------------------------------
day0_cols <- grep("^d0_(Ind|MS)", colnames(expr_data), value = TRUE)
day0_data <- expr_data[, day0_cols]

# ------------------------------
# Step 3: Min-Max Normalization (Per Gene) with Sample Names Preserved
# ------------------------------
min_max_norm <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  if (max_val == min_val) {
    return(rep(NA, length(x)))  # Cannot normalize genes with no variance
  } else {
    return((x - min_val) / (max_val - min_val))
  }
}

# Apply normalization
normalized_matrix <- t(apply(day0_data, 1, min_max_norm))

# Restore original sample names
colnames(normalized_matrix) <- colnames(day0_data)

# Convert to data frame
normalized_data <- as.data.frame(normalized_matrix)

# Remove genes with NA values
normalized_data <- normalized_data[complete.cases(normalized_data), ]

# ------------------------------
# Step 4: K-means Clustering of Samples
# ------------------------------
set.seed(123)
kmeans_result <- kmeans(t(normalized_data), centers = 2, nstart = 25)

# ------------------------------
# Step 5: Create Annotations for Heatmap
# ------------------------------
# Assign conditions based on column names
conditions <- ifelse(grepl("d0_Ind", colnames(normalized_data)), "Healthy", "MS")

annotation_col <- data.frame(
  Condition = factor(conditions),
  Cluster = factor(kmeans_result$cluster),
  row.names = colnames(normalized_data)
)

# Define annotation colors
ann_colors <- list(
  Condition = c("Healthy" = "#1f77b4", "MS" = "#d62728"),
  Cluster = c("1" = "#e41a1c", "2" = "green")
)

# Ensure factor levels match color keys
annotation_col$Condition <- factor(annotation_col$Condition, levels = names(ann_colors$Condition))
annotation_col$Cluster <- factor(annotation_col$Cluster, levels = names(ann_colors$Cluster))

# ------------------------------
# Step 6: Reorder Columns by Cluster
# ------------------------------
ordered_columns <- order(kmeans_result$cluster)
normalized_data_ordered <- normalized_data[, ordered_columns]
annotation_col_ordered <- annotation_col[ordered_columns, ]

# ------------------------------
# Step 7: Identify Gaps Between Clusters
# ------------------------------
gaps_col <- which(diff(as.numeric(kmeans_result$cluster[ordered_columns])) != 0)

# ------------------------------
# Step 8: Generate Heatmap
# ------------------------------
heatmap_obj <- pheatmap(normalized_data_ordered,
                        cluster_rows = TRUE,
                        cluster_cols = FALSE,
                        annotation_col = annotation_col_ordered,
                        annotation_colors = ann_colors,
                        show_rownames = FALSE,
                        show_colnames = TRUE,
                        scale = "none",
                        clustering_distance_rows = "euclidean",
                        clustering_method = "complete",
                        color = colorRampPalette(c("blue", "white", "red"))(100),
                        gaps_col = gaps_col,
                        main = "Day 0: Clustering of Healthy vs MS",
                        fontsize = 12,
                        annotation_names_col = TRUE,
                        annotation_legend = TRUE)

# ------------------------------
# Step 9: Export Gene Ordering & Expression Summary
# ------------------------------
# Extract row order from clustering
row_order <- heatmap_obj$tree_row$order

# Create a data frame of ordered gene names
ordered_gene_data <- data.frame(
  Gene_Symbol = rownames(normalized_data_ordered)[row_order]
)

# Calculate mean expression for each group
ordered_gene_data$Mean_Healthy <- rowMeans(normalized_data_ordered[row_order, annotation_col_ordered$Condition == "Healthy"])
ordered_gene_data$Mean_MS <- rowMeans(normalized_data_ordered[row_order, annotation_col_ordered$Condition == "MS"])

# Calculate difference (MS - Healthy)
ordered_gene_data$Difference <- ordered_gene_data$Mean_MS - ordered_gene_data$Mean_Healthy

# Export to Excel
write_xlsx(ordered_gene_data, "/Users/tanyatripathi/Desktop/MS_data/MS_Coimbra/day0_ordered_genes.xlsx")

# ------------------------------
# Step 10: Print Summary
# ------------------------------
cat("✅ Heatmap generated and gene order saved to Excel.\n")
cat("Top genes:\n")
print(head(ordered_gene_data))




















# ------------------------------
# Load Required Libraries
# ------------------------------
library(readxl)       # For reading Excel files
library(dplyr)        # For data manipulation
library(factoextra)   # For cluster visualization (fviz_nbclust)
library(ggplot2)      # For customizing plots

# ------------------------------
# Step 1: Load Expression Data
# ------------------------------
expr_data <- read_excel("/Users/tanyatripathi/Desktop/MS_data/FINAL/CPM_normalized_data.xlsx")

# Set gene names as rownames
expr_data <- expr_data %>% column_to_rownames("hgnc_symbol")

# ------------------------------
# Step 2: Select Day 0 Samples (Healthy and MS)
# ------------------------------
day0_cols <- grep("^d0_(Ind|MS)", colnames(expr_data), value = TRUE)
day0_data <- expr_data[, day0_cols]

# ------------------------------
# Step 3: Min-Max Normalization Per Gene
# ------------------------------
# Define a safe min-max function
min_max_norm <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  if (max_val == min_val) {
    return(rep(NA, length(x)))  # return NA if gene has no variance
  } else {
    return((x - min_val) / (max_val - min_val))
  }
}

# Apply normalization row-wise (per gene)
normalized_matrix <- t(apply(day0_data, 1, min_max_norm))
colnames(normalized_matrix) <- colnames(day0_data)

# Convert to data frame and remove genes with NA (due to 0 variance)
normalized_data <- as.data.frame(normalized_matrix)
normalized_data <- normalized_data[complete.cases(normalized_data), ]

# ------------------------------
# Step 4: Determine Optimal Number of Clusters (WSS Elbow Method)
# ------------------------------
# Transpose data so that samples are rows and genes are columns
# This is the input format expected by kmeans
set.seed(123)  # for reproducibility
fviz_nbclust(
  t(normalized_data),      # samples as rows
  kmeans,                  # clustering function
  method = "wss",          # within-cluster sum of squares
  k.max = 10               # try cluster numbers from 1 to 10
) +
  ggtitle("Elbow Method for Optimal Cluster Number (Day 0 Samples)") +
  theme_minimal()



