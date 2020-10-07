library(ggplot2)
library(tidyverse)
library(pheatmap)

# Import RNA-seq data (returns a data frame)
cov <- read.csv("covid19.csv") %>%
        select(-X1)



# Explore your data frame 
print(dim(cov))   # number of rows and columns
print(cov[1:10, 1:10])  # first 10 rows and 10 columns



########################################## PCA ########################################## 


# Run PCA (returns a PCA object)
pca <- prcomp(cov[, 4:ncol(cov)], 
              center = TRUE,    # Expression counts are zero-centered
              scale = TRUE)     # Expression counts are brought to the same scale 
                                # to avoid dominance by highly expressed genes 
                   

# Check your PCA result 
# Note that when your total variance is set to 1, 
# count how many PCs (dimensions) are needed to explain and cluster your data 
# effectively 
summary(pca)

# Determine minimum number of PCs (dimensions) for your data 

# Extract proportion of variance explained by PC1-10 
# and store as a data frame pca_scree
pve <- pca$sdev^2 / sum(pca$sdev^2) 
cum_pve <- cumsum(pve)
pca_df <- data.frame(prop_var = pve[1:10], 
                     cum_prop_var = cum_pve[1:10], 
                     PC = factor(1:10))

print(pca_df)

# Create a plot about proportion of variance change 
prop_var_plot <- ggplot(pca_df,
                        aes(x = PC, 
                            y = prop_var,
                            group = 1)) +
        geom_line() + 
        geom_point() +
        labs(title = "Proportion of Variance Explained by PC1-10",
             y = "Proportion of Variance Explained")

print(prop_var_plot)



# Create a plot about cumulative proportion of variance change 
cum_prop_var_plot <- ggplot(pca_df,
                            aes(x = PC, 
                                y = cum_prop_var,
                                group = 1)) +
        geom_line() + 
        geom_point() +
        labs(title = "Cumulative Proportion of Variance Explained by PC1-10",
             y = "Cumulative Proportion of Variance Explained")

print(cum_prop_var_plot)


# Extract PCA coordinates (PC1-PC3) and clean data 
# (returns a data frame)
pca_coord <- data.frame(Sample = cov$Sample,
                        X = pca$x[, 1],
                        Y = pca$x[, 2],
                        Z = pca$x[, 3]) %>% 
        inner_join(cov[, 1:3], by = "Sample")

print(pca_coord)


# Plot the PCA result in 2D space 
pca_2D_plot <- ggplot(pca_coord,
                      aes(x = X, 
                          y = Y, 
                          color = ICU,
                          shape = gender)) + 
        geom_point(size = 2, alpha = 0.5) +
        labs(title = "PCA",
             x = "PC1 (28%)",
             y = "PC2 (19%)")

print(pca_2D_plot)

########################################## Hierarchical Clustering ##########################################


# Create a coordinate matrix 
coord_matrix <- cov[, 1:3] %>%
        unite(sample_code, gender, ICU, Sample) %>%
        cbind(pca_coord[, 2:4]) %>%
        column_to_rownames(var = "sample_code") %>%
        as.matrix()




# Calculate distance: (X, Y, Z) coordinates
# (returns a distance matrix)
distance <- dist(coord_matrix, 
                 method = "euclidean")

# Perform hierarchical clustering 
Hierarchical_clustering <- hclust(distance, 
                                  method = "average")

# Create a dendrogram
plot(Hierarchical_clustering)

# Heatmap + dendrogram
meta <- cov[, 1:2] 
rownames(meta) <- rownames(coord_matrix)

pheatmap(coord_matrix,
         clustering_distance_rows = "euclidean",
         annotation_row = meta,
         fontsize_row = 5,
         main = "Heatmap")