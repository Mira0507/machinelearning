library(ggplot2)
library(tidyverse)
library(pheatmap)

# Import RNA-seq data (returns a data frame)
cov <- read.csv("covid19.csv") %>%
        select(-X1)



# Explore your data frame 
print(dim(cov))   # number of rows and columns
print(cov[1:10, 1:10])  # first 10 rows and 10 columns


# Create feature_matrix 
feature_matrix <- cov[, 4:ncol(cov)]


########################################## PCA ########################################## 


# Run PCA (returns a PCA object)
pca <- prcomp(feature_matrix, 
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


# Cluster by cutting the tree 
# (returns a series of cluster number)
cut_by_k <- cutree(Hierarchical_clustering, k = 8)
cut_by_h <- cutree(Hierarchical_clustering, h = 100)

print(cut_by_k)
print(cut_by_h)

# Combine the clustering result with the (x, y, z) coordinates
# (returns a data frame)
pca_coord$hcluster_k <- factor(cut_by_k)
pca_coord$hcluster_h <- factor(cut_by_h)

print(pca_coord)
        
# Clean the data frame before plotting
pca_coord_cleaned <- gather(pca_coord,
                            hclustering_by, 
                            hcluster_number, 
                            c("hcluster_k", "hcluster_h"))

print(pca_coord_cleaned)

# Plot hierarchical clustering results 
hierarchical_plot1 <- ggplot(pca_coord_cleaned,
       aes(x = X,
           y = Y, 
           color = hcluster_number)) +
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid(~ hclustering_by) + 
        labs(title = "Hierarchical clustering (k = 8, h = 100)",
             x = "PC1 (28%)",
             y = "PC2 (19%)")

plot(hierarchical_plot1)


hierarchical_plot2 <- ggplot(pca_coord_cleaned,
                            aes(x = X, 
                                y = Y, 
                                color = ICU, 
                                shape = gender)) + 
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid(hcluster_number ~ hclustering_by) + 
        labs(title = "Hierarchical clustering (k = 8, h = 100)",
             x = "PC1 (28%)",
             y = "PC2 (19%)")

plot(hierarchical_plot2)

# Recall original PCA plot
print(pca_2D_plot)




########################################## K-means Clustering ##########################################

# Determine optimal k (= number of clusters) by calculating total within cluster 
# sum of squares set.seed(32)
ttWithinss <- map_dbl(1:10, 
                      function(k) {
                              set.seed(32)
                              km <- kmeans(x = coord_matrix, 
                                           centers = k,
                                           nstart = 25)
                              km$tot.withinss})

print(ttWithinss)

# Create a scree plot 
ttWithinss_plot <- data.frame(tt_within_ss = ttWithinss,
                              k = factor(1:10)) %>%
        ggplot(aes(x = k, 
                   y = tt_within_ss,
                   group = 1)) +
        geom_point() + 
        geom_line() + 
        geom_vline(xintercept = 3, color = "red") + 
        labs(title = "Total Within-Cluster Sum of Squares Change",
             y = "Total Within-Cluster Sum of Squares",
             x = "Number of Clusters (k)")

print(ttWithinss_plot)

# Run k-means clustering with k = 3 
# (returns a kmeans object)
set.seed(32)
kmc <- kmeans(x = coord_matrix, 
              centers = 3,
              nstart = 25)


# Extract cluster results from the kmeans object
# (returns a series of cluster numbers)
kmc_cluster <- factor(kmc$cluster)


# Combine the kmeans clustering result with the pca coordinate table 
pca_coord$kmcluster <- kmc_cluster



# Plot the result of k-means clustering
kmeans_plot1 <- ggplot(pca_coord,
                       aes(x = X,
                           y = Y, 
                           color = kmcluster)) + 
        geom_point(size = 2, alpha = 0.5) + 
        labs(title = "K-Means Clustering", 
             x = "PC1 (28%)",
             y = "PC2 (19%)")

print(kmeans_plot1)


# Visualize characteristics of each cluster
kmeans_plot2 <- ggplot(pca_coord,
                       aes(x = X,
                           y = Y, 
                           color = ICU, 
                           shape = gender)) + 
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid(~ kmcluster) + 
        labs(title = "K-Means Clustering", 
             x = "PC1 (28%)",
             y = "PC2 (19%)")

print(kmeans_plot2)


print(pca_2D_plot)


########################################## t-SNE ##########################################


# Load the Rtsne package 
library(Rtsne)


# Set a function for running t-SNE 
tSNE_fn <- function(pp) {
        
        # tSNE
        set.seed(277)
        Rtsne(as.matrix(feature_matrix), 
              PCA = T,                 # Preliminary PCA
              perplexity = pp,         # Perplexity
              max_iter = 2000,         # Max iteration number  
              dims = 2)                # Number of output dimensions 
}
tsne1 <- tSNE_fn(1)
tsne2 <- tSNE_fn(2)
tsne5 <- tSNE_fn(5)
tsne10 <- tSNE_fn(10)

# Set a function for cleaning data 
data_clean_fn <- function(tsne_object, perplexity) { 
        
        data.frame(X = tsne_object$Y[, 1],
                   Y = tsne_object$Y[, 2],
                   gender = pca_coord$gender, 
                   ICU = pca_coord$ICU,
                   Perplexity = factor(perplexity)) 
}

# Clean your data
# (returns a data frame)
tsne_compare_df <- rbind(data_clean_fn(tsne1, 1),
                         data_clean_fn(tsne2, 2),
                         data_clean_fn(tsne5, 5),
                         data_clean_fn(tsne10, 10))

# Check out the relationship btw perplexity and the coordinates 
perplexity_plot <- ggplot(tsne_compare_df,
                          aes(x = X, 
                              y = Y, 
                              color = gender, 
                              shape = ICU)) + 
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid(~ Perplexity) + 
        labs(title = "t-SNE with Perplexity = 1, 2, 5, and 10", 
             x = "Dim1",
             y = "Dim2")

print(perplexity_plot)

# Subset with Perplexity = 2
tsne_pp <- filter(tsne_compare_df, Perplexity == 5)
rownames(tsne_pp) <- rownames(meta)

pheatmap(tsne_pp[, 1:2],
         clustering_distance_rows = "euclidean",
         annotation_row = meta,
         fontsize_row = 5,
         main = "Heatmap")


# Calculate distance: (X, Y, Z) coordinates
# (returns a distance matrix)
distance_tsne <- dist(tsne_pp[, 1:2], 
                 method = "euclidean")

# Perform hierarchical clustering 
Hierarchical_clustering_tsne <- hclust(distance, 
                                  method = "average")

# Create a dendrogram
plot(Hierarchical_clustering_tsne)

# Extract the clustering result (8 clusters)
# and clean data
hcluster_tsne <- cutree(Hierarchical_clustering_tsne,
                        k = 8)
tsne_pp$hcluster <- factor(hcluster_tsne)


# Plotting tSNE/hierarchical clustering results 
tsne_hclustering1 <- ggplot(tsne_pp,
                           aes(x = X,
                               y = Y, 
                               color = hcluster)) + 
        geom_point(size = 2, alpha = 0.5)

print(tsne_hclustering1)

tsne_hclustering2 <- ggplot(tsne_pp,
                            aes(x = X,
                                y = Y, 
                                color = ICU,
                                shape = gender)) + 
        geom_point(size = 2, alpha = 0.5) + 
        facet_grid(Perplexity ~ hcluster)

print(tsne_hclustering2)