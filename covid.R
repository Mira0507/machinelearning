library(ggplot2)
library(tidyverse)
library(pheatmap)

# Import RNA-seq data (returns a data frame)
cov <- t(read_tsv("GSE157103_genes.ec.tsv"))
cov1 <- apply(cov[2:nrow(cov), 2:ncol(cov)], 2, as.numeric)
colnames(cov1) <- cov[1, 2:ncol(cov)]
rownames(cov1) <- rownames(cov)[2:nrow(cov)]
cov2 <- cov1[, colSums(cov1) >= nrow(cov1)]

cov3 <- as.data.frame(cov2) %>% 
        rownames_to_column(var = "Sample")

# Explore your data frame 
print(dim(cov3))   # number of rows and columns
print(cov3[1:10, 1:10])  # first 10 rows and 10 columns


write.csv(cov3, "covid.csv")


# Run PCA (returns a PCA object)
pca <- prcomp(ser[, 2:ncol(ser)], 
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
pca_coord <- data.frame(Sample = ser$Sample,
                        X = pca$x[, 1],
                        Y = pca$x[, 2],
                        Z = pca$x[, 3]) %>% 
        mutate(Genotype = ifelse(str_detect(Sample, "WT"), "WT", "KO")) %>%
        separate(Sample, c("Sample", "v1", "v2"), sep = "_") %>%
        select(-v1, -v2) %>% 
        mutate(time_point = case_when(str_detect(Sample, "8") ~ "8",
                                      str_detect(Sample, "14") ~ "14",
                                      str_detect(Sample, "21") ~ "21",
                                      str_detect(Sample, "35") ~ "35",
                                      str_detect(Sample, "70") ~ "70"),
               time_point = factor(time_point,
                                   levels = c("8", 
                                              "14",
                                              "21",
                                              "35",
                                              "70")))
print(pca_coord)


# Plot the PCA result in 2D space 
ggplot(pca_coord,
       aes(x = X, 
           y = Y, 
           color = time_point,
           shape = Genotype)) + 
        geom_point(size = 2, alpha = 0.5) +
        labs(title = "PCA",
             x = "PC1 (49%)",
             y = "PC2 (46%)")
