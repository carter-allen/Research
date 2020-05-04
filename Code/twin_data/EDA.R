# packages
library(edgeR)
library(tidyverse)
library(pheatmap)

# Read count data
counts <- read.table("Data/Disease_vs_Control_Counts.txt", 
                     header = TRUE, 
                     sep = '\t')
counts_mat <- counts[,-1] # remove row name row
genes <- as.character(counts[,1]) # pull genes

# Read group data
groups_df <- read.table("Data/Disease_vs_Control_groups.txt",
                       header = TRUE,
                       sep = '\t')
# get groups
groups <- groups_df %>%
    pull(condition) %>%
    as.character()

# Get sample order
samples <- groups_df %>%
    pull(samples) %>%
    as.character()

counts_mat <- counts_mat[,samples]

colnames(counts_mat) <- paste(colnames(counts_mat),
                              groups,
                              sep = "_")

# Form DGEList
y <- DGEList(counts = counts_mat,
             genes = genes,
             group = groups)

# Order and remove duplicates
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$genes)
y <- y[!d,]

# Normalization
y <- calcNormFactors(y)

# Heatmap of correlation among samples
cor_plot <- pheatmap(cor(cpm(y)),
                     cluster_rows = F,
                     cluster_cols = F)
cor_plot

plotMDS(y)

# Plots of filtered data
n.pre <- nrow(y)
n.pre
y <- y[rowSums(cpm(y) > 1) >= 3,,keep.lib.sizes = FALSE]
n.post <- nrow(y)
n.post
mds_plot_filtered <- plotMDS(y)
cor_plot_filtered <- pheatmap(cor(cpm(y)),
                              cluster_rows = F,
                              cluster_cols = F)

y.top50 <- y[1:50,]
o <- order(rowSums(cpm(y.top50)[,4:6]+1) / rowSums(cpm(y.top50)[,1:3]+1))
p50 <- pheatmap(log10(cpm(y.top50)+1)[o,],
                scale="row", 
                cluster_rows=FALSE, 
                cluster_cols=FALSE, 
                show_rownames=TRUE)
p50

y.top100 <- y[1:100,]
o <- order(rowSums(cpm(y.top100)[,4:6]+1) / rowSums(cpm(y.top100)[,1:3]+1))
p100 <- pheatmap(log10(cpm(y.top100)+1)[o,],
                 scale="row", 
                 cluster_rows=FALSE, 
                 cluster_cols=FALSE, 
                 show_rownames=TRUE)
