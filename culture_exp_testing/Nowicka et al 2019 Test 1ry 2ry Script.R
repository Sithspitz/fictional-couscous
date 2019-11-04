############# UMAP and FlowSOM ###############
### Using Nowicka et al 2019 Method ###

library(flowCore)
library(HDCytoData)
library(readxl)
library(CATALYST)
library(ggplot2)

setwd("C:/Users/rbuch/DataShare/TSNE_A549_H2110_H460/H2110_H460_301019/Lymphocyte_and_ccr6_gated/Nowicka_test")

# Import FCS
a <- read.FCS("H2110_1ry_2ry_301019.fcs")
b <- read.FCS("H460_1ry_2ry_301019.fcs")
c <- read.FCS("RPMI_1ry_2ry_301019.fcs")

# Combine flowFrames into a flowSet
set <- flowSet(a, b, c)

# Load experiment metadata and paneldata
md <- read_excel("metadata_1ry_2ry.xlsx")
pd <- read_excel("paneldata_1ry_2ry.xlsx")

# Check that all panel columns are in the flowSet object
all(pd$fcs_colname %in% colnames(set)) 

# Specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("ctrl", "cancer"))
md$sample_id <- factor(md$sample_id,
                       levels = md$sample_id[order(md$condition)])

# Subset paneldata to select which markers looked at
## Here am removing Comp-APC-CY7-A
pd2 <- pd[1:7, ]

# construct daFrame
daf <- daFrame(set, pd2, md, cols_to_use = pd2$fcs_colname)

# Check Counts acquired in each sample
plotCounts(daf, color_by = "condition")

# Plot Expression
expression_plot <- plotExprs(daf, color_by = "condition")
expression_plot$facet$params$ncol <- 3
expression_plot

# Visualise Similarities between Patients 
## Using Multidimensional Scaling and Heatmapping
MDS_plot <- plotMDS(daf, color_by = "condition")
MDS_plot
plotExprHeatmap(daf, bin_anno = TRUE, row_anno = TRUE)

# Can determine Which markers causing variance, using PCA-based non-redundancy scoring
## Markers with higher scores explain a large proportion of variability

# UMAP to determine the number of clusters to give to FlowSOM
## Can set rows_to_use to not include all cells
set.seed(7)
dim_reduc_data <- runDR(daf, "UMAP", rows_to_use = NULL)
plotDR(dim_reduc_data, "UMAP", color_by = "sample_id")


# FlowSOM using the cluster() function
## Only using Type markers (e.g. CDs) and not State markers (e.g. pSTAT3)
dim_clusters <- cluster(dim_reduc_data, cols_to_use = type_markers(dim_reduc_data),
               xdim = 10, ydim = 10, maxK = 20, seed = 7, verbose = T) 

# Overview Plots
plotDR(dim_clusters, "UMAP", color_by = "meta20")
plotClusterHeatmap(dim_clusters,
                   hm2 = NULL, k = "meta20", m = NULL,
                   cluster_anno = TRUE, draw_freqs = TRUE) 

# Looking at the metaclustering by TSNE and PCA
## Also a Heatmap script using plotClusterHeatmap() in the paper
plotCodes(dim_clusters, k = "meta20")

# If want a rightsided heatmap of a marker
plotClusterHeatmap(dim_clusters, hm2 = "CD4", k = "meta20", draw_freqs = TRUE)

# Plotting different UMAPs
plotDR(dim_clusters, "UMAP", color_by = "meta20", facet = "sample_id")
plotDR(dim_clusters, "UMAP", color_by = "meta20", facet = "condition")

######################################################################
########## Paper then goes on to merge scripts #######################
########## But I am not doing that this time #########################
######################################################################

