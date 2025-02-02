# Script to process a single-cell RNA-Seq dataset and approximate
# a bulk RNA-Seq experiment dataset
# Goal: after this script, you can run through the scripts we have used
# in class for analyzing RNA-Seq
# Should pick up at step 3 script after this script is done

library(Seurat)
library(SeuratData)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(plotly)
library(ggrepel)
SeuratData::InstallData("pbmc3k")
library(enrichR)
library(escape)
library(enrichplot)

# set working directory to source file location
md <- read.table(gzfile("Data/GSE218743_irae_ht_integrated_metadata.csv.gz"), header = TRUE, sep = ",")
counts <- read.table(gzfile("Data/GSE218743_irae_ht_integrated_counts.csv.gz"), header = TRUE, sep = ",", row.names = 1)

seu.obj <- CreateSeuratObject(counts = counts, meta.data = md)

# filter out cells that had less than 200 genes or more than 2500 genes
seu.filtered <- subset(seu.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# look at the object and see how many cells were filtered out
seu.obj
seu.filtered

# Run Seurat workflow - need to see if the data I'm using has already had these done
seu.filtered <- NormalizeData(seu.filtered)
seu.filtered <- FindVariableFeatures(seu.filtered)
seu.filtered <- ScaleData(seu.filtered)

seu.filtered <- RunPCA(seu.filtered)
ElbowPlot(seu.filtered)
seu.filtered <- RunUMAP(seu.filtered, dims = 1:20)

# visualize the data
cell_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'CellType', label = TRUE)
cond_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'condition')

cell_plot|cond_plot

# pseudo-bulk workflow
# we are going to aggregate the counts across cells to convert between
# the single-cell data and approximate a bulk RNA-Seq experiment
# Acquire necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

# make a new column in the metadata that contains both the condition and the sample
seu.filtered$samples <- paste0(seu.filtered$condition, seu.filtered$sample)

# aggregate by cell type and sample ID
cts <- AggregateExpression(seu.filtered,
                    group.by = c("CellType", "samples"),
                    assays = "RNA",
                    slot = "counts",
                    return.seurat = FALSE)
 
# extract just the count data
cts <- cts$RNA

# transpose so that columns are genes and rows are cells
cts.t <- t(cts)
cts.t <- as.data.frame(cts.t)

# get values for where to split into cells
splitRows <- gsub('_.*', '', rownames(cts.t))

# split the data frame
cts.split <- split.data.frame(cts.t, f = factor(splitRows))

# fix column names by removing the cell types
# transpose so that genes are rows again
# the gsub expression in the function will match
# the characters after the underscore and replace
# the rowname with just the treatment type of sample number
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

# Now run differential expression analysis with DESeq2
# this will run through the analysis for just one of the cell types
# 1. Get counts matrix
counts_cell1 <- cts.split.modified$'g4'

# 2. Generate sample level metadata
colData <- read.csv("colData.csv")



#colData <- column_to_rownames(var = 'samples')



# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_cell1,
                       colData = colData,
                       design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)

# check the coefficients for the comparison
resultsNames(dds)

# generate results object
res <- results(dds, name = "condition_Hashimoto_vs_Control")
res
resOrdered <- res[order(res$pvalue),]
summary(res)

#DESeq has mechanisms for the step visualizations


# VISUALIZATIONS (MOST OF RIA'S EDITS)

#dat <- bind_rows(cts.split.modified)

#RidgePlot(cts.split.modified, features = features, ncol = 2)

ens_gene (label)
pval (y)

EnhancedVolcano(res,
            lab = rownames(res),
            x = 'log2FoldChange',
            y = 'padj',
            pCutoff = 0.05,
            labSize = 5)

dat1 <- as_tibble(res, rownames = "rownames")

#filteredparams <- res %>%
#  Filter('log2FoldChange' > 0.05)

#(plot) +
#  coord_cartesian(xlim(10, 27))

features <- c("HBG2", "IGLVI-70", "IGKV3-20", "IGHA1", "IGHV4-34", "IGKV1-9", "IGKV1-17", "IGKV3-15")
nums <- seu.filtered$CellType
seu.filtered$orig.ident <- nums

DoHeatmap(seu.filtered, features = features, group.by = "CellType",
          col.use = my_cols2)

# possible alternative to DoHeatMap, could not get working
#my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
#scales::show_col(my_cols2)
#DimHeatmap(seu.filtered, nfeatures = 100, dims = 1:13, cells = 500, balanced = TRUE, col.use = "red")

RidgePlot(seu.filtered, features = "HBG2", group.by = "orig.ident") # could not get completely right
VlnPlot(seu.filtered, features = features, group.by = "orig.ident") # should come with colored density metrics in violin shape, but could not get completely right
FeaturePlot(seu.filtered, features = features) # shows where certain genes are expressed, comparable with PCA-like plot from earlier
DotPlot(seu.filtered, features = features, group.by = "orig.ident") + RotatedAxis() # expression levels of different genes in different cells

# everything below is my attempt to make enrichment plots, could not get working

#dbs <- listEnrichrDbs() %>%
#  arrange(-geneCoverage)

#DEenrichRPlot(seu.filtered, balanced = TRUE, max.genes = 100, test.use = "t", enrich.database = "Enrichr_Users_Contributed_Lists_2020")

#gseaplot(colData, geneSetID = 1, by = "runningScore", title = edo2$Description[1])


#DEenrichRPlot(
#  seu.filtered,
#  ident.1 = colData,
#  ident.2 = barsCtrl,
#  balanced = T,
#  logfc.threshold = 0.25,
#  assay = NULL,
#  max.genes = 100,
#  test.use = "wilcox",
#  p.val.cutoff = 0.05,
#  cols = NULL,
#  enrich.database = "GO_Molecular_Function_2023",
#  num.pathway = 30,
#  return.gene.list = F)


