## Load Packages ####
library(data.table)
library(Seurat)
library(dplyr)

##-----------Read and Format Data---------------------------------####

## Read in data, separate into phenotypic and feature data tables ####
dat <- fread("data/GSE72056_melanoma_single_cell_revised_v2.txt")
edat <- as.matrix(dat[4:nrow(dat), 2:ncol(dat)])
fdat <- data.table(genes=unlist(dat[4:nrow(dat),1]))
fdat$genes[17008] <- paste0(fdat$genes[17008], "_1")
fdat$genes[23202] <- paste0(fdat$genes[23202], "_2")
fdat$genes[8027] <- paste0(fdat$genes[8027], "_1")
fdat$genes[23637] <- paste0(fdat$genes[23637], "_2")
pdat <- data.table(
  cellID=colnames(dat)[2:length(colnames(dat))],
  tumor=unlist(dat[1,2:ncol(dat)]),
  malignant=unlist(dat[2,2:ncol(dat)]), #malignant(1=no, 2=yes, 0=unresolved)
  cell_type=unlist(dat[3,2:ncol(dat)]) #non-malignant cell type (1=T, 2=B, 3=Macro., 4=Endo., 5=CAF; 6=NK)
)
rownames(edat) <- fdat$genes

malignant <- edat[,pdat$cellID[pdat$malignant == 2]]
nonMalignant <- edat[,pdat$cellID[pdat$malignant != 2]]

# does not include unresovled (0s)
nonmalignant <- edat[,pdat$cellID[pdat$malignant == 1]]

## Take a look at the data ####
dat[1:4, 1:4]
edat[1:3, 1:3]
head(fdat)
head(pdat)


##----------Analysis with Seurat----------------------------------####

## Create Seurat Object ####
sobjM <- CreateSeuratObject(raw.data = malignant)
sobjN <- CreateSeuratObject(raw.data = nonMalignant)

## Expression QC - visualize and plot gene and molecular counts ####
VlnPlot(sobjM, features.plot = c("nGene", "nUMI"))
GenePlot(sobjM, gene1 = "nGene", gene2 = "nUMI")
VlnPlot(sobjN, features.plot = c("nGene", "nUMI"))
GenePlot(sobjN, gene1 = "nGene", gene2 = "nUMI")

## Normalization ####
sobjM <- NormalizeData(sobjM, normalization.method = 'LogNormalize', scale.factor = 10000)
sobjN <- NormalizeData(sobjN, normalization.method = 'LogNormalize', scale.factor = 10000)

## Finding Highly Variable Genes ###
sobjM <- FindVariableGenes(sobjM,
                          mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125,
                          x.high.cutoff = 3,
                          y.cutoff = 0.5)
sobjN <- FindVariableGenes(sobjN,
                          mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125,
                          x.high.cutoff = 3,
                          y.cutoff = 0.5)


## Scaling data to regress out confounders (detected molecules per cell) ####
sobjM <- ScaleData(sobjM, vars.to.regress = c('nUMI'))
sobjN <- ScaleData(sobjN, vars.to.regress = c('nUMI'))

## Dimensionality Reduction ####
sobjM <- RunPCA(sobjM, pc.genes = sobjM@var.genes, do.print = T, pcs.print = 1:5, genes.print = 5)
PCAPlot(sobjM, dim.1 = 1, dim.2 = 2)
PCHeatmap(sobjM, pc.use = 1:6, cells.use = 500, do.balanced = T, label.columns = F, use.full = F)
PCElbowPlot(sobjM)

sobjN <- RunPCA(sobjN, pc.genes = sobjN@var.genes, do.print = T, pcs.print = 1:5, genes.print = 5)
PCAPlot(sobjN, dim.1 = 1, dim.2 = 2)
PCHeatmap(sobjN, pc.use = 1:6, cells.use = 500, do.balanced = T, label.columns = F, use.full = F)
PCElbowPlot(sobjN)


## Clustering and tSNE Plot ####
sobjM <- FindClusters(sobjM, reduction.type = "pca", dims.use = 1:8, resolution = 1.0, print.output = 0, save.SNN = T)
sobjM <- RunTSNE(sobjM, dims.use = 1:8)
TSNEPlot(sobjM)

sobjN <- FindClusters(sobjN, reduction.type = "pca", dims.use = 1:8, resolution = 1.0, print.output = 0, save.SNN = T)
sobjN <- RunTSNE(sobjN, dims.use = 1:8)
TSNEPlot(sobjN)

## Find Marker Genes ####
# markers2 <- FindMarkers(sobj, 2)
# VlnPlot(object = sobj, features.plot = rownames(markers2)[1:6])
# FeaturePlot(
#   sobj, 
#   head(rownames(markers2)), 
#   cols.use = c("lightgrey", "blue"), 
#   nCol = 3
# )

# markers <- FindAllMarkers(
#   object = sobj,
#   only.pos = TRUE,
#   min.pct = 0.25,
#   thresh.use = 0.25
# )
# top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# DoHeatmap(
#   object = sobj,
#   genes.use = top10$gene,
#   slim.col.label = TRUE,
#   remove.key = TRUE
# )
# FeaturePlot(
#   sobj,
#   head(rownames(markers)),
#   cols.use = c("lightgrey", "blue"),
#   nCol = 3
# )
