## Load Packages ####
library(data.table)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scran)

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
  mal=unlist(dat[2,2:ncol(dat)]), #malignant(1=no, 2=yes, 0=unresolved)
  cell_type=unlist(dat[3,2:ncol(dat)]) #non-malignant cell type (1=T, 2=B, 3=Macro., 4=Endo., 5=CAF; 6=NK)
)
rownames(edat) <- fdat$genes

## Tumors with at least 100 cells
filteredTumors1 <- pdat[,.N, by = tumor] %>%
  filter(., N >= 100) %>%
  select(., tumor) %>%
  pull(., tumor)

## Tumors with >50 malignant cells
filteredTumors2 <- pdat[mal ==2 & tumor %in% filteredTumors1, .N, by = tumor] %>%
  filter(., N > 50) %>%
  select(., tumor) %>%
  pull(., tumor)

## Define malignant and non-Malignant Cells ####
nonMalignant <- edat[,pdat[tumor %in% filteredTumors1 & mal == 1, cellID]]
malignant <- edat[,pdat[tumor %in% filteredTumors2 & mal == 2, cellID]]

## Take a look at the data ####
dat[1:4, 1:4]
edat[1:3, 1:3]
head(fdat)
head(pdat)
factor(pdat$tumor)

factor(dat[1,])

##----------Analysis with Seurat----------------------------------####

## For Convenience these steps can be condensed into a single function ####
run_seurat <- function(exprsMatrix){
  ## Create Seurat Object ####
  sobj <- CreateSeuratObject(raw.data = exprsMatrix)
  ## Normalization ####
  sobj <- NormalizeData(sobj, normalization.method = 'LogNormalize', scale.factor = 10000)
  ## Finding Highly Variable Genes ###
  sobj <- FindVariableGenes(sobj,
                            mean.function = ExpMean,
                            dispersion.function = LogVMR,
                            x.low.cutoff = 0.0125,
                            x.high.cutoff = 3,
                            y.cutoff = 0.5)
  ## Scaling data to regress out confounders (detected molecules per cell) ####
  sobj <- ScaleData(sobj, vars.to.regress = c('nUMI'))
  ## Dimensionality Reduction ####
  sobj <- RunPCA(sobj, pc.genes = sobj@var.genes, do.print = T, pcs.print = 1:5, genes.print = 5)
  ## Clustering and tSNE Plot ####
  sobj <- FindClusters(sobj, reduction.type = "pca", dims.use = 1:8, resolution = 1.0, print.output = 0, save.SNN = T)
  sobj <- RunTSNE(sobj, dims.use = 1:8)
  TSNEPlot(sobj)
  return(sobj)
}

## Create Seurat Object to Assess QC ####
sobjM <- CreateSeuratObject(raw.data = malignant)
sobjN <- CreateSeuratObject(raw.data = nonMalignant)

## Before Expression QC - visualize and plot gene and molecular counts ####
VlnPlot(sobjM, features.plot = c("nGene", "nUMI"))
GenePlot(sobjM, gene1 = "nGene", gene2 = "nUMI")
VlnPlot(sobjN, features.plot = c("nGene", "nUMI"))
GenePlot(sobjN, gene1 = "nGene", gene2 = "nUMI")

## Run Seurat Pipeline ####
sobjM <- run_seurat(malignant)
sobjN <- run_seurat(nonMalignant)
sobjC <- run_seurat(edat)

## After Expression QC - visualize and plot gene and molecular counts ####
VlnPlot(sobjM, features.plot = c("nGene", "nUMI"))
GenePlot(sobjM, gene1 = "nGene", gene2 = "nUMI")
VlnPlot(sobjN, features.plot = c("nGene", "nUMI"))
GenePlot(sobjN, gene1 = "nGene", gene2 = "nUMI")

## Dimensionality Reduction Plots ####
PCAPlot(sobjM, dim.1 = 1, dim.2 = 2)
PCHeatmap(sobjM, pc.use = 1:6, cells.use = 500, do.balanced = T, label.columns = F, use.full = F)
PCElbowPlot(sobjM)
PCAPlot(sobjN, dim.1 = 1, dim.2 = 2)
PCHeatmap(sobjN, pc.use = 1:6, cells.use = 500, do.balanced = T, label.columns = F, use.full = F)
PCElbowPlot(sobjN)

## Color TSNE by tumor ####
par(mfrow=c(1,2))
plot(sobjM@dr$tsne@cell.embeddings, col=factor(pdat$tumor[pdat$cellID %in% colnames(malignant)]), main="Malignant Cells")
plot(sobjN@dr$tsne@cell.embeddings, col=factor(pdat$tumor[pdat$cellID %in% colnames(nonMalignant)]), main="Non-Malignant Cells")
mtext("Colored by Tumor", side = 3, line = -1, outer = TRUE)
par(mfrow=c(1,1))

## Extract tSNE Data for Malignant Cells ####
mtSNE <- data.table(
  cellID = rownames(sobjM@dr$tsne@cell.embeddings),
  tSNE_1 = sobjM@dr$tsne@cell.embeddings[,1],
  tSNE_2 = sobjM@dr$tsne@cell.embeddings[,2],
  tumor = factor(pdat$tumor[pdat$cellID %in% colnames(malignant)]),
  cell_type = factor(pdat$cell_type[pdat$cellID %in% colnames(malignant)])
)

## Extract tSNE Data for non-Malignant Cells ####
ntSNE <- data.table(
  cellID = rownames(sobjN@dr$tsne@cell.embeddings),
  tSNE_1 = sobjN@dr$tsne@cell.embeddings[,1],
  tSNE_2 = sobjN@dr$tsne@cell.embeddings[,2],
  tumor = factor(pdat$tumor[pdat$cellID %in% colnames(nonMalignant)]),
  cell_type = factor(pdat$cell_type[pdat$cellID %in% colnames(nonMalignant)])
)

## Extract tSNE Data for All Cells ####
ctSNE <- data.table(
  cellID = rownames(sobjC@dr$tsne@cell.embeddings),
  tSNE_1 = sobjC@dr$tsne@cell.embeddings[,1],
  tSNE_2 = sobjC@dr$tsne@cell.embeddings[,2],
  tumor = factor(pdat$tumor[pdat$cellID %in% colnames(edat)]),
  cell_type = factor(pdat$cell_type[pdat$cellID %in% colnames(edat)])
)

## tSNE Plots - colored by tumor and cell type ####
plot_grid(
  ggplot(data = mtSNE, aes(x = tSNE_1, y = tSNE_2, col=tumor))+
    ggtitle(label = "Malignant Cells")+
    scale_color_manual(
      labels = gsub("^", "Mel", levels(mtSNE$tumor)),
      values = brewer.pal(6, "Dark2"))+
    geom_point(),
  
  ggplot(data = ntSNE, aes(x = tSNE_1, y = tSNE_2, col=cell_type))+
    ggtitle(label = "Non-Malignant Cells")+
    scale_color_manual(
      labels = c("Malignant", "T-Cells", "B-Cells", "Macrophages", "Endothelial Cells", "CAF", "NK"), 
      values = brewer.pal(7, "Dark2"))+
    geom_point()
)
# ggsave("plots/tSNE_tumor_celltype.pdf", width = 15, height = 6)

## tSNE Plots - colored by cell type
plot_grid(
  ggplot(data = mtSNE, aes(x = tSNE_1, y = tSNE_2, col=cell_type))+
    ggtitle(label = "Malignant Cells")+
    scale_color_manual(
      labels = c("Malignant", "T-Cells", "Macrophages"), 
      values = brewer.pal(3, "Dark2"))+
    geom_point(),
  
  ggplot(data = ntSNE, aes(x = tSNE_1, y = tSNE_2, col=cell_type))+
    ggtitle(label = "Non-Malignant Cells")+
    scale_color_manual(
      labels = c("Malignant", "T-Cells", "B-Cells", "Macrophages", "Endothelial Cells", "CAF", "NK"), 
      values = brewer.pal(7, "Dark2"))+
    geom_point()
)
# ggsave("plots/tSNE_celltype.pdf", width = 15, height = 6)

## tSNE Plots - colored by tumor
plot_grid(
  ggplot(data = mtSNE, aes(x = tSNE_1, y = tSNE_2, col=tumor))+
    ggtitle(label = "Malignant Cells")+
    geom_point(),
  
  ggplot(data = ntSNE, aes(x = tSNE_1, y = tSNE_2, col=tumor))+
    ggtitle(label = "Non-Malignant Cells")+
    geom_point()
)
# ggsave("plots/tSNE_tumor.pdf", width = 15, height = 6)

## Combined tSNE Plot
plot_grid(
  ggplot(data = ctSNE, aes(x = tSNE_1, y = tSNE_2, col=tumor))+
    ggtitle(label = "All Cells")+
    geom_point(),
  
  ggplot(data = ctSNE, aes(x = tSNE_1, y = tSNE_2, col=cell_type))+
    ggtitle(label = "All Cells")+
    scale_color_manual(
      labels = c("Malignant", "T-Cells", "B-Cells", "Macrophages", "Endothelial Cells", "CAF", "NK"), 
      values = brewer.pal(7, "Dark2"))+
    geom_point()
)
# ggsave("plots/combined_tSNE_tumor_celltype.pdf", width = 15, height = 6)

## Identifying Clusters with DBScan ####
kNNdistplot(sobjN@dr$tsne@cell.embeddings, k=10)
abline(h=2)
db <- fpc::dbscan(sobjN@dr$tsne@cell.embeddings, eps = 2, MinPts = 10)

db_plot <- ggplot(ntSNE, aes(x = tSNE_1, y = tSNE_2, col = factor(db$cluster)))+
  ggtitle(label = "After Clustering with DBSCAN")+
  geom_point()

plot_grid(
  ggplot(ntSNE, aes(x = tSNE_1, y = tSNE_2))+
    ggtitle(label = "Before Clustering")+
    geom_point(),
  
  ggplot(ntSNE, aes(x = tSNE_1, y = tSNE_2, col = factor(db$cluster)))+
    ggtitle(label = "After Clustering with DBSCAN")+
    geom_point()
)
# ggsave("plots/dbscan_clusters.pdf", width = 15, height = 6)

## Compile Gene Signature Markers and Make Feature Plots ####
markers <- list(
  tcells=c('CD2','CD3D','CD3E','CD3G','CD8A','SIRPG','TIGIT','GZMK','ITK','SH2D1A','CD247','PRF1','NKG7',
           'IL2RB','SH2D2A','KLRK1','ZAP70','CD7','CST7','LAT','PYHIN1','SLA2','STAT4','CD6','CCL5','CD96',
           'TC2N','FYN','LCK','TCF7','TOX','IL32','SPOCK2','SKAP1','CD28','CBLB','APOBEC3G','PRDM1'),
  bcells=c("CD19", "CD79A", "CD79B", "BLK", "MS4A1", "BANK1", "IGLL3P", 
           "FCRL1", "PAX5", "CLEC17A", "CD22", "BCL11A", "VPREB3", "HLA-DOB", 
           "STAP1", "FAM129C", "TLR10", "RALGPS2", "AFF3", "POU2AF1", "CXCR5", 
           "PLCG2", "HVCN1", "CCR6", "P2RX5", "BLNK", "KIAA0226L", "POU2F2", 
           "IRF8", "FCRLA", "CD37"),
  macro=c("CD163", "CD14", "CSF1R", "C1QC", "VSIG4", "C1QA", "FCER1G", 
          "F13A1", "TYROBP", "MSR1", "C1QB", "MS4A4A", "FPR1", "S100A9", 
          "IGSF6", "LILRB4", "FPR3", "SIGLEC1", "LILRA1", "LYZ", "HK3", 
          "SLC11A1", "CSF3R", "CD300E", "PILRA", "FCGR3A", "AIF1", "SIGLEC9", 
          "FCGR1C", "OLR1", "TLR2", "LILRB2", "C5AR1", "FCGR1A", "MS4A6A", 
          "C3AR1", "HCK", "IL4I1", "LST1", "LILRA5", "CSTA", "IFI30", "CD68", 
          "TBXAS1", "FCGR1B", "LILRA6", "CXCL16", "NCF2", "RAB20", "MS4A7", 
          "NLRP3", "LRRC25", "ADAP2", "SPP1", "CCR1", "TNFSF13", "RASSF4", 
          "SERPINA1", "MAFB", "IL18", "FGL2", "SIRPB1", "CLEC4A", "MNDA", 
          "FCGR2A", "CLEC7A", "SLAMF8", "SLC7A7", "ITGAX", "BCL2A1", "PLAUR", 
          "SLCO2B1", "PLBD1", "APOC1", "RNF144B", "SLC31A2", "PTAFR", "NINJ1", 
          "ITGAM", "CPVL", "PLIN2", "C1orf162", "FTL", "LIPA", "CD86", 
          "GLUL", "FGR", "GK", "TYMP", "GPX1", "NPL", "ACSL1"),
  endo=c("PECAM1", "VWF", "CDH5", "CLDN5", "PLVAP", "ECSCR", "SLCO2A1", 
         "CCL14", "MMRN1", "MYCT1", "KDR", "TM4SF18", "TIE1", "ERG", "FABP4", 
         "SDPR", "HYAL2", "FLT4", "EGFL7", "ESAM", "CXorf36", "TEK", "TSPAN18", 
         "EMCN", "MMRN2", "ELTD1", "PDE2A", "NOS3", "ROBO4", "APOLD1", 
         "PTPRB", "RHOJ", "RAMP2", "GPR116", "F2RL3", "JUP", "CCBP2", 
         "GPR146", "RGS16", "TSPAN7", "RAMP3", "PLA2G4C", "TGM2", "LDB2", 
         "PRCP", "ID1", "SMAD1", "AFAP1L1", "ELK3", "ANGPT2", "LYVE1", 
         "ARHGAP29", "IL3RA", "ADCY4", "TFPI", "TNFAIP1", "SYT15", "DYSF", 
         "PODXL", "SEMA3A", "DOCK9", "F8", "NPDC1", "TSPAN15", "CD34", 
         "THBD", "ITGB4", "RASA4", "COL4A1", "ECE1", "GFOD2", "EFNA1", 
         "PVRL2", "GNG11", "HERC2P2", "MALL", "HERC2P9", "PPM1F", "PKP4", 
         "LIMS3", "CD9", "RAI14", "ZNF521", "RGL2", "HSPG2", "TGFBR2", 
         "RBP1", "FXYD6", "MATN2", "S1PR1", "PIEZO1", "PDGFA", "ADAM15", 
         "HAPLN3", "APP"),
  caf=c("FAP", "THY1", "DCN", "COL1A1", "COL1A2", "COL6A1", "COL6A2", 
        "COL6A3", "CXCL14", "LUM", "COL3A1", "DPT", "ISLR", "PODN", "CD248", 
        "FGF7", "MXRA8", "PDGFRL", "COL14A1", "MFAP5", "MEG3", "SULF1", 
        "AOX1", "SVEP1", "LPAR1", "PDGFRB", "TAGLN", "IGFBP6", "FBLN1", 
        "CA12", "SPOCK1", "TPM2", "THBS2", "FBLN5", "TMEM119", "ADAM33", 
        "PRRX1", "PCOLCE", "IGF2", "GFPT2", "PDGFRA", "CRISPLD2", "CPE", 
        "F3", "MFAP4", "C1S", "PTGIS", "LOX", "CYP1B1", "CLDN11", "SERPINF1", 
        "OLFML3", "COL5A2", "ACTA2", "MSC", "VASN", "ABI3BP", "C1R", 
        "ANTXR1", "MGST1", "C3", "PALLD", "FBN1", "CPXM1", "CYBRD1", 
        "IGFBP5", "PRELP", "PAPSS2", "MMP2", "CKAP4", "CCDC80", "ADAMTS2", 
        "TPM1", "PCSK5", "ELN", "CXCL12", "OLFML2B", "PLAC9", "RCN3", 
        "LTBP2", "NID2", "SCARA3", "AMOTL2", "TPST1", "MIR100HG", "CTGF", 
        "RARRES2", "FHL2"),
  
  ## Natural Killer gene signature from BioRxiv
  ## https://www.biorxiv.org/content/biorxiv/suppl/2018/07/23/375253.DC1/375253-1.pdf
  ## https://www.biorxiv.org/content/10.1101/375253v2
  
  # nk=c("CCL5", "CD2", "CD244", "CD247", "CD7", "CRTAM", "FASLG", "GZMA", 
  #      "GZMB", "GZMH", "GZMK", "GZMM", "IKZF3", "IL18RAP", "IL2RB", 
  #      "KIR2DL1", "KIR2DL4", "KIR2DS4", "KIR3DL1", "KIR3DL2", "KLRB1", 
  #      "KLRC1", "KLRC3", "KLRD1", "KLRK1", "LCK", "LTA", "NCR1", "NCR3", 
  #      "NKG7", "PRF1", "PTPRCAP", "PYHIN1", "SAMD3", "SH2D1A", "TBX21", 
  #      "TNFSF14", "XCL1", "XCL2", "ZAP70")
  
  nk=c("KLRB1", "KLRC1", "KLRD1", "KLRF1"),
  cd8T=c("CD8A", "CD8B")
)

temp <- nonMalignant
tcell_sig <- sapply(1:ncol(temp), function(x) temp[which(fdat$genes %in% markers$tcells),x] %>% mean)
bcell_sig <- sapply(1:ncol(temp), function(x) temp[which(fdat$genes %in% markers$bcells),x] %>% mean)
macro_sig <- sapply(1:ncol(temp), function(x) temp[which(fdat$genes %in% markers$macro),x] %>% mean)
endo_sig <- sapply(1:ncol(temp), function(x) temp[which(fdat$genes %in% markers$endo),x] %>% mean)
caf_sig <- sapply(1:ncol(temp), function(x) temp[which(fdat$genes %in% markers$caf),x] %>% mean)
nk_sig <- sapply(1:ncol(temp), function(x) temp[which(fdat$genes %in% markers$nk),x] %>% mean)
cd8T_sig <- sapply(1:ncol(temp), function(x) temp[which(fdat$genes %in% markers$cd8T),x] %>% mean)
temp <- rbind(temp, tcell_sig, bcell_sig, macro_sig, endo_sig, caf_sig, nk_sig, cd8T_sig)
rownames(temp)[nrow(temp)]
## Run tSNE Pipeline ####
temp_obj <- run_seurat(temp)

## Feature Plots of 'Cell Markers' ####
p <- FeaturePlot(temp_obj, c("tcell_sig", "bcell_sig", "macro_sig", "endo_sig", "caf_sig", "nk_sig", "cd8T_sig"),
                 cols.use = c("lightgrey", "blue"), nCol = 4, no.legend = F, do.return = T)

p$tcell_sig <- p$tcell_sig + ggtitle("T-Cells")
p$bcell_sig <- p$bcell_sig + ggtitle("B-Cells")
p$macro_sig <- p$macro_sig + ggtitle("Macrophages")
p$endo_sig <- p$endo_sig + ggtitle("Endothelial Cells")
p$caf_sig <- p$caf_sig + ggtitle("CAFs")
p$nk_sig <- p$nk_sig + ggtitle("NK Cells")
p$cd8T_sig <- p$cd8T_sig + ggtitle("CD8T Cells")

## Feature Plots ####
plot_grid(p$tcell_sig, p$bcell_sig, p$macro_sig, p$endo_sig, p$caf_sig, p$nk_sig, p$cd8T_sig, db_plot, ncol = 4)
# ggsave("plots/feature_plots.pdf", width = 18, height = 8)

