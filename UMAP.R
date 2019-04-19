library(data.table)
library(Seurat)
library(sctransform)
library(umap)
library(ggplot2)

#-----------------------------------------------------------------------------#
# Using Eric's manipulated data...
# malignant, nonmalignant, pdat

tumor_M <- pdat[malignant==2, tumor]
tumor_N <- pdat[malignant==1, tumor]
celltype_N <- pdat[malignant==1, cell_type]

#-----------------------------------------------------------------------------#
# Seurat Transform and UMAP with 3000 top variable genes
#-----------------------------------------------------------------------------#

# Seurat normalization, no scaling, selection of highly variable genes (3000)
sM <- CreateSeuratObject(counts = malignant)
sM <- SCTransform(object = sM, return.only.var.genes = TRUE, verbose = FALSE)
sN <- CreateSeuratObject(counts = nonmalignant)
sN <- SCTransform(object = sN, return.only.var.genes = TRUE, verbose = FALSE)

# Store results from Seurat object
adj_M <- sM[["SCT"]]@scale.data
adj_N <- sN[["SCT"]]@scale.data

# Prepare data for UMAP
adj_M.t <- t(adj_M)
adj_N.t <- t(adj_N)

dat.M <- cbind(adj_M.t, tumor_M)
dat.N <- cbind(adj_N.t, tumor_N, celltype_N)

# select only tumors with >50 malignant cells
m.tab = table(tumor_M)>50
which.tumors.m = as.numeric(names(m.tab[which(m.tab==T)])) #8 tumors selected
dat.M = dat.M[tumor_M %in% which.tumors.m,]

# select only tumors with >100 non-malignant cells
n.tab = table(tumor_N)>100
which.tumors.n = as.numeric(names(n.tab[which(n.tab==T)])) #12 tumors selected
dat.N = dat.N[tumor_N %in% which.tumors.n,]
dat.N = subset(dat.N, dat.N[,3002]>0)

# UMAP
M_umap <- umap(dat.M[,1:3000])
N_umap <- umap(dat.N[,1:3000])
write.table(M_umap$layout, "./data/M_umap.txt", sep="\t")
write.table(N_umap$layout, "./data/N_umap.txt", sep="\t")


# set up data for plotting
M.plot.dat <- as.data.frame(cbind(M_umap$layout, dat.M[,3001]))
colnames(M.plot.dat) <- c("UMAP1", "UMAP2", "Tumor")
M.plot.dat$Tumor <- factor(M.plot.dat$Tumor)

N.plot.dat <- as.data.frame(cbind(N_umap$layout, dat.N[,3001], dat.N[,3002]))
colnames(N.plot.dat) <- c("UMAP1", "UMAP2", "Tumor", "CellType")
N.plot.dat$Tumor <- factor(N.plot.dat$Tumor)
N.plot.dat$CellType <- factor(N.plot.dat$CellType)

# plot feature set up
m.colors = c("darkgreen", "yellow", "black", "red",
             "limegreen", "blue", "pink", "orange")

n.colors = c("darkgreen", "darkred", "lightblue", "lightgrey",
              "yellow", "seashell3", "red","limegreen",
              "blue", "pink", "orange", "purple")

blanktheme = theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"),
                   panel.background = element_blank(),
                   legend.key = element_blank())

# UMAP plots
ggplot(M.plot.dat, aes(UMAP1, UMAP2,col=Tumor))+geom_point()+
  ggtitle("UMAP Visualization of Malignant Cells\n Top 3000 Variable Genes")+
  scale_color_manual(values=m.colors)+blanktheme
ggsave("malignant3000.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

ggplot(N.plot.dat, aes(UMAP1, UMAP2,col=Tumor))+geom_point()+
  ggtitle("UMAP Visualization of Non-Malignant Cells\n Top 3000 Variable Genes")+
  scale_color_manual(values=n.colors)+blanktheme
ggsave("nonmalignant3000_tumor.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

# color by cell-type
ggplot(N.plot.dat, aes(UMAP1, UMAP2,col=CellType))+geom_point()+blanktheme+
  ggtitle("UMAP Visualization of Non-Malignant Cells\n Top 3000 Variable Genes") + 
  scale_color_manual(values=c("red", "limegreen", "blue", "purple", "grey", "orange"), 
                      labels=c("T-cells", "B-cells", "Macrophages", "Endo.", "CAFs", "NK"))
ggsave("nonmalignant3000_celltype.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")


#-----------------------------------------------------------------------------#
# Seurat Transform and UMAP with 300 top variable genes
#-----------------------------------------------------------------------------#

# Seurat normalization, no scaling, selection of top 300 variable genes
sM300 <- CreateSeuratObject(counts = malignant)
sM300 <- SCTransform(object = sM300, verbose = FALSE,
                     return.only.var.genes = TRUE, variable.features.n = 300 )

sN300 <- CreateSeuratObject(counts = nonmalignant)
sN300 <- SCTransform(object = sN300, verbose = FALSE,
                     return.only.var.genes = TRUE, variable.features.n = 300)

# Store results from Seurat object
adj_M300 <- sM300[["SCT"]]@scale.data
adj_N300 <- sN300[["SCT"]]@scale.data

# Prepare data for UMAP
adj_M300.t <- t(adj_M300)
adj_N300.t <- t(adj_N300)

dat.M300 <- cbind(adj_M300.t, tumor_M)
dat.N300 <- cbind(adj_N300.t, tumor_N, celltype_N)

# select only tumors with >50 malignant cells
dat.M300 = dat.M300[tumor_M %in% which.tumors.m,]

# select only tumors with >100 non-malignant cells
dat.N300 = dat.N300[tumor_N %in% which.tumors.n,]
dat.N300 = subset(dat.N300, dat.N300[,302]>0)

# UMAP
M300_umap <- umap(dat.M300[,1:300])
N300_umap <- umap(dat.N300[,1:300])
write.table(M300_umap$layout, "./data/M_umap.txt", sep="\t")
write.table(N300_umap$layout, "./data/N_umap.txt", sep="\t")


# set up data for plotting
M300.plot.dat <- as.data.frame(cbind(M300_umap$layout, dat.M300[,301]))
colnames(M300.plot.dat) <- c("UMAP1", "UMAP2", "Tumor")
M300.plot.dat$Tumor <- factor(M300.plot.dat$Tumor)

N300.plot.dat <- as.data.frame(cbind(N300_umap$layout, dat.N300[,301], dat.N300[,302]))
colnames(N300.plot.dat) <- c("UMAP1", "UMAP2", "Tumor", "CellType")
N300.plot.dat$Tumor <- factor(N300.plot.dat$Tumor)
N300.plot.dat$CellType <- factor(N300.plot.dat$CellType)


# UMAP plots
ggplot(M300.plot.dat, aes(UMAP1, UMAP2,col=Tumor))+geom_point()+
  ggtitle("UMAP Visualization of Malignant Cells\n Top 300 Variable Genes")+
  scale_color_manual(values=m.colors)+blanktheme
ggsave("malignant300.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

ggplot(N300.plot.dat, aes(UMAP1, UMAP2,col=Tumor))+geom_point()+
  ggtitle("UMAP Visualization of Non-Malignant Cells\n Top 300 Variable Genes")+
  scale_color_manual(values=n.colors)+blanktheme
ggsave("nonmalignant300_tumor.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

# color by cell-type
ggplot(N300.plot.dat, aes(UMAP1, UMAP2,col=CellType))+geom_point()+blanktheme+
  ggtitle("UMAP Visualization of Non-Malignant Cells\n Top 300 Variable Genes") + 
  scale_color_manual(values=c("red", "limegreen", "blue", "purple", "grey", "orange"), 
                     labels=c("T-cells", "B-cells", "Macrophages", "Endo.", "CAFs", "NK"))
ggsave("nonmalignant300_celltype.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")





