library(SIMLR)
library(parallel)
library(Matrix)
library(igraph)

# read in Seurat created data
# top 300 variable genes
adjM300_data <- readRDS("./data/adjM300_data.rds")
adjN300_data <- readRDS("./data/adjN300_data.rds")

simlrM300 <- SIMLR(X = t(adjM300_data[,1:300]), c=8, cores.ratio=0, normalize=F) #16 minutes
simlrM300_large <- SIMLR_Large_Scale(X = t(adjM300_data[,1:300]), c=8, normalize=F) #3 seconds
simlrN300_large <- SIMLR_Large_Scale(X = t(adjN300_data[,1:300]), c=6, normalize=F) #6 seconds

# To assess the performance of our method, we compute the normalized mutual
# information (NMI) between the clusters inferred by SIMLR and the ground truth
# clusters. NMI takes values between 0 and 1, with higher values reflecting better
# performance.
nmi_M300 = compare(adjM300_data[,301], simlrM300$y$cluster, method="nmi")
nmi_M300

nmi_M300_large = compare(adjM300_data[,301], simlrM300_large$y$cluster, method="nmi")
nmi_M300_large

nmi_N300_large = compare(adjN300_data[,302], simlrN300_large$y$cluster, method="nmi")
nmi_N300_large

# confusion matrix of clustering
table(simlrM300$y$cluster,adjM300_data[,301])
table(simlrM300_large$y$cluster,adjM300_data[,301])
table(simlrN300_large$y$cluster,adjN300_data[,302])

# prepare data for plotting
simlrM300_dat <- cbind.data.frame(simlrM300$ydata, simlrM300$y$cluster, adjM300_data[,301])
colnames(simlrM300_dat) <- c("SIMLR1", "SIMLR2", "Cluster", "Tumor")
simlrM300_dat$Tumor <- factor(simlrM300_dat$Tumor)
simlrM300_dat$Cluster <- factor(simlrM300_dat$Cluster)

simlrM300_large_dat <- cbind.data.frame(simlrM300_large$ydata, simlrM300_large$y$cluster, adjM300_data[,301])
colnames(simlrM300_large_dat) <- c("SIMLR1", "SIMLR2", "Cluster", "Tumor")
simlrM300_large_dat$Tumor <- factor(simlrM300_large_dat$Tumor)
simlrM300_large_dat$Cluster <- factor(simlrM300_large_dat$Cluster)

simlrN300_large_dat <- cbind.data.frame(simlrN300_large$ydata, simlrN300_large$y$cluster,
                                        adjN300_data[,301], adjN300_data[,302])
colnames(simlrN300_large_dat) <- c("SIMLR1", "SIMLR2", "Cluster", "Tumor", "CellType")
simlrN300_large_dat$Tumor <- factor(simlrN300_large_dat$Tumor)
simlrN300_large_dat$Cluster <- factor(simlrN300_large_dat$Cluster)
simlrN300_large_dat$CellType <- factor(simlrN300_large_dat$CellType)


# plot set up (same as UMAP)
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


# SIMLR plots

# Malignant colored by tumor
ggplot(simlrM300_dat, aes(SIMLR1,SIMLR2,col=Tumor))+geom_point()+blanktheme+
  scale_color_manual(values=m.colors)+
  ggtitle("SIMLR Clustering of Malignant Cells\nTop 300 Variable Genes\nColored by Tumor")
ggsave("simlrM300_tumor.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

# Malignant colored by cluster
ggplot(simlrM300_dat, aes(SIMLR1,SIMLR2,col=Cluster))+geom_point()+blanktheme+
  scale_color_manual(values=m.colors)+
  ggtitle("SIMLR Clustering of Malignant Cells\nTop 300 Variable Genes\nColored by Cluster")
ggsave("simlrM300_cluster.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

# Malignant colored by tumor (fast pca)
ggplot(simlrM300_large_dat, aes(SIMLR1,SIMLR2,col=Tumor))+geom_point()+blanktheme+
  scale_color_manual(values=m.colors)+
  ggtitle("SIMLR_Large_Scale Clustering of Malignant Cells\nTop 300 Variable Genes\nColored by Tumor")
ggsave("simlrM300_large_tumor.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

# Malignant colored by cluster (fast pca)
ggplot(simlrM300_large_dat, aes(SIMLR1,SIMLR2,col=Cluster))+geom_point()+blanktheme+
  scale_color_manual(values=m.colors)+
  ggtitle("SIMLR_Large_Scale Clustering of Malignant Cells\nTop 300 Variable Genes\nColored by Cluster")
ggsave("simlrM300_large_cluster.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

# Non-Malignant colored by cell type (fast pca)
ggplot(simlrN300_large_dat, aes(SIMLR1,SIMLR2,col=CellType))+geom_point()+blanktheme+
  scale_color_manual(values=m.colors)+
  ggtitle("SIMLR_Large_Scale Clustering of Non-Malignant Cells\nTop 300 Variable Genes\nColored by Cell Type")
ggsave("simlrN300_large_celltype.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

# Non-Malignant colored by tumor (fast pca)
ggplot(simlrN300_large_dat, aes(SIMLR1,SIMLR2,col=Tumor))+geom_point()+blanktheme+
  scale_color_manual(values=n.colors)+
  ggtitle("SIMLR_Large_Scale Clustering of Non-Malignant Cells\nTop 300 Variable Genes\nColored by Tumor")
ggsave("simlrN300_large_tumor.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")

# Non-Malignant colored by cluster (fast pca)
ggplot(simlrN300_large_dat, aes(SIMLR1,SIMLR2,col=Cluster))+geom_point()+blanktheme+
  scale_color_manual(values=m.colors)+
  ggtitle("SIMLR_Large_Scale Clustering of Non-Malignant Cells\nTop 300 Variable Genes\nColored by Cluster")
ggsave("simlrN300_large_cluster.jpg", path = "./plots/",
       width = 6, height = 4, units = "in")
