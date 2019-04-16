library(umap)
library(data.table)
#-----------------------------------------------------------------------------#

# Data set up
setwd("/Users/tmlagler/OneDrive/Third Year 2018-2019/785 - Gene Assoc/Project/")
data = fread("GSE72056_melanoma_single_cell_revised_v2.txt", header=T)

data[2,1] = "malignant"; data[3,1] = "non_malignant"

data.t = transpose(data[,-1])
colnames(data.t) = unlist(data[,1])

data.m = data.t[malignant==2]
data.nm = data.t[malignant==1]

#-----------------------------------------------------------------------------#
# select only tumors with >50 malignant cells
m.tab = table(data.m$tumor)>50
which.tumors.m = as.numeric(names(m.tab[which(m.tab==T)])) #8 tumors selected
data.m = data.m[tumor %in% which.tumors.m,]

data.m.tumor = as.factor(unlist(data.m[,1]))
data.m.cells = data.m[,-(1:3)]

# select only tumors with >50 malignant cells
nm.tab = table(data.nm$tumor)>100
which.tumors.nm = as.numeric(names(nm.tab[which(nm.tab==T)])) #12 tumors selected
data.nm = data.nm[tumor %in% which.tumors.nm,]

data.nm.celltype = as.factor(unlist(data.nm[,3]))
data.nm.tumor = as.factor(unlist(data.nm[,1]))
data.nm.cells = data.nm[,-(1:3)]

#-----------------------------------------------------------------------------#
# Run UMAP
m.umap = umap(data.m.cells)
nm.umap = umap(data.nm.cells)

write.table(m.umap$layout, "m_umap.txt", sep="\t")
write.table(nm.umap$layout, "nm_umap.txt", sep="\t")
# unadj_m.umap <- read.csv("m_umap.txt", sep="\t")
# unadj_nm.umap <- read.csv("nm_umap.txt", sep="\t")

#-----------------------------------------------------------------------------#
# UMAP visualizaiton using plot features defined in UMAP.R
unadj_M_plot <- cbind(unadj_m.umap, data.m.tumor)
colnames(unadj_M_plot) <- c("UMAP1", "UMAP2", "Tumor")

unadj_N_plot <- cbind(unadj_nm.umap, data.nm.tumor, data.nm.celltype)
colnames(unadj_N_plot) <- c("UMAP1", "UMAP2", "Tumor", "CellType")


ggplot(unadj_M_plot, aes(UMAP1, UMAP2,col=Tumor))+geom_point()+
  ggtitle("UMAP Visualization of Malignant Cells\n Unadjusted")+
  scale_color_manual(values=m.colors)+blanktheme
ggsave("malignant_unadj.jpg",
       width = 6, height = 4, units = "in")

ggplot(unadj_N_plot, aes(UMAP1, UMAP2,col=Tumor))+geom_point()+
  ggtitle("UMAP Visualization of Non-Malignant Cells\n Unadjusted")+
  scale_color_manual(values=n.colors)+blanktheme
ggsave("nonmalignant_unadj_tumor.jpg",
       width = 6, height = 4, units = "in")

# color by cell-type
ggplot(unadj_N_plot, aes(UMAP1, UMAP2,col=CellType))+geom_point()+blanktheme+
  ggtitle("UMAP Visualization of Non-Malignant Cells\n Unadjusted") + 
  scale_color_manual(values=c("yellow", "red", "limegreen", "blue", "purple", "grey", "orange"), 
                     labels=c("Unknown", "T-cells", "B-cells", "Macrophages", "Endo.", "CAFs", "NK"))
ggsave("nonmalignant_unadj_celltype.jpg",
       width = 6, height = 4, units = "in")












# 59 71 78 79
# 80 81 88 89
m.colors = c("darkgreen", "yellow", "black", "red",
             "green", "blue", "pink", "orange")

# 53 58 60 72
# 74 75 79 80
# 84 88 89 94
nm.colors = c("darkgreen", "darkred", "lightblue", "lightgrey",
              "yellow", "seashell3", "red","green",
              "blue", "pink", "orange", "purple")

plot.umap(m.umap, labels = data.m.tumor,
          color = m.colors, 
          main="A UMAP visualization of the Melanoma dataset: Malignant Cells")

plot.umap(nm.umap, data.nm.celltype,
          color = nm.colors,
          main="A UMAP visualization of the Melanoma dataset: Non-malignant Cells")

plot.umap = function(x, labels, color,
                     main="A UMAP visualization of the Melanoma dataset",
                     pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
                     cex.main=1, cex.legend=1) {
  
  layout = x
  if (class(x)=="umap") {
    layout = x$layout
  }
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=color[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    #legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend("bottomright", legend=legend.text,
         col=color[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}
#https://cran.r-project.org/web/packages/umap/vignettes/umap.html
