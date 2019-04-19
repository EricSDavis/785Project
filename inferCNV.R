# inferCNV Code
# Author: Yue Pan
# Last edited: 04/19/19
########################################
# Notes:
# This will be a directly application of inferCNV in current package
########################################

library("devtools")
# install.packages("rjags", repos=c("http://cran.rstudio.com"))
library(rjags)
library(data.table)
# install.packages("rjags")
# devtools::install_github("broadinstitute/infercnv", ref="RELEASE_0_99_4")
library(infercnv)


########################################
## Start real code

# prepare for raw count matrix
data <- fread('GSE72056_melanoma_single_cell_revised_v2.txt',
              sep = "\t",header = TRUE,stringsAsFactors=FALSE)
data <- as.data.frame(data)
for (i in 2:ncol(data)){
  if (data[2,i] == 1) data[2,i] = 'non-malignant'
  if (data[2,i] == 2) data[2,i] = 'malignant'
  if (data[2,i] == 0) data[2,i] = 'unresolved'
}

# count <- data[-c(1,2,3),]
count <- data[-c(1,2,3),c(1,which(data[1,]==78))] # change here for subject
rownames <- count$Cell
count <- count[,-1]
rownames(count) <- make.names(rownames, unique=TRUE)


cell_name <- colnames(count)
# cell_type <- as.factor(data[2,])[2:ncol(data)]
cell_type <- data[2,][which(data[1,]==78)] # change here for subject


# prepare for annotation file
sample_annotation <- matrix(c(cell_name,cell_type), ncol = 2)


# final raw counts matrix
count <- as.data.frame(lapply(count, as.numeric))
colnames(count) <- cell_name
rownames(count) <- make.names(rownames, unique=TRUE)
count <- (2^count - 1)*10



# inferCNV object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(count),
                                    annotations_file="cellAnnotations_78.txt",
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names=c("non-malignant"),
                                    #ref_group_names = NULL,
                                    chr_exclude=c('chrY', 'chrM'))



# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             # max_centered_threshold = NA,
                             out_dir="output_dir_80_test",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             # HMM=T
                            
)


# final plot
plot_cnv(infercnv_obj,
         out_dir="output_dir_78plot",
         obs_title="malignant/unknown",
         ref_title="non-malignant",
         cluster_by_groups=TRUE,
         #k_obs_groups = 3,
         contig_cex=2.5,
         x.center=mean(infercnv_obj@expr.data),
         x.range="auto", #NA,
         hclust_method='ward.D',
         color_safe_pal=FALSE,
         output_filename="infercnv_plot",
         output_format="png", #pdf, png, NA
         png_res=300,
         dynamic_resize=0,
         ref_contig = NULL,
         write_expr_matrix=FALSE)


