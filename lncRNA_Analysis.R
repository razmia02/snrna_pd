
########### Identifying long Non-Coding RNAs in Parkinsons Disease ############

####################### Dataset ID: GSE157783 #################################

########## Single nucleus RNA-Seq from PD & control Substantia Nigra ###########

########## Setup libraries & packages ########################

library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)


############ STEP-1: LOAD THE DATA #####################
?readMM

########## As counts is in tsv format, so use read.table function ##############

counts <- read.table("Data/GSE157783_UMI/midbrain_UMI.tsv", sep = "\t", header = TRUE)

########### Convert to matrix ########################

counts_matrix <- as(as.matrix(counts), "dgCMatrix")

head(counts_matrix) ##### Check the structure of matrix

dim(counts_matrix) ###### Check dimensions: rows as genes & columns as barcodes

###### 26737 genes & 41435 cells ######

########## Read in barcodes info ###############

barcodes <- read.table("Data/GSE157783_cell/midbrain_cell.tsv", sep = "\t", header = TRUE)

head(barcodes)

dim(barcodes)

########## Readin gene information #############

features <- read.table("Data/GSE157783_genes/midbrain_genes.tsv", sep="\t", header=TRUE)

head(features)

dim(features)

########### Rename matrix rows as per gene names ##############

rownames(counts_matrix) <- make.unique(as.character(features[,1]))

head(rownames(counts_matrix))

############## Rename matrix columns to barcodes #################

colnames(counts_matrix) <- barcodes[,1]

head(colnames(counts_matrix))

######################## STEP-2: CREATE SEURAT OBJECT #########################

seurat <- CreateSeuratObject(counts_matrix, project="PD")

seurat

###################### STEP-3: QUALITY CONTROL ##############################

########## Remove cells with too few or too many genes ############

######### Filter out high mitochondrial percentages ###############

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

summary(seurat$percent.mt) ##### Returns 0 mitochondrial genes. 

####### Since features are in ensemble format, calculate percent.mt manually ########


######### Obtain complete list of mitochondrial genes from biomart ensemble ################

mt_ensembl <- read.table("Data/mito_id.txt", sep="\t", header=TRUE)

mt_ensembl

mt_ensembl <- as.character(mt_ensembl[, 1])

mt_ensembl

######## Match these IDs against rownames of counts matrix ################

mt_features <- intersect(mt_ensembl, rownames(seurat))

############ Calculate percent.mt #####################

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mt_features)

summary(seurat$percent.mt) ###### Returns percent.mt 

################ Plot the features, counts and percent.mt ##################

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
