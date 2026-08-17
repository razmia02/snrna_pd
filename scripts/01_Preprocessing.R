
############### R SCRIPT FOR PRE PROCESSING THE DATA ###########################

#### Filtering, QC, Doublet removal, Normalization, Scaling, UMAP ####

########### Identifying long Non-Coding RNAs in Parkinsons Disease ############

####################### Dataset ID: GSE157783 #################################

########## Single nucleus RNA-Seq from PD & control Substantia Nigra ###########

########## Setup libraries & packages ########################

library(remotes)
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(scDblFinder)
library(SingleCellExperiment)
library(harmony)
library(ggplot2)
library(EnhancedVolcano)
library(DESeq2)
library(muscat)



################################################################################

############ STEP-1: LOAD THE DATA #####################
?readMM

########## As counts is in tsv format, so use read.table function ##############

counts <- read.table("Data/GSE157783_UMI/midbrain_UMI.tsv", sep = "\t", header = TRUE)

rownames(counts)

dim(counts)

########### Convert to matrix ########################

counts_matrix <- as(as.matrix(counts), "dgCMatrix")

head(counts_matrix) ##### Check the structure of matrix

rownames(counts_matrix)

colnames(counts_matrix)

dim(counts_matrix) ###### Check dimensions: rows as genes & columns as barcodes

###### 26737 genes & 41435 cells ######

########## Reading barcodes info ###############

barcodes <- read.table("Data/GSE157783_cell/midbrain_cell.tsv", sep = "\t", header = TRUE)

head(barcodes)

dim(barcodes)

barcodes$barcode <- gsub("-", ".", barcodes$barcode)

rownames(barcodes)

colnames(barcodes)

########## Reading gene information #############

features <- read.table("Data/GSE157783_genes/midbrain_genes.tsv", sep="\t", header=TRUE)

head(features)

dim(features)

rownames(features)

colnames(features)

########### Notice that gene column contains ensemble ids ###############


################################################################################

################ STEP-2: ADD GENE SYMBOLS AS ROWNAMES ##############

###### The current rownames/gene ids are ensemble ids #########

features$gene

###### Convert the ensemble ids to gene symbols ########

######## Add gene symbols for better understanding ###########

gene_map <- ensembldb::select(EnsDb.Hsapiens.v86, 
                              keys = as.character(features$gene),
                              keytype = "GENEID", columns = "SYMBOL")

gene_map <- gene_map[!duplicated(gene_map$GENEID), ]  #### keep mapping per Ensembl ID


symbols <- gene_map$SYMBOL[match(features$gene, gene_map$GENEID)]

symbols[is.na(symbols)] <- features$gene[is.na(symbols)] 

rownames(counts_matrix) <- make.unique(symbols) 

rownames(counts_matrix)

colnames(counts_matrix)

head(counts_matrix)


################################################################################

######################## STEP-3: CREATE SEURAT OBJECT #########################

seurat <- CreateSeuratObject(counts_matrix, project="PD")

seurat

dim(seurat)

rownames(seurat)

Layers(seurat)


#################################################################################

################## STEP-4: ADD METADATA INFORMATION TO SEURAT OBJECT ###########

######### Compare column names of seurat & barcodes ############

head(colnames(seurat)) ### Check column names of seurat

head(barcodes$barcode) ### Should match the seurat column names

sum(barcodes$barcode %in% colnames(seurat)) ### 41435

ncol(seurat) ### 41435

rownames(barcodes) <- barcodes$barcode #### Change the rownames 

#### Add metadata directly from barcodes ########

seurat <- AddMetaData(seurat, metadata = barcodes)

###### Add condition column to metadata ##########

table(seurat$patient, seurat$condition)

seurat$condition <- ifelse(grepl("PD", seurat$patient), "PD", "Control")

table(seurat$patient, seurat$condition)


################################################################################

###################### STEP-5: QUALITY CONTROL ##############################

########## Remove cells with too few or too many genes ############

######### Filter out high mitochondrial percentages ###############

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

summary(seurat$percent.mt) ##### Returns 0 mitochondrial genes. 


################ Plot the features, counts and percent.mt ##################

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

#### Violin plot & percent.mt shows that no true mitochondrial genes exist in dataset #########

####### Check if the no of genes & no of transcripts are correlated across cells ########


corr.plot <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

corr.plot


##### Corr.plot shows that as number of genes increases, the transcript no. also increase ######


################################################################################

################### STEP-6: FILTERING ##################################

######### Filtering cells with > 1500 UMIs per cell, > 1000 genes per cell & percent.mt < 10% #######

###### These filtering thresholds are obtained from original paper ############

#### Any cell that shows < 1500 UMIs & 1000 genes & > 10% mtDNA will be filtered #######

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nCount_RNA > 1500 & percent.mt < 10)

seurat@meta.data$patient


################################################################################

################ STEP-7: REMOVE DOUBLETS ################################

########### Using scDblFinder ######################

###### Convert seurat object to SingleCellExperiment ###################

sce <- as.SingleCellExperiment(seurat)

saveRDS(seurat, "seurat.rds")

############# Calculate doublets #########################
?scDblFinder

sce <- scDblFinder(sce, samples = "orig.ident")

sce$scDblFinder.score ##### Check score; higher score, greater chance of doublet

sce$scDblFinder.class ### Check class; singlet or doublet 

########## Assign score and class back to seurat object ###############

seurat$scDblFinder.score <- sce$scDblFinder.score

seurat$scDblFinder.class <- sce$scDblFinder.class


########### Visualize ################

#### See total Singlets vs Doublets count ####

table(seurat$scDblFinder.class)

##### Breakdown of doublets per patient sample #####

table(seurat$patient, seurat$scDblFinder.class)

###### Barplot of Doublets per Patient ######

ggplot(as.data.frame(seurat@meta.data), aes(x = patient, fill = scDblFinder.class)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(y = "Proportion", title = "Doublet Proportion by Patient")

######### Subset the singlets only ##################

seurat <- subset(seurat, subset = scDblFinder.class == "singlet")


################################################################################

################# STEP-8: NORMALIZATION ##############################

#### Make gene expression levels comparable between cells ################

seurat <- NormalizeData(seurat)


################################################################################

#################### STEP-9: FIND TOP VARIABLE FEATURES ######################

####### Identify genes with most variable expression across cells ############

seurat <- FindVariableFeatures(seurat, nfeatures = 2000)

seurat


################################################################################

##################### STEP-10: SCALING DATA ############################

seurat <- ScaleData(seurat)

seurat


################################################################################

################### STEP-11: PRINCIPAL COMPONENT ANALYSIS ####################

########## Dimensionality reduction #################

seurat <- RunPCA(seurat, npcs = 50)


ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))


########## Elbow plot shows that first 20 PCs explain the most variance #######


################################################################################

################# STEP-12: BATCH CORRECTION ##########################

######## Using Harmony ####################

########### Run UMAP before and after Harmony to check batch correction ####

seurat <- RunUMAP(seurat, reduction = "pca", 
                  dims = 1:20, 
                  reduction.name = "umap.unintegrated")

plot1 <- DimPlot(seurat, reduction = "umap.unintegrated", 
                 group.by = "patient") + ggtitle("Before Harmony")

plot1

############ Run Harmony for Batch Correction ###############
?RunHarmony
seurat <- RunHarmony(seurat, 
                     group.by.vars = "patient", 
                     dims.use = 1:20, max_iter = 50)


################################################################################

##################### STEP-13: NON-LINEAR DIMENSIONALITY REDUCTION ##############

######## UMAP ###########

seurat <- RunUMAP(seurat, reduction = "harmony", 
                  dims = 1:20, reduction.name = "umap")


####### Plot to see the difference #################

plot2 <- DimPlot(seurat, reduction = "umap", 
                 group.by = "patient") + ggtitle("After Harmony")

plot2


plot1 + plot2 #### See if the batch correction was applied correctly


################################################################################

################ STEP-14: CLUSTER THE CELLS ########################

seurat <- FindNeighbors(seurat, dims = 1:20)

seurat <- FindClusters(seurat, resolution = 1) ### Yields 24 clusters 


########## Change the resolution to check the number of clusters ##########

seurat <- FindClusters(seurat, resolution = 0.2) #### Yields 15 clusters 

seurat <- FindClusters(seurat, resolution = 0.1) ####### Yields 13 clusters 


DimPlot(seurat, reduction = "umap", label = TRUE)

######### Final resolution: 0.1 ##################

saveRDS(seurat, "seurat.rds")

################################################################################
