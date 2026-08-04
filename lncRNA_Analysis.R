
########### Identifying long Non-Coding RNAs in Parkinsons Disease ############

####################### Dataset ID: GSE157783 #################################

########## Single nucleus RNA-Seq from PD & control Substantia Nigra ###########

########## Setup libraries & packages ########################

library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(ggplot2)




############ STEP-1: LOAD THE DATA #####################
?readMM

########## As counts is in tsv format, so use read.table function ##############

counts <- read.table("Data/GSE157783_UMI/midbrain_UMI.tsv", sep = "\t", header = TRUE)

########### Convert to matrix ########################

counts_matrix <- as(as.matrix(counts), "dgCMatrix")

head(counts_matrix) ##### Check the structure of matrix

dim(counts_matrix) ###### Check dimensions: rows as genes & columns as barcodes

###### 26737 genes & 41435 cells ######

########## Reading barcodes info ###############

barcodes <- read.table("Data/GSE157783_cell/midbrain_cell.tsv", sep = "\t", header = TRUE)

head(barcodes)

dim(barcodes)

########## Reading gene information #############

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

######## Match these IDs against rownames of counts matrix ################

mt_features <- intersect(mt_ensembl$Gene.stable.ID, rownames(seurat))

############ Calculate percent.mt #####################

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mt_features)

summary(seurat$percent.mt) ###### Returns percent.mt 

################ Plot the features, counts and percent.mt ##################

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

#### Violin plot & percent.mt shows that no true mitochondrial genes exist in dataset #########

####### Check if the no of genes & no of transcripts are correlated across cells ########


corr.plot <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

corr.plot

##### Corr.plot shows that as number of genes increases, the transcript out also increase ######


################### STEP-4: FILTERING ##################################

######### Filtering cells with > 1500 UMIs per cell, > 1000 genes per cell & percent.mt < 10% #######

###### These filtering thresholds are obtained from original paper ############

#### Any cell that shows < 1500 UMIs & 1000 genes & > 10% mtDNA will be filtered #######

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nCount_RNA > 1500 & percent.mt < 10)


################# STEP-5: NORMALIZATION ##############################

#### Make gene expression levels comparable between cells ################

seurat <- NormalizeData(seurat)


#################### STEP-6: FIND TOP VARIABLE FEATURES ######################

####### Identify genes with most variable expression across cells ############

seurat <- FindVariableFeatures(seurat, nfeatures = 2000)

seurat


##################### STEP-7: SCALING DATA ############################

seurat <- ScaleData(seurat)

seurat

################### STEP-8: PRINCIPAL COMPONENT ANALYSIS ####################

########## Dimensionality reduction #################

seurat <- RunPCA(seurat, npcs = 50)


ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))

########## Elbow plot shows that first 15-20 PCs explain the most variance #######


##################### STEP-8: NON-LINEAR DIMENSIONALITY REDUCTION ##############

seurat <- RunTSNE(seurat, dims = 1:20)

seurat <- RunUMAP(seurat, dims = 1:20)

######### Visualize UMAP & tSNE ##################

plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2

#### Selecting UMAP as it preserves local and global structure of dataset ####


################ STEP-9: CLUSTER THE CELLS ########################

seurat <- FindNeighbors(seurat, dims = 1:20)

seurat <- FindClusters(seurat, resolution = 1)

########## Change the resolution to check the number of clusters ##########

seurat <- FindClusters(seurat, resolution = 0.2)

seurat <- FindClusters(seurat, resolution = 0.1) ####### Yields 13 clusters 


seurat <- FindClusters(seurat, resolution = 0.1) ####### Yields 15 clusters 


DimPlot(seurat, reduction = "umap", label = TRUE)

######### Final resolution: 0.1 ##################

################## STEP-10: ADD GENE SYMBOLS TO THE DATA ##################

######### The current dataset has ensemble ids #############

######## Add gene symbols for better understanding ###########

gene_map <- ensembldb::select(EnsDb.Hsapiens.v86, keys = rownames(seurat),
                              keytype = "GENEID", columns = "SYMBOL")

gene_map <- gene_map[!duplicated(gene_map$GENEID), ]  #### keep mapping per Ensembl ID

seurat[["RNA"]]@meta.data$symbol <- gene_map$SYMBOL[match(rownames(seurat), gene_map$GENEID)]


############## Some genes have NA symbol #####################

sum(is.na(seurat[["RNA"]]@meta.data$symbol)) ##### Check the number of NA ids

########## Get the list of NA ids #########################

na_ids <- rownames(seurat)[is.na(seurat[["RNA"]]@meta.data$symbol)]

length(na_ids) ##### 192 Genes have NA symbol 

na_ids

########## Assign the ensemble id back to the NA symbol ##########

symbol_col <- seurat[["RNA"]]@meta.data$symbol

na_idx <- is.na(symbol_col)

symbol_col[na_idx] <- rownames(seurat)[na_idx]

seurat[["RNA"]]@meta.data$symbol <- symbol_col

sum(is.na(seurat[["RNA"]]@meta.data$symbol))  # should now be 0

seurat[["RNA"]]@meta.data$symbol


################### STEP-11: REMOVE DOUBLETS #######################
saveRDS(seurat, file = "seurat.rds")


symbols <- ifelse(is.na(seurat[["RNA"]]@meta.data$symbol), rownames(seurat), seurat[["RNA"]]@meta.data$symbol)

valid.symbols <- make.unique(symbols)

######### Assin symbols to rownames ##################

rownames(seurat[["RNA"]]) <- valid.symbols

rownames(seurat@assays$RNA) <- valid.symbols

################ STEP-11: REMOVE DOUBLETS ################################
########### Using scDblFinder ######################

###### Convert seurat object to SingleCellExperiment ###################

sce <- as.SingleCellExperiment(seurat)

sce$seurat_clusters

############# Calculate doublets #########################

sce <- scDblFinder(sce, clusters = "seurat_clusters")

sce$scDblFinder.score ##### Check score; higher score, greater chance of doublet

sce$scDblFinder.class ### Check class; singlet or doublet 

########## Assign score and class back to seurat object ###############

seurat$scDblFinder.score <- sce$scDblFinder.score

seurat$scDblFinder.class <- sce$scDblFinder.class


########### Visualize ################

DimPlot(seurat, reduction = "umap", group.by = "scDblFinder.class")

######### Subset the singlets only ##################

seurat <- subset(seurat, subset = scDblFinder.class == "singlet")


################## STEP-12: RE-CLUSTER THE CELLS ####################

seurat <- FindVariableFeatures(seurat, nfeatures = 2000) ### Find Variable features

seurat <- ScaleData(seurat) ### Scale the data

seurat <- RunPCA(seurat, npcs = 50) ### PCA

ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca"))) ### Plot PCs

seurat <- RunUMAP(seurat, dims = 1:20) ### Non-Linear dimensionality reduction; UMAP

seurat <- FindNeighbors(seurat, dims = 1:20) ### Find neighbors

seurat <- FindClusters(seurat, resolution = 0.1) ### Find clusters

DimPlot(seurat, reduction = "umap", label = TRUE) ### Plot the clusters

### 12 Clusters identified

################ STEP-13: FIND MARKERS OF EACH CLUSTER ##################

#### Clusters 3 & 11 are almost merged together #############

#### First check if there is an actual biological difference between clusters 3 & 11 ######

?FindMarkers
diff_3_11 <- FindMarkers(seurat,
                         ident.1 = 3,
                         ident.2 = 11,
                         min.pct = 0.25,
                         logfc.threshold = 0.25)
head(diff_3_11, 10)

##### The p-values are significant, showing that they are different clusters #####

########### Find Markers of All Clusters ##################

cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(cl_markers, "Results/markers.csv")


#################### STEP-14: ANNOTATE CELL CLUSTERS ########################

###### Visualize some cell markers to see theri expression in clusters ##########

######## Check the expression of each marker with known literature ########

####### Original Study: https://pmc.ncbi.nlm.nih.gov/articles/PMC9050543/ ####

##### Another: https://doi.org/10.1101/2025.08.26.672469 #######

######### Astrocytes: AQP4 & GFAP ###################
#### AQP4: Cluster 2 & 8
#### GFAP: Cluster 2,4,6,7,8 & 10

######### Oligodendrocytes: MBP, MOBP, MOG & OPALIN ###########
#### MOG & MOBP: Clusters 0,1 & 11
#### OPALIN: Cluster 0

########### Oligodendrocyte Precursor Cells (OPCs): VCAN, PDGFRA, PCDH15 ###########
##### VCAN: Clusters 2,4 & 8
##### PDGFRA: Clusters 4
##### PCDH15: Clusters 4,5 & 9

############ Ependymal Cells: FOXJ1 ###############
#### FOXJ1: Cluster 8

############# Microglia: CD74 & APBB1IP #################
#### CD74: Clusters 3,6 & 11
#### APBB1IP: Clusters 3 & 11

################ Endothelial Cells: CLDN5 ###################
##### CLDN5: Clusters 6 & 7

################## Excitatory Neurons: SLC17A6 ##############
##### SLC17A6: Cluster 5

################### Inhibitory Neurons: GAD1 & GAD2 ################
##### GAD1: Clusters 4,5 & 9
##### GAD2: Clusters 5 & 9

#################### GABAergic Neurons: GAD2 & GRIK1 ###############
##### GAD2: Clusters 5 & 9
##### GRIK1: Clusters 4,5 & 9

################## Dopaminergic Neurons: TH #################
##### TH: Not detectable in any clusters. 


##########  Define a vector of canonical markers ############

markers <- c(
  "AQP4", "GFAP",          # Astrocytes
  "MOBP", "MBP", "MOG", "OPALIN",    # Oligodendrocytes
  "PDGFRA", "VCAN", "PCDH15",  # OPCs
  "FOXJ1",                  # Ependymal cells
  "CD74", "APBB1IP",       # Microglia
  "CLDN5",                 # Endothelial
  "SLC17A6",  # Excitatory Neurons
  "GAD1", "GAD2", # Inhibitory Neurons
  "GRIK1", # GABAergic Neurons
  "TH", "SLC6A3", "ALDH1A1", "NR4A2", "SLC18A2", #Dopaminergic Neurons
  "PDGFRB", "RGS5", "NOTCH3", "ABCC9", # Pericytes
  "CADPS2" # Mentioned in paper
)

#######  Generate DotPlot to check the no of cells expressing markers ######

DotPlot(seurat, features = markers) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######## Violin Plot ##########

VlnPlot(seurat, features = markers, pt.size = 0)

table(Idents(seurat))

########### Find markers of unambiguous clusters ##################

cluster5.markers <- FindMarkers(seurat,
                                ident.1 = 5,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)

head(cluster5.markers, 50)

cluster9.markers <- FindMarkers(seurat,
                                ident.1 = 9,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)

head(cluster9.markers, 50)

cluster10.markers <- FindMarkers(seurat,
                                 ident.1 = 10,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)

head(cluster10.markers, 50)
View(cluster10.markers)

cluster11.markers <- FindMarkers(seurat,
                                 ident.1 = 11,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)

head(cluster11.markers, 50)

cluster8.marker <- FindMarkers(seurat,
                               ident.1 = 8,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

head(cluster8.marker, 50)


############## ANNOTATE CLUSTERS BASED ON MARKER EXPRESSION #############

markers <- c(
  "AQP4", "GFAP",          # Astrocytes
  "MOBP", "MBP", "MOG", "OPALIN",    # Oligodendrocytes
  "PDGFRA", "VCAN", "PCDH15",  # OPCs
  "FOXJ1",                  # Ependymal cells
  "CD74", "APBB1IP",       # Microglia
  "CLDN5",                 # Endothelial
  "SLC17A6",  # Excitatory Neurons
  "GAD1", "GAD2", # Inhibitory Neurons
  "GRIK1", # GABAergic Neurons
  "TH", "SLC6A3", "ALDH1A1", "NR4A2", "SLC18A2", #Dopaminergic Neurons
  "PDGFRB", "RGS5", "NOTCH3", "ABCC9", # Pericytes
  "CADPS2",  # Mentioned in paper
  "ITK", "IL7R", "CD96" # T-cells Cluster 11 markers
)

#######  Generate DotPlot to check the no of cells expressing markers ######

DotPlot(seurat, features = markers) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


########### Final Annotation based on cell markers ############

# 0 : OPALIN+ Oligodendrocytes
# 1 : Oligodendrocytes
# 2 : Astrocytes 
# 3 : Microglia
# 4 : Oligodendrocytes precursor cells (OPCs)
# 5 : Excitatory Neurons
# 6 : Endothelial Cells
# 7 : Pericytes
# 8 : Ependymal Cells
# 9 : Inhibitory Neurons
# 10 : CADPS2
# 11 : T-Lymphocyte. 


new_ident <- setNames(c("OPALIN+ Oligodendrocytes",
                        "Oligodendrocytes",
                        "Astrocytes",
                        "Microglia",
                        "Oligodendrocytes precursor cells (OPCs)",
                        "Excitatory Neurons",
                        "Endothelial Cells",
                        "Pericytes",
                        "Ependymal Cells",
                        "Inhibitory Neurons",
                        "CADPS2+",
                        "T cells"),
                      levels(seurat))

seurat <- RenameIdents(seurat, new_ident)

DimPlot(seurat, reduction = "umap", label = TRUE)

