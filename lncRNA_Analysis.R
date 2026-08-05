
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
library(harmony)
library(scDblFinder)
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


################## STEP-3: ADD METADATA INFORMATION TO SEURAT OBJECT ###########

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

seurat@meta.data

################ STEP-4: ADD GENE SYMBOLS AS ROWNAMES ##############

###### The current rownames/gene ids are ensemble ids #########

rownames(seurat)

###### Convert the ensemble ids to gene symbols ########

######## Add gene symbols for better understanding ###########

gene_map <- ensembldb::select(EnsDb.Hsapiens.v86, keys = rownames(seurat),
                              keytype = "GENEID", columns = "SYMBOL")

gene_map <- gene_map[!duplicated(gene_map$GENEID), ]  #### keep mapping per Ensembl ID


seurat[["RNA"]]@meta.data$symbol <- gene_map$SYMBOL[match(rownames(seurat), gene_map$GENEID)]


############## Some genes have NA symbol #####################

sum(is.na(seurat[["RNA"]]@meta.data$symbol)) ##### Check the number of NA ids


########## Assign the ensemble id back to the NA symbol ##########

symbol_col <- seurat[["RNA"]]@meta.data$symbol

na_idx <- is.na(symbol_col)

symbol_col[na_idx] <- rownames(seurat)[na_idx]

seurat[["RNA"]]@meta.data$symbol <- symbol_col

sum(is.na(seurat[["RNA"]]@meta.data$symbol))  # should now be 0

seurat[["RNA"]]@meta.data$symbol
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


################### STEP-6: FILTERING ##################################

######### Filtering cells with > 1500 UMIs per cell, > 1000 genes per cell & percent.mt < 10% #######

###### These filtering thresholds are obtained from original paper ############

#### Any cell that shows < 1500 UMIs & 1000 genes & > 10% mtDNA will be filtered #######

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nCount_RNA > 1500 & percent.mt < 10)

seurat

################ STEP-7: REMOVE DOUBLETS ################################

########### Using scDblFinder ######################

###### Convert seurat object to SingleCellExperiment ###################

sce <- as.SingleCellExperiment(seurat)


############# Calculate doublets #########################

sce <- scDblFinder(sce)

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

################# STEP-8: NORMALIZATION ##############################

#### Make gene expression levels comparable between cells ################

seurat <- NormalizeData(seurat)


#################### STEP-9: FIND TOP VARIABLE FEATURES ######################

####### Identify genes with most variable expression across cells ############

seurat <- FindVariableFeatures(seurat, nfeatures = 2000)

seurat

##################### STEP-10: SCALING DATA ############################

seurat <- ScaleData(seurat)

seurat

################### STEP-11: PRINCIPAL COMPONENT ANALYSIS ####################

########## Dimensionality reduction #################

seurat <- RunPCA(seurat, npcs = 50)


ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))


########## Elbow plot shows that first 20 PCs explain the most variance #######

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


##################### STEP-13: NON-LINEAR DIMENSIONALITY REDUCTION ##############

######## UMAP ###########

seurat <- RunUMAP(seurat, reduction = "harmony", 
                  dims = 1:20, reduction.name = "umap")


####### Plot to see the difference #################

plot2 <- DimPlot(seurat, reduction = "umap", 
              group.by = "patient") + ggtitle("After Harmony")

plot2


plot1 + plot2 #### See if the batch correction was applied correctly

################ STEP-14: CLUSTER THE CELLS ########################

seurat <- FindNeighbors(seurat, dims = 1:20)

seurat <- FindClusters(seurat, resolution = 1) ### Yields 25 clusters 


########## Change the resolution to check the number of clusters ##########

seurat <- FindClusters(seurat, resolution = 0.2) #### Yields 16 clusters 

seurat <- FindClusters(seurat, resolution = 0.1) ####### Yields 13 clusters 


DimPlot(seurat, reduction = "umap", label = TRUE)

######### Final resolution: 0.1 ##################


################ STEP-15: FIND MARKERS OF EACH CLUSTER ##################

cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10_cl_markers

write.csv(cl_markers, "Results/markers.csv")

saveRDS(seurat, "seurat.rds")

############# Convert ensemble ids to gene symbols ###########

##### Direct vector match to map Ensembl IDs -> Symbols

cl_markers$symbol <- seurat[["RNA"]]@meta.data$symbol[match(cl_markers$gene, rownames(seurat))]

####### Get top 10 markers per cluster (now including the symbol column) #####

top10_cl_markers <- cl_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

###### View results with readable symbols ########

top10_cl_markers


########### Update the rownames of seurat object as per gene symbols #########



#################### STEP-16: ANNOTATE CELL CLUSTERS ########################

###### Visualize some cell markers to see their expression in clusters ##########

######## Check the expression of each marker with known literature ########

####### Original Study: https://pmc.ncbi.nlm.nih.gov/articles/PMC9050543/ ####

##### Another: https://doi.org/10.1101/2025.08.26.672469 #######

######### Astrocytes: AQP4 & GFAP ###################
DefaultAssay(seurat) <- "RNA"
FeaturePlot(seurat, c("AQP4","GFAP"), ncol = 1)

VlnPlot(seurat, features = c("AQP4","GFAP"), pt.size = 0)
grep("AQP", rownames(seurat), value = TRUE)
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

colnames(seurat)

