
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

######################## STEP-3: CREATE SEURAT OBJECT #########################

seurat <- CreateSeuratObject(counts_matrix, project="PD")

seurat

dim(seurat)

rownames(seurat)

Layers(seurat)

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

seurat@meta.data$patient

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

seurat <- FindClusters(seurat, resolution = 1) ### Yields 24 clusters 


########## Change the resolution to check the number of clusters ##########

seurat <- FindClusters(seurat, resolution = 0.2) #### Yields 15 clusters 

seurat <- FindClusters(seurat, resolution = 0.1) ####### Yields 13 clusters 


DimPlot(seurat, reduction = "umap", label = TRUE)

######### Final resolution: 0.1 ##################


################ STEP-15: FIND MARKERS OF EACH CLUSTER ##################

cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10_cl_markers

write.csv(cl_markers, "Results/markers.csv")

saveRDS(seurat, "seurat.rds")


############## Some genes have NA symbol #####################
#################### STEP-16: ANNOTATE CELL CLUSTERS ########################

###### Visualize some cell markers to see their expression in clusters ##########

######## Check the expression of each marker with known literature ########

####### Original Study: https://pmc.ncbi.nlm.nih.gov/articles/PMC9050543/ ####

##### Another: https://doi.org/10.1101/2025.08.26.672469 #######

######### Astrocytes: AQP4 & GFAP ###################

FeaturePlot(seurat, c("AQP4","GFAP"), ncol = 1)

VlnPlot(seurat, features = c("SLC6A3", "ALDH1A1", "NR4A2", "SLC18A2"), pt.size = 0)

#### AQP4: Cluster 2 & 9
#### GFAP: Cluster 2,5,6,7,9 & 11

######### Oligodendrocytes: MBP, MOBP, MOG & OPALIN ###########
#### MOG & MOBP: Clusters 0,1 & 12
#### MBP: Clusters all except 9 & 11
#### OPALIN: Cluster 0

########### Oligodendrocyte Precursor Cells (OPCs): VCAN, PDGFRA, PCDH15 ###########
##### VCAN: Clusters 2,5 & 9
##### PDGFRA: Clusters 5
##### PCDH15: Clusters 3,5,8 & 10

############ Ependymal Cells: FOXJ1 ###############
#### FOXJ1: Cluster 9

############# Microglia: CD74 & APBB1IP #################
#### CD74: Clusters 4,6 & 12
#### APBB1IP: Clusters 4 & 12

################ Endothelial Cells: CLDN5 ###################
##### CLDN5: Clusters 6 & 7

################## Excitatory Neurons: SLC17A6 ##############
##### SLC17A6: Cluster 3 & 10

################### Inhibitory Neurons: GAD1 & GAD2 ################
##### GAD1: Clusters 3,5 & 8
##### GAD2: Clusters 3 & 8

#################### GABAergic Neurons: GAD2 & GRIK1 ###############
##### GAD2: Clusters 3 & 8
##### GRIK1: Clusters 3,5,8 & 10

################## Dopaminergic Neurons: TH #################
##### TH: Not detectable in any clusters. 



########### Find markers of unambiguous clusters ##################

cluster3.markers <- FindMarkers(seurat,
                                ident.1 = 3,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)

head(cluster3.markers, 50)


cluster8.markers <- FindMarkers(seurat,
                                ident.1 = 8,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)

head(cluster8.markers, 50) ### Inhibitory neurons markers

markers_3v10 <- FindMarkers(seurat, ident.1 = 3, ident.2 = 10, 
                            min.pct = 0.25, logfc.threshold = 0.25)
head(markers_3v10, 30)

nrow(markers_3v10[markers_3v10$p_val_adj < 0.05, ])  # how many genes actually distinguish them

########## Cluster 3 & 11 are both excitatory neurons but different subtypes ######

############## ANNOTATE CLUSTERS BASED ON MARKER EXPRESSION #############

markers <- c(
  "AQP4", "GFAP",          # Astrocytes
  "MOBP", "MBP", "MOG", "OPALIN",    # Oligodendrocytes
  "PDGFRA", "VCAN", "PCDH15",  # OPCs
  "FOXJ1",                  # Ependymal cells
  "CD74", "APBB1IP",       # Microglia
  "CLDN5",                 # Endothelial
  "SLC17A6",  # Excitatory Neurons
  "GAD1", "GAD2", "OTX2-AS1",  # Inhibitory Neurons
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
# 3 : Excitatory Neurons (SubtypeA)
# 4 : Microglia
# 5 : OPCs
# 6 : Endothelial Cells
# 7 : Pericytes
# 8 : Inhibitory Neurons
# 9 : Ependymal cells
# 10 : Excitatory Neurons (SubtypeB)
# 11 : CADPS2+
# 12 : T-lymphocyte


new_ident <- setNames(c("OPALIN+ Oligodendrocytes",
                        "Oligodendrocytes",
                        "Astrocytes",
                        "Excitatory Neurons (SubtypeA)",
                        "Microglia",
                        "Oligodendrocytes precursor cells (OPCs)",
                        "Endothelial Cells",
                        "Pericytes",
                        "Inhibitory Neurons",
                        "Ependymal Cells",
                        "Excitatory Neurons (SubtypeB)",
                        "CADPS2+",
                        "T lymphocytes"),
                      levels(seurat))

seurat$cell_ontology <- Idents(seurat)

DimPlot(seurat, reduction = "umap", label = TRUE)


##################### STEP-17: IDENTIFYING LncRNAS #######################

######## Get lncRNA gene list ####################

biotype_map <- ensembldb::select(EnsDb.Hsapiens.v86, keys = rownames(seurat),
                                 keytype = "SYMBOL", columns = "GENEBIOTYPE")

table(biotype_map$GENEBIOTYPE) ##### Check which type of molecules are in data

######### Define the types of lncRNAs present in data by a vector #############

lnc_biotypes <- c("lincRNA", "antisense",
                  "sense_overlapping", "processed_transcript")


######## Filter out the selected lncRNAs ############

lncrna_genes <- unique(biotype_map$SYMBOL[biotype_map$GENEBIOTYPE %in% lnc_biotypes])

######### Check the lncRNAs in dataset ##################

lncrna_genes <- intersect(lncrna_genes, rownames(seurat))

length(lncrna_genes) ### 9073 lnRNAs in total

######## Average expression per cluster, restricted to lncRNAs ########
?AverageExpression
remove(avg_expr)
avg_expr <- AverageExpression(seurat, features = lncrna_genes, 
                              group.by = "seurat_clusters")

write.csv(avg_expr, "Results/lncRNA_expression.csv")

sum(duplicated(rownames(avg_expr))) 

typeof(avg_expr) ### Check the type of avg_expr

sum(duplicated(rownames(avg_expr))) ### Check for duplicated gene names

avg_expr_num <- do.call(rbind, avg_expr)


# avg_expr: rows = lncRNA genes, columns = clusters, linear-scale mean expression

######## STEP 3: Compare expression — standard DE, restricted to lncRNAs ########

lncrna_markers <- FindAllMarkers(seurat, features = lncrna_genes, 
                                 only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

head(lncrna_markers)

write.csv(lncrna_markers, "Results/lncrna_markers.csv")


######## STEP 4: tau index — cluster-specificity score per lncRNA ########

tau_calc <- function(x) {
  if (max(x) == 0) return(0)
  n <- length(x)
  sum(1 - (x / max(x))) / (n - 1)
}


tau_scores <- apply(avg_expr_num, 1, tau_calc)

tau_df <- data.frame(gene = names(tau_scores), tau = tau_scores)

tau_df <- tau_df[order(-tau_df$tau), ]

head(tau_df, 20)   # most cluster-specific lncRNAs, ranked

write.csv(tau_df, "Results/tau_index.csv")



################# STEP-18: FIND DIFFERENTIALLY EXPRESSED LNCRNAS ###########

####### Set main identity to condition column ######

Idents(seurat) <- "condition"

######### Find differentially expressed lncRNAs between PD and Control ######

pd_vs_ctrl_markers <- FindMarkers(
  seurat,
  ident.1 = "PD",
  ident.2 = "Control",
  features = rownames(avg_expr_num), # Restrict to your lncRNAs
  logfc.threshold = 0.25,
  min.pct = 0.1
)

###### Look at top PD-upregulated lncRNAs #######

head(pd_vs_ctrl_markers[order(-pd_vs_ctrl_markers$avg_log2FC), ])

pd_vs_ctrl_markers

###### Volcano Plot to visualize DE lncRNAs ############

EnhancedVolcano(
  pd_vs_ctrl_markers,
  lab = rownames(pd_vs_ctrl_markers),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 4.0,
  title = 'lncRNAs: PD vs. Control',
  subtitle = 'Volcano plot DE lncRNAs (PD vs Control)'
)

write.csv(pd_vs_ctrl_markers, "Results/DE_lncRNAs_PD.csv")


############ CLUSTER TYPE SPECIFIC DE ANALYSIS ##############

######## Combine Cell Type and Condition in metadata #########

seurat$celltype_condition <- paste0(seurat$cell_ontology, "_", seurat$condition)

Idents(seurat) <- "celltype_condition"

######### PD vs Control specifically in Microglia ##########

pd_microglia_deg <- FindMarkers(
  seurat,
  ident.1 = "Microglia_PD",
  ident.2 = "Microglia_Control",
  features = rownames(avg_expr_num)
)

######### Volcano Plot #############

EnhancedVolcano(
  pd_microglia_deg,
  lab = rownames(pd_microglia_deg),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 4.0,
  title = 'lncRNAs: PD vs. Control',
  subtitle = 'Volcano plot showing DE lncRNAs in Microglia (PD vs Control)'
)


###### DE lncrnas in Astrocytes (PD vs. Control) ########

astrocyte_de <- FindMarkers(
  seurat,
  ident.1 = "Astrocytes_PD",
  ident.2 = "Astrocytes_Control",
  features = rownames(avg_expr_num),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

######## DE lncRNAs in Pericytes (PD vs. Control) ############

pericyte_de <- FindMarkers(
  seurat,
  ident.1 = "Pericytes_PD",
  ident.2 = "Pericytes_Control",
  features = rownames(avg_expr_num),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

EnhancedVolcano(
  astrocyte_de,
  lab = rownames(astrocyte_de),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 4.0,
  title = 'lncRNAs: PD vs. Control',
  subtitle = 'Volcano plot showing DE lncRNAs in Astrocytes (PD vs Control)'
)


EnhancedVolcano(
  pericyte_de,
  lab = rownames(pericyte_de),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 4.0,
  title = 'lncRNAs: PD vs. Control',
  subtitle = 'Volcano plot showing DE lncRNAs in Pericytes (PD vs Control)'
)

write.csv(pericyte_de, "Results/DElncRNAs_Pericytes.csv")
write.csv(pd_microglia_deg, "Results/DElncRNAs_Microglia.csv")
write.csv(astrocyte_de, "Results/DElncRNAs_Astrocytes.csv")

################### STEP-19: VISUALIZATION ##########################

######### DotPlot to visualize lncRNAs from each cluster ##########

###### Identify top 2-3 specific lncRNAs per cluster #########

top_lncs <- apply(specific_expr_matrix, 2, function(col) {
  names(sort(col, decreasing = TRUE, na.last = TRUE))[1:5]
})

target_genes <- unique(as.vector(top_lncs))

######## Map Ensembl IDs to Gene Symbols #########

symbol_lookup <- setNames(biotype_map$SYMBOL, biotype_map$GENEID)
gene_symbols <- symbol_lookup[target_genes]
gene_symbols[is.na(gene_symbols)] <- target_genes[is.na(gene_symbols)]

######### Dot Plot showing top lncRNAs per cluster/cell type #############

DotPlot(seurat, features = target_genes, group.by = "cell_ontology") +
  scale_x_discrete(labels = gene_symbols) +
  coord_flip() + # Flips plot so cell types are on X-axis and lncRNAs on Y-axis for easier reading
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic")
  ) +
  labs(
    title = "Top Specific lncRNAs Across Annotated Cell Types",
    x = "lncRNAs",
    y = "Cluster/Cell Type"
  )

############### DotPlot as per tau index per cluster ################

###### Select top 3-5 lncRNAs per cluster STRICTLY by highest Tau index ######

tau_order <- order(tau_scores, decreasing = TRUE, na.last = TRUE)
sorted_genes <- rownames(avg_expr_num)[tau_order]

top_tau_genes <- sapply(colnames(avg_expr_num), function(cluster_name) {
  cluster_active <- sorted_genes[avg_expr_num[sorted_genes, cluster_name] > 0]
  cluster_active[1:5] 
})

target_genes <- unique(as.vector(top_tau_genes))

######## Map Ensembl IDs to Gene Symbols ##########

symbol_lookup <- setNames(biotype_map$SYMBOL, biotype_map$GENEID)
gene_symbols <- symbol_lookup[target_genes]
gene_symbols[is.na(gene_symbols)] <- target_genes[is.na(gene_symbols)]

######### Dot Plot showing specific lncRNAs per cluster per tau index #######

DotPlot(seurat, features = target_genes, group.by = "cell_ontology") +
  scale_x_discrete(labels = gene_symbols) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic")
  ) +
  labs(
    title = "Top Cell-Type Specific lncRNAs (Strictly by Tau Index)",
    x = "Annotated Cell Type",
    y = "lncRNA Symbol"
  )


