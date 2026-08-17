############### R SCRIPT FOR CLUSTER ANALYSIS  ###########################

#### Clustering, Marker Identification, Annotation ####

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

################ STEP-1: FIND MARKERS OF EACH CLUSTER ##################

cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10_cl_markers

write.csv(cl_markers, "Results/markers.csv")

saveRDS(seurat, "seurat.rds")


############## Some genes have NA symbol #####################


################################################################################

#################### STEP-2: ANNOTATE CELL CLUSTERS ########################

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

seurat <- RenameIdents(seurat, new_ident)

seurat$cell_ontology <- Idents(seurat)

DimPlot(seurat, reduction = "umap", label = TRUE)

seurat@meta.data$cell_ontology

saveRDS(seurat, "seurat.rds")