
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

##################### STEP-1: IDENTIFYING LncRNAS #######################

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

avg_expr <- AverageExpression(seurat, features = lncrna_genes, 
                              group.by = "seurat_clusters")

write.csv(avg_expr, "Results/lncRNA_expression.csv")

sum(duplicated(rownames(avg_expr))) ### Check for duplicated gene names

typeof(avg_expr) ### Check the type of avg_expr


avg_expr_num <- do.call(rbind, avg_expr)


# avg_expr: rows = lncRNA genes, columns = clusters, linear-scale mean expression

######## STEP 3: Compare expression — standard DE, restricted to lncRNAs ########

lncrna_markers <- FindAllMarkers(seurat, features = lncrna_genes, 
                                 only.pos = TRUE, min.pct = 0.1, 
                                 logfc.threshold = 0.25)

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


################################################################################

################# STEP-2: FIND DIFFERENTIALLY EXPRESSED LNCRNAS ###########

######### Aggregate counts by Patient ONLY ##########

pseudo_global <- AggregateExpression(
  seurat,
  group.by = "patient",
  return.seurat = TRUE
)

####### Add Condition Information ##########

pseudo_global$condition <- ifelse(grepl("PD", pseudo_global$patient), "PD", "Control")
Idents(pseudo_global) <- "condition"

######### Global DESeq2 DE on lncRNAs ########
global_de <- FindMarkers(
  object = pseudo_global,
  ident.1 = "PD",
  ident.2 = "Control",
  test.use = "DESeq2",
  features = rownames(avg_expr_num),
  min.pct = 0,
  logfc.threshold = 0
)

EnhancedVolcano(
  global_de,
  lab = rownames(global_de),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 4.0,
  title = 'lncRNAs: PD vs. Control',
  subtitle = 'Volcano plot DE lncRNAs (PD vs Control)'
)

write.csv(global_de, "Results/DElncRNAs_PD_vs_Control.csv")



################ DE Analysis per cluster #################


######## Convert Seurat to SingleCellExperiment ########

sce_pb <- as.SingleCellExperiment(seurat, assay = "RNA")
sce_pb$cluster_id <- seurat$cell_ontology
sce_pb$sample_id  <- seurat$patient
sce_pb$group_id   <- seurat$condition

######## Restrict to lncRNAs  ########

sce_pb <- sce_pb[intersect(lncrna_genes, rownames(sce_pb)), ]

######## Prep + aggregate to donor x cluster pseudobulk counts ########

sce_pb <- prepSCE(sce_pb, kid = "cluster_id", sid = "sample_id", gid = "group_id", drop = TRUE)
pb <- aggregateData(sce_pb, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

######## Run DS analysis: DESeq2 backend, one call, all clusters at once ########

res_muscat <- pbDS(pb, method = "DESeq2")

######## Extract per-cluster results as a single tidy table ########

muscat_results <- resDS(sce_pb, res_muscat)

write.csv(muscat_results, "Results/DElncRNAs_all_clusters_muscat.csv")

######## Filter muscat results to significant hits only, keep cluster info ########

sig_de_by_cluster <- muscat_results[!is.na(muscat_results$p_adj.loc) & 
                                      muscat_results$p_adj.loc < 0.05, ]

######## Keep just the useful columns, sorted by cluster then significance ########

sig_de_by_cluster <- sig_de_by_cluster[order(sig_de_by_cluster$cluster_id, 
                                             sig_de_by_cluster$p_adj.loc), 
                                       c("gene", "cluster_id", "logFC", "p_val", "p_adj.loc")]

head(sig_de_by_cluster, 20)

######## How many significant lncRNAs per cluster ########

table(sig_de_by_cluster$cluster_id)

######## Save it ########

write.csv(sig_de_by_cluster, "Results/DElncRNAs_significant_muscat.csv", row.names = FALSE)

########## DE lncRNAs per cluster #############

de_counts_per_cluster <- as.data.frame(table(sig_de_by_cluster$cluster_id))
colnames(de_counts_per_cluster) <- c("cluster", "n_significant_lncRNAs")
de_counts_per_cluster <- de_counts_per_cluster[order(-de_counts_per_cluster$n_significant_lncRNAs), ]

de_counts_per_cluster

######## Bar plot of DE lncRNA counts per cluster ########

ggplot(de_counts_per_cluster, aes(x = reorder(cluster, -n_significant_lncRNAs),
                                  y = n_significant_lncRNAs)) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell type", y = "Number of significant DE lncRNAs (padj < 0.05)",
       title = "DE lncRNAs per Cluster (PD vs Control)")


######### Check the expression of specific lncRNAs known to be involved in PD ####

###### NEAT1 Expression across clusters ###########

neat1_expr <- AverageExpression(seurat, features = "NEAT1", group.by = "cell_ontology")$RNA

###### Sorting ############

sort(neat1_expr["NEAT1", ], decreasing = TRUE)


df_neat1 <- data.frame(CellType = colnames(neat1_expr), NEAT1 = as.vector(neat1_expr))

# 2. Plot
ggplot(df_neat1, aes(x = reorder(CellType, -NEAT1), y = NEAT1)) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type", y = "Average NEAT1 Expression", 
       title = "NEAT1 Expression Across Cell Types")



######## Check the expression of multiple lncRNA ############

target_lncs <- intersect(c("NEAT1", "MALAT1", "MEG3", "MIAT"), rownames(seurat))

#### Get average expression and reshape into long format for ggplot #####

df_multi <- as.data.frame(AverageExpression(seurat, features = target_lncs, group.by = "cell_ontology")$RNA)
df_multi$Gene <- rownames(df_multi)
df_long <- reshape2::melt(df_multi, id.vars = "Gene", variable.name = "CellType", value.name = "Expression")

####### Plot Grouped Bar Plot ########

ggplot(df_long, aes(x = reorder(CellType, -Expression), y = Expression, fill = Gene)) +
  geom_col(show.legend = FALSE) +
  
  # Separate panel per lncRNA with individual Y-axis scales
  facet_wrap(~Gene, scales = "free_y", ncol = 2) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"), # Header panel style
    strip.text = element_text(face = "bold.italic", size = 11)        # Bold italic lncRNA names
  ) +
  labs(
    x = "Cell Type", 
    y = "Average Expression", 
    title = "Expression of Canonical lncRNAs Across Cell Types"
  )


################################################################################

################### STEP-3: VISUALIZATION ##########################

######### DotPlot to visualize lncRNAs from each cluster ##########

###### Identify top 2-3 specific lncRNAs per cluster #########

specific_expr_matrix <- avg_expr_num * tau_scores

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
  cluster_active[1:10] 
})

target_genes <- unique(as.vector(top_tau_genes))

######## Map Ensembl IDs to Gene Symbols ##########

symbol_lookup <- setNames(biotype_map$SYMBOL, biotype_map$GENEID)
gene_symbols <- symbol_lookup[target_genes]
gene_symbols[is.na(gene_symbols)] <- target_genes[is.na(gene_symbols)]

######### Dot Plot showing specific lncRNAs per cluster per tau index #######

DotPlot(seurat, features = target_genes, group.by = "cell_ontology") +
  scale_x_discrete(labels = gene_symbols) +
  
  # Add horizontal and vertical dotted guidelines
  geom_hline(yintercept = seq_len(length(unique(seurat$cell_ontology))) + 0.5, 
             color = "grey80", linetype = "dotted") +
  geom_vline(xintercept = seq_len(length(target_genes)) + 0.5, 
             color = "grey80", linetype = "dotted") +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "italic"),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed") # Panel grid fallback
  ) +
  labs(
    title = "Top Cell-Type Specific lncRNAs (Strictly by Tau Index)",
    x = "lncRNAs",
    y = "Cell Types/Clusters"
  )


saveRDS(seurat, "seurat.rds")

# install.packages("renv")
# 
# renv::init()
# 
renv::snapshot()

########################### END OF ANALYSIS ####################################
