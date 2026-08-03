# Identifying Non-Coding RNAs in Parkinsons Disease (PD)

![Status: In Progress](https://img.shields.io/badge/status-in--progress-orange)

Analysis of snRNA-Seq to profile long non-coding RNAs (lncRNAs) and predict cell-type-specific circular RNA (circRNA) splicing in Parkinson's Disease (PD) and Parkinson's Disease Dementia (PDD) from NCBI GEO dataset [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783).

## Workflow Overview
1. **Metadata Curation:** Classify 11 human brain donor samples.
2. **Preprocessing & QC:** Merge matrices, apply strict nuclear QC filters (mitochondrial reads < 1%), and cluster using Seurat.
3. **Clustering:** Cell type clustering and annotation.
4. **lncRNA Profiling:** Directly extract and analyze physically captured, nuclear-enriched lncRNA expression across cell clusters and calculate cell-type specificity index. 
5. **circRNA Prediction:** Leverage the same cell clusters to extract RNA-Binding Protein (RBP) profiles, using CIRI-deepA to predict differential circRNA splicing.