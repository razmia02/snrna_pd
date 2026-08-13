# Identifying Non-Coding RNAs in Parkinsons Disease (PD)

![Status: In Progress](https://img.shields.io/badge/status-in--progress-orange)

Analysis of snRNA-Seq to profile long non-coding RNAs (lncRNAs) and predict cell-type-specific circular RNA (circRNA) splicing in Parkinson's Disease (PD) and Parkinson's Disease Dementia (PDD) from NCBI GEO dataset [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783).

## Workflow Overview
1. **Metadata Curation:** Classify 11 human brain donor samples.
2. **Preprocessing & QC:** Merge matrices, apply strict nuclear QC filters (mitochondrial reads < 1%), and cluster using Seurat.
3. **Clustering:** Cell type clustering and annotation.
4. **lncRNA Profiling:** Directly extract and analyze physically captured, nuclear-enriched lncRNA expression across cell clusters and calculate cell-type specificity index. 
5. **circRNA Prediction:** Leverage the same cell clusters to extract RNA-Binding Protein (RBP) profiles, using CIRI-deepA to predict differential circRNA splicing.

## 🔬 Methodology

### 1. Data Acquisition & Preprocessing
* **Dataset**: Single-nucleus RNA sequencing (snRNA-seq) data from human post-mortem Substantia Nigra (`GSE157783`).
* **Matrix Formatting**: Raw expression count matrices were converted into sparse matrices (`dgCMatrix`) in R.
* **Gene Annotation**: Ensembl gene IDs were converted to official HGNC Gene Symbols using the Ensembl database (`EnsDb.Hsapiens.v86`). Duplicate mappings were resolved using `make.unique()`.

---

### 2. Quality Control & Filtering
* **Filtering Thresholds**: High-quality nuclei were retained based on published quality control criteria:
  * Unique Molecular Identifiers (UMIs / `nCount_RNA`) > 1,500
  * Detected Genes (`nFeature_RNA`) > 1,000
  * Mitochondrial gene transcript fraction (`percent.mt`) < 10%

---

### 3. Doublet Identification & Removal
* Doublets were detected using `scDblFinder` by converting the Seurat object to a `SingleCellExperiment` format.
* Doublet scores and class labels (`singlet` vs. `doublet`) were computed per patient sample, and all predicted doublets were removed prior to downstream analysis.

---

### 4. Normalization, Feature Selection & Dimensionality Reduction
* **Normalization**: Expression counts were log-normalized using standard Seurat library-size normalization (`NormalizeData`).
* **Variable Features**: The top 2,000 highly variable features were identified using `FindVariableFeatures`.
* **Scaling & PCA**: Data were scaled (`ScaleData`), and Principal Component Analysis (PCA) was performed for the top 50 PCs.
* **Dimensionality**: Based on elbow plot heuristics, the top 20 Principal Components were selected for downstream integration and visualization.

---

### 5. Batch Correction & Clustering
* **Batch Correction**: Integration across individual patient samples was performed using `Harmony` (`RunHarmony`) on the top 20 PCs to mitigate sample-to-sample batch effects.
* **Clustering**: Graph-based clustering was applied (`FindNeighbors` on Harmony embeddings) across multiple resolution parameters ($0.1, 0.2, 1.0$). A final resolution of **0.1** was selected, yielding 13 robust cell clusters.
* **Visualization**: Uniform Manifold Approximation and Projection (UMAP) was used for non-linear dimensional reduction.

---

### 6. Cell Type Annotation
* Clusters were manually annotated using canonical cell-type specific marker genes from published midbrain literature:
  * **Astrocytes**: *AQP4*, *GFAP*
  * **Oligodendrocytes**: *MOBP*, *MBP*, *MOG*, *OPALIN*
  * **Oligodendrocyte Precursor Cells (OPCs)**: *PDGFRA*, *VCAN*, *PCDH15*
  * **Microglia**: *CD74*, *APBB1IP*
  * **Endothelial Cells**: *CLDN5*
  * **Pericytes**: *PDGFRB*, *RGS5*, *NOTCH3*, *ABCC9*
  * **Excitatory Neurons**: *SLC17A6*
  * **Inhibitory / GABAergic Neurons**: *GAD1*, *GAD2*, *GRIK1*
  * **Ependymal Cells**: *FOXJ1*
  * **T-lymphocytes**: *ITK*, *IL7R*, *CD96*

---

### 7. Long Non-Coding RNA (lncRNA) Identification & Cell-Type Specificity
* **lncRNA Filtering**: Transcripts were annotated via `EnsDb.Hsapiens.v86` and filtered for lncRNA biotypes (`lincRNA`, `antisense`, `sense_overlapping`, `processed_transcript`), isolating **9,073 lncRNAs**.
* **Cell-Type Specificity Score ($\tau$ Index)**: Calculated to quantify lncRNA cell-type specificity across clusters using average linear expression values:
  $$\tau = \frac{\sum_{i=1}^{N} \left(1 - \frac{x_i}{x_{\max}}\right)}{N - 1}$$
* **Specific Expression Scoring**: LncRNAs were ranked by multiplying mean expression values by their calculated Tau ($\tau$) indices.

---

### 8. Differential Expression Analysis (PD vs. Control)
* Differential expression was calculated using Wilcoxon Rank-Sum tests (`FindMarkers`) with thresholds set to $\text{log2FC} > 0.25$ and expression present in $\ge 10\%$ of cells (`min.pct = 0.1`).
* Analyses were conducted across two scopes:
  1. **Global DE**: Overall expression shifts in Parkinson's Disease vs. Healthy Controls across all lncRNAs.
  2. **Cell-Type Specific DE**: Condition-wise comparison (e.g., `Microglia_PD` vs. `Microglia_Control`, `Astrocytes_PD` vs. `Astrocytes_Control`, `Pericytes_PD` vs. `Pericytes_Control`).