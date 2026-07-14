# Predicting Circular RNAs in Parkinsons Disease (PD)

![Status: In Progress](https://img.shields.io/badge/status-in--progress-orange)

Predicting cell-type-specific circular RNA splicing in Parkinson's Disease (PD) and Parkinson's Disease Dementia (PDD) from single-nucleus RNA-seq data using Seurat, Harmony, and CIRI-deepA.

## Workflow Overview
1. **Metadata Curation:** Classify 34 human brain donor samples.
2. **Preprocessing & QC:** Merge matrices, apply strict nuclear QC filters, and cluster using Seurat.
3. **Integration:** Address potential batch effects using Harmony.
4. **Prediction:** Extract RNA-Binding Protein (RBP) profiles and run CIRI-deep to predict differential circRNA splicing.
