# scATACseq Analysis in ArchR

This repository contains the analysis pipeline and results for a single-cell ATAC-seq (scATAC-seq) project focused on chromatin accessibility profiling across different neuronal cell types.

## Project Overview

This project involves comprehensive analysis of scATAC-seq data to characterize chromatin accessibility landscapes across neuronal populations. The analysis includes quality control, dimensionality reduction, clustering, peak calling, motif enrichment, and integration with RNA-seq data to identify cell type-specific regulatory elements.

## Data Summary

After quality control:
- Total cells analyzed: 8,149
- Samples: dc1r3_r1, dc2r2_r1, dc2r2_r2
- Median TSS enrichment: 19.838
- Median fragments per cell: 10,432
- Peak matrix dimensions: 6,062,095 × 8,149

## Analysis Workflow

### 1. Quality Control and Preprocessing

Initial dataset contained 9,827 cells, which was reduced to 8,149 after QC filtering. QC metrics included:
- TSS enrichment scores
- Fragment counts
- Fragment length distribution

The dc2r2_r2 sample showed potential quality issues with lower cell counts and some cells with very low TSS enrichment scores.

### 2. Dimensionality Reduction

- Used Latent Semantic Indexing (LSI) instead of PCA due to the extreme sparsity of scATAC-seq data
- Applied Harmony batch correction to integrate data across samples
- Generated UMAP embeddings for visualization

### 3. Clustering Analysis

- Identified 24 distinct clusters
- Characterized each cluster by cell count and sample composition
- Generated cluster-specific accessibility profiles

### 4. Peak Calling and Analysis

- Created aggregate accessibility profiles for each cluster
- Generated reproducible peak sets using a tile-based method
- Constructed cells × peaks matrix for differential accessibility analysis
- Parameters used for differential analysis:
  - Matrix: PeakMatrix
  - Grouping: Clusters
  - Bias correction: TSSEnrichment, log10(nFrags)
  - Test method: Wilcoxon rank-sum test
  - Max cells: 500 per group

### 5. Gene Activity Scoring

- Generated gene activity scores based on chromatin accessibility
- Identified marker genes for each cluster using:
  - FDR ≤ 0.01
  - Log2FC ≥ 1.25
- Applied MAGIC imputation to enhance visualization of gene activity patterns

### 6. Transcription Factor Motif Analysis

- Used the Cis-BP motif database for TF binding site annotation
- Scanned peak regions for potential TF binding motifs
- Identified cluster-specific enriched motifs
- Visualized motif activities across clusters

### 7. Integration with RNA-seq Data

- Integrated ATAC-seq with matched RNA-seq data
- Evaluated correlation between chromatin accessibility and gene expression
- Key findings:
  - Major populations like nIPC/GluN1 were well preserved in integration
  - Rare populations (GluN5, SP) were underrepresented after integration
  - Highest agreement genes: DLX6-AS1, DLX6, SLCO1C1, SATB2, ST8SIA5
  - Lowest agreement genes: KLC1, TACR1, TMED4, LINC01103, DLGAP1-AS2

### 8. Cell Type-Specific Analyses

Special focus on two populations:
- GluN5 cells
- Cycling Progenitor cells

For these populations:
- Performed differential accessibility analysis
- Identified enriched TF motifs:
  - GluN5: GLI2, WT1, TFAP2C, SP4, XBP1
  - Cycling Progenitors: BHLHA15, TAL1, TAL2, NEUROG1, ATOH7
- Generated TF motif activity visualizations

### 9. TF Footprinting Analysis

- Applied Tn5 bias correction using "Divide" normalization method
- Generated aggregate footprints for GluN5 and Cycling progenitor cells
- Observed distinct binding patterns between cell types, with Cycling Progenitors showing more pronounced peaks

### 10. Co-accessibility and Peak-to-Gene Linkage

- Computed co-accessibility between peaks
- Identified potential enhancers linked to cell type-specific marker genes
- Visualized genomic tracks showing peak-to-gene linkages

## Technical Implementation

All analyses were performed using the ArchR package in R, which provides a comprehensive framework for scATAC-seq data analysis. Key ArchR functions used include:
- `addGroupCoverages`
- `addReproduciblePeakSet`
- `addPeakMatrix`
- `getMarkerFeatures`
- `addBgdPeaks`
- `addCoAccessibility`

## Future Directions

Potential next steps include:
1. Integration with additional omics data (e.g., scMultiome)
2. Validation of putative enhancers using reporter assays
3. Deeper analysis of cell type-specific TF networks
4. Targeted CRISPR experiments to validate regulatory elements

## Contact

For questions or collaborations, please contact Zain at zazi00001@stud.uni-saarland.de.
