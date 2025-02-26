
# setting seed for reproducibility
set.seed(11)

# loading necessary libraries
library(ArchR)
library(pheatmap)
library(BSgenome.Hsapiens.UCSC.hg38)

# adding genome
addArchRGenome("hg38")

# read the data into an appropriate 
# data structure and apply filtering
sampleFiles <- c("hft_ctx_w21_dc1r3_r1_atac_fragments.tsv.gz",
                 "hft_ctx_w21_dc2r2_r1_atac_fragments.tsv.gz",
                 "hft_ctx_w21_dc2r2_r2_atac_fragments.tsv.gz")

sampleNames <- c("dc1r3_r1", "dc2r2_r1", "dc2r2_r2")

ArrowFiles <- createArrowFiles(
  inputFiles = file.path("data", sampleFiles),
  sampleNames = sampleNames,
  minTSS = 4,
  minFrags = 500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# identify doublets
doubletScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# create a project that includes all samples. 
# Inspect the cell metadata.
archrproj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "save",
  copyArrows = TRUE
)

# 1.4 questions
print(archrproj)
tilematrix <- getMatrixFromProject(
    ArchRProj = archrproj,
    useMatrix="TileMatrix",
    binarize = TRUE)
print(dim(tilematrix))

# 1.5 questions
sample_cell_counts <- table(getCellColData(archrproj, select = "Sample"))
write.csv(sample_cell_counts, "sample_cell_counts.csv", row.names = FALSE)

# plot the fragment length distribution 
# of all samples in a single plot
flPlot <- plotFragmentSizes(
    ArchRProj = archrproj,
    groupBy = "Sample")

plotPDF(flPlot, 
    name = "Fragment-Length-Distribution.pdf",
    width = 8,
    height = 6,
    ArchRProj = archrproj,
    addDOC = FALSE)

# Plot the distribution of TSS 
# enrichment scores in each sample
TSSenrchplot <- plotTSSEnrichment(ArchRProj = archrproj)

plotPDF(TSSenrchplot, 
    name = "TSS-Enrichment-Distribution", 
    width = 8, 
    height = 6,
    ArchRProj = archrproj,
    addDOC = FALSE)

# stricter QC
archrproj <- filterDoublets(archrproj)

make_qcplot <- function(archrproj, name) {
    df <- getCellColData(archrproj, select = c("log10(nFrags)", "TSSEnrichment"))
    qc_p <- ggPoint(
        x = df[,1], 
        y = df[,2], 
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
    )

    plotPDF(qc_p, name = name, ArchRProj = archrproj, addDOC = FALSE)
}

make_qcplot(archrproj, "noQC_ArchRProj")

idxPass <- which(
    archrproj$TSSEnrichment >= 12 &
    archrproj$TSSEnrichment <= 27 &
    log10(archrproj$nFrags) >= 3.3 &
    log10(archrproj$nFrags) <= 4.6
)

cellsPass <- archrproj$cellNames[idxPass]
archrproj <- archrproj[cellsPass, ]
make_qcplot(archrproj, "QC_ArchRProj")

# 1.4 questions
print(archrproj)
tilematrix_afterqc <- getMatrixFromProject(
    ArchRProj = archrproj,
    useMatrix="TileMatrix",
    binarize = TRUE)
print(dim(tilematrix_afterqc))

saveArchRProject(
    ArchRProj = archrproj,
    outputDirectory = "save",
    load = FALSE
)

# Iterative LSI
archrproj <- addIterativeLSI(
    ArchRProj = archrproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 3, 
    clusterParams = list(
        resolution = c(0.4, 0.8, 1.2),
        sampleCells = 10000, 
        n.start = 20), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

# add umap
archrproj <- addUMAP(
    ArchRProj = archrproj,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = archrproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", size = 1.25)
p2 <- plotEmbedding(ArchRProj = archrproj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", size = 1.25)
p3 <- plotEmbedding(ArchRProj = archrproj, colorBy = "cellColData", name = "nFrags", embedding = "UMAP", size = 1.25)

plotPDF(p1, p2, p3, name = "UMAP_Sample_TSS_nfrags_NoBC.pdf", ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)

# dealing with batch effects
archrproj <- addHarmony(
    ArchRProj = archrproj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

archrproj <- addUMAP(
    ArchRProj = archrproj,
    reducedDims = "Harmony",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
)

p1 <- plotEmbedding(ArchRProj = archrproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", size = 1.25)
p2 <- plotEmbedding(ArchRProj = archrproj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", size = 1.25)
p3 <- plotEmbedding(ArchRProj = archrproj, colorBy = "cellColData", name = "nFrags", embedding = "UMAP", size = 1.25)

plotPDF(p1, p2, p3, name = "UMAP_Sample_TSS_nfrags.pdf", ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
    ArchRProj = archrproj,
    outputDirectory = "save",
    load = FALSE
)

# clustring
archrproj <- addClusters(
    input = archrproj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = c(0.8, 1.2, 1.6),
    force = TRUE
)

p1 <- plotEmbedding(ArchRProj = archrproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.25)
plotPDF(p1, name = "Clusters_UMAP.pdf", ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)


# How many cells does each cluster contain?
cell_count <- table(archrproj$Clusters)
write.csv(cell_count, "cell_count.csv", row.names = FALSE)

# What are the sample proportions in each cluster?
cM <- confusionMatrix(paste0(archrproj$Clusters), paste0(archrproj$Sample))
cM <- cM / Matrix::rowSums(cM)
p1 <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteRed"), 
    border_color = "black"
)
plotPDF(p1, name = "Sample_Proportions.pdf", ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)

cluster_props_pct <- cM * 100
write.csv(cluster_props_pct, "cluster_props_pct.csv", row.names = TRUE)

# peak calling using ArchR native peak caller
# as MACS2 is not supported on windows
archrproj <- addGroupCoverages(
    ArchRProj = archrproj,
    groupBy = "Clusters"
)

dir_path <- "save/PeakCalls"
if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, mode = "0777")}

archrproj <- addReproduciblePeakSet(
    ArchRProj = archrproj, 
    groupBy = "Clusters",
    peakMethod = "Tiles",
    method = "p"
)

archrproj <- addPeakMatrix(
    ArchRProj = archrproj
)

# cluster marker peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = archrproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    maxCells = 500
)

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = TRUE,
  labelRows = FALSE,
  nLabel = 0,
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = archrproj, addDOC = FALSE)

p <- plotBrowserTrack(
    ArchRProj = archrproj, 
    groupBy = "Clusters", 
    geneSymbol = c('TOP2A', 'MKI67', 'AURKA', 'SATB2', 'SLC12A7'),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)

plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = archrproj, addDOC = FALSE)

# gene activity
# identiy marker genes
markerGenes <- getMarkerFeatures(
    ArchRProj = archrproj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    maxCells = 500
)

heatmapGenes <- markerHeatmap(
  seMarker = markerGenes, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  nLabel = 1,
  transpose = TRUE
)

plotPDF(heatmapGenes, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = archrproj, addDOC = FALSE)

p1 <- plotEmbedding(
    ArchRProj = archrproj,
    colorBy = "GeneScoreMatrix",
    name = c('TOP2A', 'MKI67', 'AURKA', 'SATB2', 'SLC12A7'),
    embedding = "UMAP",
    size = 1.25
)

plotPDF(
    p1,
    name = "Gene-Score-UMAPs",
    width = 5,
    height = 5,
    ArchRProj = archrproj,
    addDOC = FALSE
)

# magic imputation
archrproj <- addImputeWeights(archrproj, reducedDims = 'Harmony')

p1 <- plotEmbedding(
    ArchRProj = archrproj,
    colorBy = "GeneScoreMatrix",
    name = c('TOP2A', 'MKI67', 'AURKA', 'SATB2', 'SLC12A7'),
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archrproj)
)

plotPDF(
    p1,
    name = "Gene-Score-UMAPs_withMAGIC",
    width = 5,
    height = 5,
    ArchRProj = archrproj,
    addDOC = FALSE
)

# compute TF motif activity
archrproj <- addMotifAnnotations(ArchRProj = archrproj, motifSet = "cisbp", name = "Motif")
motifsUp <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = archrproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
p1 <- plotEnrichHeatmap(motifsUp, n = 6, transpose = TRUE)
plotPDF(p1, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = archrproj, addDOC = FALSE)

df <- data.frame(
    TF = rownames(motifsUp),
    score = rowSums(assay(motifsUp))
)

df <- df[order(df$score, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
top2_df <- head(df, 2)
top2_motifs <- top2_df$TF
print(top2_df)

archrproj <- addBgdPeaks(archrproj)
archrproj <- addDeviationsMatrix(
    ArchRProj = archrproj,
    peakAnnotation = "Motif",
    force = TRUE
)

p1 <- plotEmbedding(
    ArchRProj = archrproj,
    colorBy = "MotifMatrix",
    name = paste0("z:", top2_motifs[1]),
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archrproj)
)

p2 <- plotEmbedding(
    ArchRProj = archrproj,
    colorBy = "MotifMatrix",
    name = paste0("z:", top2_motifs[2]),
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archrproj)
)

plotPDF(
    p1, p2,
    name = "top2markermotifs",
    width = 5,
    height = 5,
    ArchRProj = archrproj,
    addDOC = FALSE
)

saveArchRProject(
    ArchRProj = archrproj,
    outputDirectory = "save",
    load = FALSE
)

#integration
seRNA <- readRDS("data\\new_pbmc.rds")

archrproj <- addGeneIntegrationMatrix(
    ArchRProj = archrproj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seRNA,
    addToArrow = TRUE,
    groupRNA = "Cluster.Name",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE,
    nGenes = 2000,
    dimsToUse = 1:30
    )

p <- plotEmbedding(
    archrproj, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un",
    size = 0.8
)

plotPDF(p, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = archrproj, addDOC = FALSE, embedding = "UMAP", width = 5, height = 5)
plotPDF(p, name = "GluN5.pdf", ArchRProj = archrproj, addDOC = FALSE, embedding = "UMAP", width = 5, height = 5)

# plot marker genes in UMAP
genes <- c("ID4", "NFIA", "EGR1", "FOS", "ASCL1", "HES5", 
           "OLIG2", "MEIS2", "NHLH1", "SOX21", "SOX10", 
           "PBX1", "NEUROD1", "NEUROG2", "EOMES")

p1 <- plotEmbedding(
    ArchRProj = archrproj, 
    colorBy = "GeneIntegrationMatrix", 
    name = genes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archrproj)
)

plotPDF(p1, name = "highlighted-marker-genes.pdf", ArchRProj = archrproj, addDOC = FALSE, width = 5, height = 5)

# Correlation Coefficients
corr_mat <- correlateMatrices(
    ArchRProj = archrproj,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "GeneIntegrationMatrix",
    reducedDims = "Harmony"
)

cmat <- corr_mat[!is.na(corr_mat$pval) & corr_mat$pval < 0.05, ]
cmat <- cmat[order(cmat$cor, decreasing = TRUE), ]
cmat <- cmat[-2,] #duplicated entry?
cmat <- cmat[-2,] #duplicated entry?
top5 <- head(cmat, 5)
bottom5 <- tail(cmat, 5)
combined_results <- rbind(top5, bottom5)
write.csv(combined_results, "top_bottom_5_correlations.csv", row.names = FALSE)

# Cluster labels from gene expression
cM <- confusionMatrix(archrproj$Clusters, archrproj$predictedGroup_Un)
labelOld <- rownames(cM)
labels <- colnames(cM)[apply(cM, 1, which.max)]
labels[3] <- "GluN5"
archrproj$Clusters2 <- mapLabels(archrproj$Clusters, newLabels = labels, oldLabels = labelOld)
write.csv(cM, "cM.csv", row.names = TRUE)

p1 <- plotEmbedding(
    ArchRProj = archrproj,
    colorBy = "cellColData",
    name = "Clusters2",
    embedding = "UMAP",
    size = 0.8
    )

plotPDF(p1, name = "UMAP-Clusters2.pdf", ArchRProj = archrproj, embedding = "UMAP", addDOC = FALSE, width = 5, height = 5)

#peak2gene
archrproj <- addPeak2GeneLinks(
    ArchRProj = archrproj,
    reducedDims = "Harmony"
)

p2g <- getPeak2GeneLinks(
    ArchRProj = archrproj,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)

p <- plotPeak2GeneHeatmap(ArchRProj = archrproj, groupBy = "Clusters2")

plotPDF(p, name = "Peak2GeneHeatmap2.pdf", ArchRProj = archrproj, addDOC = FALSE, width = 8, height = 6)

#  differential peak accessibility
diff_peakTest <- getMarkerFeatures(
  ArchRProj = archrproj, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cyc. Prog.",
  bgdGroups = "GluN5"
)

p <- plotMarkers(seMarker = diff_peakTest, name = "Cyc. Prog.", cutOff = "FDR <= 0.43 & abs(Log2FC) >= 1", plotAs = "MA")
p2 <- plotMarkers(seMarker = diff_peakTest, name = "Cyc. Prog.", cutOff = "FDR <= 0.43 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(p, p2, name = "GluN5-vs-CycProg-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = archrproj, addDOC = FALSE)


# TF motif enrichment
motifsUp <- peakAnnoEnrichment(
    seMarker = diff_peakTest,
    ArchRProj = archrproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.2 & Log2FC >= 0.25"
  )

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
up_tfs <- df$TF[1:5]
print("Upregulated TFs in Cyc. Prog. >>>>>>>>>>>>")
print(up_tfs)

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  ggtitle("Upregulated TFs in Cyc. Prog. Cells") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsDo <- peakAnnoEnrichment(
    seMarker = diff_peakTest,
    ArchRProj = archrproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.3 & Log2FC <= -0.25"
  )

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
down_tfs <- df$TF[1:5]
print("Upregulated TFs in GluN5. >>>>>>>>>>>>")
print(down_tfs)

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  ggtitle("Upregulated TFs in GluN5 Cells") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

selected_tfs <- c(up_tfs, down_tfs)

p <- plotEmbedding(
    ArchRProj = archrproj,
    colorBy = "MotifMatrix",
    name = c(paste0("z:", up_tfs), paste0("z:", down_tfs)),
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archrproj)
)

plotPDF(p, name = "TF_UMAP_Plots.pdf", width = 10, height = 8, ArchRProj = archrproj)

plotPDF(ggUp, ggDo, name = "GluN5vsCycProg-Markers-Motifs-Enriched", width = 8, height = 8, ArchRProj = archrproj, addDOC = FALSE)

# footprinting
motifPositions <- getPositions(archrproj)
motifs <- up_tfs[1:3]

archrproj <- addGroupCoverages(ArchRProj = archrproj, groupBy = "Clusters2")

seFoot <- getFootprints(
  ArchRProj = archrproj, 
  positions = motifPositions[motifs], 
  groupBy = "Clusters2"
)

seFootSubset <- getFootprints(
    ArchRProj = archrproj,
    positions = motifPositions[motifs],
    groupBy = "Clusters2",
    useGroups = c("GluN5", "Cyc. Prog.")
)

plotFootprints(
  seFoot = seFootSubset,
  ArchRProj = archrproj, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

# co access
archrproj <- addCoAccessibility(
    ArchRProj = archrproj,
    reducedDims = "Harmony"
)

topGenes <- c("ID4", "OLIG2")

p <- plotBrowserTrack(
    ArchRProj = archrproj, 
    groupBy = "Clusters2", 
    geneSymbol = topGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(archrproj)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = archrproj, 
    addDOC = FALSE, width = 5, height = 5)

genes <- c("ID4", "NFIA", "EGR1", "FOS", "ASCL1", "HES5", "OLIG2", "MEIS2",
          "NHLH1", "SOX21", "SOX10", "PBX1", "NEUROD1", "NEUROG2", "EOMES")

p <- plotBrowserTrack(
    ArchRProj = archrproj, 
    groupBy = "Clusters2", 
    geneSymbol = genes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(archrproj)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-ALL-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = archrproj, 
    addDOC = FALSE, width = 5, height = 5)