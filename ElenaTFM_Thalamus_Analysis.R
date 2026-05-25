# CHICKEN THALAMUS SINGLE-NUCLEI ATLAS QUALITY CONTROL, PREPROCESSING, CLUSTERING AND ANNOTATION
# Elena del Rocío Iglesias - Neuroscience Master Thesis (EHU), 2025-26


# Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(dplyr)


# Import data
gg.E14.1.counts <- ReadMtx(mtx = "Sample11.scRNA.filtered.matrix.mtx.gz", cells = "Sample11.scRNA.filtered.barcodes.tsv.gz",
                                features = "Sample11.scRNA.filtered.features.tsv.gz")

# Build Seurat object
gg.E14.1 <- CreateSeuratObject(counts=gg.E14.1.counts, project = "FGM-E14.1")

# Add metadata
gg.E14.1 <- AddMetaData(gg.E14.1, "E14.1", col.name = "Sample_Name")
gg.E14.1 <- AddMetaData(gg.E14.1, "Chick", col.name = "Species")
gg.E14.1 <- AddMetaData(gg.E14.1, "E14", col.name = "Stage")
gg.E14.1 <- AddMetaData(gg.E14.1, "Embryonic immature neuron NUCLEI", col.name = "Cell_type")
gg.E14.1 <- AddMetaData(gg.E14.1, "Brain", col.name = "Organ")
gg.E14.1 <- AddMetaData(gg.E14.1, "Thalamus", col.name = "Area")
gg.E14.1 <- AddMetaData(gg.E14.1, "Fluent/Ilumina", col.name = "Protocol")
gg.E14.1 <- AddMetaData(gg.E14.1, "T10", col.name = "Kit_version")

# %MT calculation
gg.E14.1[["Percent_mt"]] <- PercentageFeatureSet(gg.E14.1, features = c("COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB", "ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6", "ENSGALG00010000002", "ENSGALG00010000003", "ENSGALG00010000004", "ENSGALG00010000005", "ENSGALG00010000006", "ENSGALG00010000008", "ENSGALG00010000009", "ENSGALG000100000010", "ENSGALG000100000012", "ENSGALG000100000013", "ENSGALG000100000014", "ENSGALG000100000015", "ENSGALG000100000016", "ENSGALG000100000018", "ENSGALG000100000019", "ENSGALG000100000021", "ENSGALG000100000025", "ENSGALG000100000027", 
                                                                          "ENSGALG000100000030", "ENSGALG000100000031", "ENSGALG000100000032", "ENSGALG000100000035", "ENSGALG000100000036", "ENSGALG000100000038"))

# Quality control

# Visualization and QC filtering
VlnPlot(gg.E14.1, features = c("nFeature_RNA", "nCount_RNA", "Percent_mt"), group.by = "orig.ident")
FeatureScatter(gg.E14.1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") 
summary(gg.E14.1@meta.data)
high.nF <- quantile(gg.E14.1$nFeature_RNA, 0.98)
high.nC <- quantile(gg.E14.1$nCount_RNA, 0.98)
high.MT <- quantile(gg.E14.1$Percent_mt, 0.97)
gg.E14.1 <- subset(gg.E14.1, subset = nFeature_RNA > 600 & nFeature_RNA < high.nF & nCount_RNA < high.nC & Percent_mt < high.MT)

# Cell cycle scoring
fgm.s.genes <- c("GINS1", "RFC1", "WXO1", "ENSGALG00000011747", "H2AFJ", "ENSGALG00000042491", "HIST1H110")
fgm.g2m.genes <- c("BUB1B", "NCAPG", "NCAPH", "CDC25B", "SMC2", "SPC25")
s.genes <- toupper(c(cc.genes$s.genes))
g2m.genes <- toupper(c(cc.genes$g2m.genes))

options(future.globals.maxSize = 6*1024^3)
gg.E14.1 <- SCTransform(gg.E14.1)
gg.E14.1 <- CellCycleScoring(gg.E14.1, s.features = c(s.genes, fgm.s.genes), g2m.features = c(g2m.genes, fgm.g2m.genes), set.ident = TRUE, nbin = 10)

# Preprocessing
# Normalization, scaling, out-regression...
gg.E14.1 <- SCTransform(gg.E14.1, vars.to.regress = c("Percent_mt", "S.Score", "G2M.Score"))

# Dimensional reduction
gg.E14.1 <- RunPCA(gg.E14.1)
ElbowPlot(gg.E14.1, ndims = 50) # 30 PCs
gg.E14.1 <- FindNeighbors(gg.E14.1, dims = 1:30)
gg.E14.1 <- FindClusters(gg.E14.1, resolution = c(0.5, 1, 1.5, 2))
set.seed(123)
gg.E14.1 <- RunUMAP(gg.E14.1, dims = 1:30)

# Clustering
DimPlot(gg.E14.1, group.by = "SCT_snn_res.0.5", label = T)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.1", label = TRUE)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.1.5", label = T)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.2", label = T)

# Visualization - Which is the best resolution for our data?
DimPlot(gg.E14.1, reduction="umap", group.by = "Phase", pt.size = 0.5)
FeaturePlot(gg.E14.1, c("S.Score", "G2M.Score"), col = c("white","#8D8680"), pt.size = 0.7)

# Save processed object
#saveRDS(gg.E14.1, "elenaTFM_processed_object_E14_1.rds")
#gg.E14.1 <- readRDS("elenaTFM_processed_object_E14_1.rds")


# Annotations

# Markers
gg.E14.1 <- PrepSCTFindMarkers(gg.E14.1, assay = "SCT", verbose = TRUE)
options(future.globals.maxSize = 7 * 1024^3)
gg.E14.1.markers <- FindAllMarkers(object = gg.E14.1,
                                            only.pos = TRUE,
                                            verbose = FALSE,
                                            group.by = "SCT_snn_res.0.5",
                                            min.pct = 0.25,
                                            logfc.threshold = 0.5)

top <- gg.E14.1.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top <- top %>%
  filter(!grepl("^ENSGAL", gene))
DotPlot(gg.E14.1, features = unique(top$gene), group.by = "SCT_snn_res.0.5") + coord_flip() +
  theme(axis.text.y = element_text(size = 8))
#DoHeatmap(gg.E14.1.markers, features = top3$gene, group.by = "SCT_snn_res.1")

# Known gene markers
gene.list <- c("VWF", "FLT1",
               "TFEC",
               "TOP2A",
               "AQP4",
               "OLIG2", "PDGFRA", "PLP1",
               "SOX2", "RFX4", "SNAP25",
               "GAD2", "SLC17A6",
               "LMO3", "PAX6",
               "LHX9", "GBX2", "PROX1",
               "ISL1", "MEIS2", "GATA3", "SOX14",
               "AVP")
for(gene in gene.list){
  print(FeaturePlot(gg.E14.1, gene))
}

# Filter non-neural cells
gg.E14.1 <- subset(gg.E14.1, SCT_snn_res.1 != "15" & SCT_snn_res.1 != "16" & SCT_snn_res.1 != "23")

gg.E14.1 <- SCTransform(gg.E14.1, vars.to.regress = c("Percent_mt", "S.Score", "G2M.Score"))
gg.E14.1 <- RunPCA(gg.E14.1)
ElbowPlot(gg.E14.1, ndims = 50) # 30 PCs
gg.E14.1 <- FindNeighbors(gg.E14.1, dims = 1:30)
gg.E14.1 <- FindClusters(gg.E14.1, resolution = c(0.5, 1, 1.5, 2))
set.seed(123)
gg.E14.1 <- RunUMAP(gg.E14.1, dims = 1:30)


#saveRDS(gg.E14.1, "ElenaTFM_gg-E14_1.rds")
#gg.E14.1 <- readRDS("ElenaTFM_gg-E14_1.rds")

DimPlot(gg.E14.1, group.by = "SCT_snn_res.0.5", label = T)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.1", label = TRUE)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.1.5", label = T)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.2", label = T)

gene.list <- c("VWF", "FLT1", # Endothelial
               "TFEC", # Microglia
               "TOP2A", "MKI67", # Progenitors
               "AQP4", # Astrocytes
               "OLIG2", "PDGFRA", "PLP1", # Oligos
               
               "SOX2", "RFX4", "SNAP25", # Immatures vs matures (neurons)
               "GAD2", "SLC17A6", # GABAergic vs glutamatergic
               "LMO3", "PAX6", # pretectum
               "LHX9", "GBX2", "PROX1", # thalamus
               "ISL1", "MEIS2", "GATA3", "SOX14", # prethalamus
               "AVP") # hypothalamus
for(gene in gene.list){
  print(FeaturePlot(gg.E14.1, gene))}


# Filter glia
gg.E14.1 <- subset(gg.E14.1, SCT_snn_res.0.5 != "10" & SCT_snn_res.0.5 != "5" & SCT_snn_res.0.5 != "17" 
                   & SCT_snn_res.0.5 != "0" & SCT_snn_res.0.5 != "13")

gg.E14.1 <- SCTransform(gg.E14.1, vars.to.regress = c("Percent_mt", "S.Score", "G2M.Score"))
gg.E14.1 <- RunPCA(gg.E14.1)
ElbowPlot(gg.E14.1, ndims = 50) # 40 PCs
gg.E14.1 <- FindNeighbors(gg.E14.1, dims = 1:40)
gg.E14.1 <- FindClusters(gg.E14.1, resolution = c(0.5, 1, 1.5, 2))
set.seed(123)
gg.E14.1 <- RunUMAP(gg.E14.1, dims = 1:40)


#saveRDS(gg.E14.1, "ElenaTFM_gg-E14_1_neurons.rds")
gg.E14.1 <- readRDS("ElenaTFM_gg-E14_1_neurons.rds")

# Clusters
DimPlot(gg.E14.1, group.by = "SCT_snn_res.0.5", label = T)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.1", label = TRUE)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.1.5", label = T)
DimPlot(gg.E14.1, group.by = "SCT_snn_res.2", label = T)

gene.list <- c("VWF", "FLT1", # Endothelial
               "TFEC", # Microglia
               "TOP2A", "MKI67", # Progenitors
               "AQP4", # Astrocytes
               "OLIG2", "PDGFRA", "PLP1", # Oligos
               
               "SOX2", "RFX4", "SNAP25", # Immatures vs matures (neurons)
               "GAD2", "SLC17A6", # GABAergic vs glutamatergic
               "LMO3", "PAX6", # pretectum
               "LHX9", "GBX2", "PROX1", # thalamus
               "ISL1", "MEIS2", "GATA3", "SOX14", # prethalamus
               "AVP") # hypothalamus

for(gene in gene.list){
  print(FeaturePlot(gg.E14.1, gene))}

# Markers
options(future.globals.maxSize = 8*1024^3)
gg.E14.markers <- FindAllMarkers(object = gg.E14.1,
                                          only.pos = TRUE,
                                          verbose = FALSE,
                                          group.by = "SCT_snn_res.0.5",
                                          min.pct = 0.25,
                                          logfc.threshold = 0.5)
top <- gg.E14.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top <- top %>%
  filter(!grepl("^ENSGAL", gene))
DotPlot(gg.E14.1, features = unique(top$gene), group.by = "SCT_snn_res.0.5") + coord_flip() +
  theme(axis.text.y = element_text(size = 8))


# Rename Idents
Idents(gg.E14.1) <- "SCT_snn_res.1" # Resolution you are working with
gg.E14.1 <- RenameIdents(gg.E14.1,
                              c("0" = "Low quality cells",
                                "1" = "Thalamus",
                                "2" = "Prethalamus",
                                "3" = "Low quality cells",
                                "4" = "Low quality cells",
                                "5" = "Low quality cells",
                                "6" = "Thalamus",
                                "7" = "Pretectum",
                                "8" = "Prethalamus",
                                "9" = "Prethalamus",
                                "10" = "Prethalamus",
                                "11" = "Prethalamus",
                                "12" = "Astrocytes",
                                "13" = "Thalamus",
                                "14" = "Astrocytes",
                                "15" = "Prethalamus",
                                "16" = "Hypothalamus",
                                "17" = "Prethalamus",
                                "18" = "Thalamus",
                                "19" = "Pretectum",
                                "20" = "Prethalamus",
                                "21" = "Prethalamus",
                                "22" = "Low quality cells",
                                "23" = "Astrocytes",
                                "24" = "Prethalamus",
                                "25" = "Pretectum",
                                "26" = "Prethalamus")) # Match cluster number with the new cluster name


# Definitive plots
DimPlot(gg.E14.1, cols = c("Low quality cells" = "#E0E0E0",
                           "Astrocytes" = "#9DDFE8",
                           "Thalamus" = "#FF80FF",
                           "Prethalamus" = "#FFC5D9",
                           "Pretectum" = "#FFE0B3",
                           "Hypothalamus" = "#B1E0B7"), pt.size = 0.7)
ggsave("E14_Elena.png", get_last_plot())

FeaturePlot(gg.E14.1, "LHX9", pt.size = 0.7) # Thalamus
ggsave("LHX9.png", get_last_plot())
FeaturePlot(gg.E14.1, "PROX1", pt.size = 0.7) # Thalamus
ggsave("PROX1.png", get_last_plot())
FeaturePlot(gg.E14.1, "ISL1", pt.size = 0.7) # Prethalamus
ggsave("ISL1.png", get_last_plot())
FeaturePlot(gg.E14.1, "MEIS2", pt.size = 0.7) # Prethalamus
ggsave("MEIS2.png", get_last_plot())
FeaturePlot(gg.E14.1, "GATA3", pt.size = 0.7) # Prethalamus
ggsave("GATA3.png", get_last_plot())
FeaturePlot(gg.E14.1, "LMO3", pt.size = 0.7) # Pretectum
ggsave("LMO3.png", get_last_plot())
FeaturePlot(gg.E14.1, "AQP4", pt.size = 0.7) # Astrocytes
ggsave("AQP4.png", get_last_plot())
FeaturePlot(gg.E14.1, "OTP", pt.size = 0.7) # Hypothalamus
ggsave("OTP.png", get_last_plot())

VlnPlot(gg.E14.1, features = c("nFeature_RNA", "nCount_RNA", "Percent_mt"), group.by = "SCT_snn_res.0.5")
ggsave("low_quality.png", get_last_plot(), width = 16)
