########## Loading packages ##########

library(Seurat)
library(SoupX)
library(DoubletFinder)
library(sctransform)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(VISION)
library(parallel)
library(pheatmap)

########## Removal of ambient RNA ##########

toc = Seurat::Read10X(data.dir = "Y_1/filtered_feature_bc_matrix/")  # data.dir = "Y_2/filtered_feature_bc_matrix/" for sample 2
tod = Seurat::Read10X(data.dir = "Y_1/raw_feature_bc_matrix/")  # data.dir = "Y_2/raw_feature_bc_matrix/" for sample 2
sn = SoupChannel(tod, toc)
srat <- CreateSeuratObject(counts = toc)
srat <- SCTransform(srat, verbose = F)
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, verbose = T)
meta <- srat@meta.data
umap <- srat@reductions$umap@cell.embeddings
sn <- setClusters(sn, setNames(meta$seurat_clusters, rownames(meta)))
sn <- setDR(sn, umap)
sn = autoEstCont(sn)
out = adjustCounts(sn, roundToInt = T)
srat = CreateSeuratObject(out)

########## Quality control ##########

srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
srat <- subset(srat, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 1)

########## Doublet detection ##########

srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
srat <- FindNeighbors(srat, dims = 1:20)
srat <- FindClusters(srat, resolution = 0.5)
srat <- RunUMAP(srat, dims = 1:20)
sweep.res.list_kidney <- paramSweep(srat, PCs = 1:20, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
pK <- bcmvn_kidney %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))
annotations <- srat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.056*nrow(srat@meta.data))  # Assuming 5.6% doublet formation rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
srat_df <- doubletFinder(srat, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
srat_df2 <- doubletFinder(srat_df, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_328", sct = FALSE)  # reuse.pANN = "pANN_0.25_0.01_425" for sample2 data
srat_df2 = subset(srat_df2, DF.classifications_0.25_0.01_300 == "Singlet")  # DF.classifications_0.25_0.01_390 == "Singlet" for sample2 data

saveRDS(srat_df2, file = "sample1.rds")  # file = "sample2.rds" for sample 2

########## Data integration ##########

srat1 <- readRDS("sample1.rds")
srat2 <- readRDS("sample2.rds")

srat <- merge(x = srat1, y = srat2, project = "CNI")
srat@meta.data$stim <- c(rep("Veh", ncol(srat1)), rep("CsA", ncol(srat2)))
srat$stim <- factor(srat$stim, levels = c("Veh", "CsA"))
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat)
srat <- ScaleData(srat)
srat <- RunPCA(srat)
srat <- RunUMAP(srat, dims = 1:20)
srat <- IntegrateLayers(srat, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

########## Clustering ##########

srat <- RunUMAP(srat, dims = 1:20, reduction = "integrated.cca")
srat <- FindNeighbors(srat, reduction = "integrated.cca", dims = 1:20)
srat <- FindClusters(srat, resolution = 0.2)
Idents(srat) = srat$RNA_snn_res.0.2

########## Annotation of clusters ##########

srat <- JoinLayers(srat)
srat.markers <- FindAllMarkers(srat, only.pos = TRUE)
srat.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(srat, ident.1 = 0)
cluster1.markers <- FindMarkers(srat, ident.1 = 1)
cluster2.markers <- FindMarkers(srat, ident.1 = 2)
cluster3.markers <- FindMarkers(srat, ident.1 = 3)
cluster4.markers <- FindMarkers(srat, ident.1 = 4)
cluster5.markers <- FindMarkers(srat, ident.1 = 5)
cluster6.markers <- FindMarkers(srat, ident.1 = 6)
cluster7.markers <- FindMarkers(srat, ident.1 = 7)
cluster8.markers <- FindMarkers(srat, ident.1 = 8)
cluster9.markers <- FindMarkers(srat, ident.1 = 9)
cluster10.markers <- FindMarkers(srat, ident.1 = 10)
cluster11.markers <- FindMarkers(srat, ident.1 = 11)
cluster12.markers <- FindMarkers(srat, ident.1 = 12)
cluster13.markers <- FindMarkers(srat, ident.1 = 13)
cluster14.markers <- FindMarkers(srat, ident.1 = 14)
cluster15.markers <- FindMarkers(srat, ident.1 = 15)
cluster16.markers <- FindMarkers(srat, ident.1 = 16)
cluster17.markers <- FindMarkers(srat, ident.1 = 17)

write.csv(cluster0.markers, "cluster0markers.csv")
write.csv(cluster1.markers, "cluster1markers.csv")
write.csv(cluster2.markers, "cluster2markers.csv")
write.csv(cluster3.markers, "cluster3markers.csv")
write.csv(cluster4.markers, "cluster4markers.csv")
write.csv(cluster5.markers, "cluster5markers.csv")
write.csv(cluster6.markers, "cluster6markers.csv")
write.csv(cluster7.markers, "cluster7markers.csv")
write.csv(cluster8.markers, "cluster8markers.csv")
write.csv(cluster9.markers, "cluster9markers.csv")
write.csv(cluster10.markers, "cluster10markers.csv")
write.csv(cluster11.markers, "cluster11markers.csv")
write.csv(cluster12.markers, "cluster12markers.csv")
write.csv(cluster13.markers, "cluster13markers.csv")
write.csv(cluster14.markers, "cluster14markers.csv")
write.csv(cluster15.markers, "cluster15markers.csv")
write.csv(cluster16.markers, "cluster16markers.csv")
write.csv(cluster17.markers, "cluster17markers.csv")

new.cluster.ids <- c("PCT", "Endo", "TAL", "PST1", "KIM-1+ PT", "DCT", "Fibro", "Uro", "PST2", "CD PC", "Immune", "DTL", "ATL", "CD ICA", "CD ICB", "Podo", "PEC", "Doublets")
names(new.cluster.ids) <- levels(srat)
srat <- RenameIdents(srat, new.cluster.ids)
DimPlot(srat, reduction = "umap", label = TRUE) + NoLegend()
levels(srat) <- c("PCT", "PST1", "PST2", "KIM-1+ PT", "DTL", "ATL", "TAL", "DCT", "CD PC", "CD ICA", "CD ICB", "PEC", "Podo", "Endo", "Fibro", "Uro", "Immune", "Doublets")

features <- c("Slc5a2", "Slc22a6", "Slc22a13", "Havcr1", "Cryab", "Slc4a11", "Slc7a6", "Slc12a1", "Slc12a3", "Aqp2", "Slc26a7", "Slc26a4", "Ncam1", "Nphs1", "Flt1", "Col1a2", "Upk1b", "Ptprc")

levels(srat) <- rev(levels(srat))  # Level change for data visualization purpose only
DotPlot(srat, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat) <- rev(levels(srat))

########## Removal of residual doublets ##########

srat <- subset(srat, idents = "Doublets", invert = TRUE)

########## Data visualization ##########

levels(srat) <- c("PCT", "Endo", "TAL", "PST1", "KIM-1+ PT", "DCT", "Fibro", "Uro", "PST2", "CD PC", "Immune", "DTL", "ATL", "CD ICA", "CD ICB", "Podo", "PEC")  # Level change for data visualization (coloring) purpose only
DimPlot(srat, reduction = "umap", label = TRUE,  pt.size = .1) + NoLegend()
DimPlot(srat, reduction = "umap", label = TRUE,  pt.size = .1, split.by = "stim") + NoLegend()
levels(srat) <- c("PCT", "PST1", "PST2", "KIM-1+ PT", "DTL", "ATL", "TAL", "DCT", "CD PC", "CD ICA", "CD ICB", "PEC", "Podo", "Endo", "Fibro", "Uro", "Immune")

levels(srat) <- rev(levels(srat))  # Level change for data visualization purpose only
DotPlot(srat, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(srat) <- rev(levels(srat))

srat[["clusters"]] <- Idents(srat)
table(srat@meta.data$RNA_snn_res.0.2, srat@meta.data$stim)

FeaturePlot(srat, features = "Havcr1", max.cutoff = 3, cols = c("gray", "red"), split.by = "stim", label = TRUE)
FeaturePlot(srat, features = "Vcam1", max.cutoff = 3, cols = c("gray", "red"), split.by = "stim", label = TRUE)

saveRDS(srat, file = "integrated.rds")

########## Gene set enrichment analysis ##########

srat <- readRDS("integrated.rds")
options(mc.cores=8)

signatures <- "mh.all.v2023.2.Mm.symbols.gmt"
vision.obj <- Vision(srat[["RNA"]]$counts, signatures = signatures, pool = F)
vision.obj <- analyze(vision.obj)
sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)
srat_hallmark <- srat
srat_hallmark[['GSEA']] <- CreateAssayObject(counts = t(sigScores))
DefaultAssay(srat_hallmark) <- "GSEA"

FeaturePlot(srat_hallmark, "HALLMARK-OXIDATIVE-PHOSPHORYLATION", cols =c("#E4F0F0","#F33D0C"), max.cutoff = "q99",min.cutoff = "q1", split.by = "stim")
VlnPlot(srat_hallmark, "HALLMARK-OXIDATIVE-PHOSPHORYLATION", split.plot = TRUE, split.by = "stim", pt.size = 0.3)
FeaturePlot(srat_hallmark, "HALLMARK-GLYCOLYSIS", cols =c("#E4F0F0","#F33D0C"), max.cutoff = "q99",min.cutoff = "q1", split.by = "stim")
VlnPlot(srat_hallmark, "HALLMARK-GLYCOLYSIS", split.plot = TRUE, split.by = "stim", pt.size = 0.3)
FeaturePlot(srat_hallmark, "HALLMARK-TGF-BETA-SIGNALING", cols =c("#E4F0F0","#F33D0C"), max.cutoff = "q99",min.cutoff = "q1", split.by = "stim")
VlnPlot(srat_hallmark, "HALLMARK-TGF-BETA-SIGNALING", split.plot = TRUE, split.by = "stim", pt.size = 0.3)

signatures <- "m5.all.v2023.2.Mm.symbols.gmt"
vision.obj <- Vision(srat[["RNA"]]$counts, signatures = signatures, pool = F)
vision.obj <- analyze(vision.obj)
sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)
srat_go <- srat
srat_go[['GSEA']] <- CreateAssayObject(counts = t(sigScores))
DefaultAssay(srat_go) <- "GSEA"
FeaturePlot(srat_go, "GOBP-CELLULAR-SENESCENCE", cols =c("#E4F0F0","#F33D0C"), max.cutoff = "q99",min.cutoff = "q1", split.by = "stim")
VlnPlot(srat_go, "GOBP-CELLULAR-SENESCENCE", split.plot = TRUE, split.by = "stim", pt.size = 0.3)

########## Differentially expressed genes (DEG) between proximal tubule cells in CNI nephrotoxicity mouse model vs proximal tubule cells in its control ##########

srat <- readRDS("integrated.rds")
srat$celltype.stim <- paste(srat$RNA_snn_res.0.2, srat$stim, sep = "_")
Idents(srat) <- "celltype.stim"

levels(srat) <- c("0_Veh", "1_Veh", "2_Veh", "3_Veh", "4_Veh", "5_Veh", "6_Veh", "7_Veh", "8_Veh", "9_Veh", "10_Veh", "11_Veh", "12_Veh", "13_Veh", "14_Veh", "15_Veh", "16_Veh", "0_CsA", "1_CsA", "2_CsA", "3_CsA", "4_CsA", "5_CsA", "6_CsA", "7_CsA", "8_CsA", "9_CsA", "10_CsA", "11_CsA", "12_CsA", "13_CsA", "14_CsA", "15_CsA", "16_CsA")
new.cluster.ids <- c("PT_Veh", "1_Veh", "2_Veh", "PT_Veh", "PT_Veh", "5_Veh", "6_Veh", "7_Veh", "PT_Veh", "9_Veh", "10_Veh", "11_Veh", "12_Veh", "13_Veh", "14_Veh", "15_Veh", "16_Veh", "PT_CsA", "1_CsA", "2_CsA", "PT_CsA", "PT_CsA", "5_CsA", "6_CsA", "7_CsA", "PT_CsA", "9_CsA", "10_CsA", "11_CsA", "12_CsA", "13_CsA", "14_CsA", "15_CsA", "16_CsA")
names(new.cluster.ids) <- levels(srat)
srat <- RenameIdents(srat, new.cluster.ids)
srat[["clusters.PT"]] <- Idents(srat)

clusterPT.markers <- FindMarkers(srat, ident.1 = "PT_CsA", ident.2 = "PT_Veh")
write.csv(clusterPT.markers, "clusterPTmarkers.csv")

########## Differentially expressed genes (DEG) between KIM-1–positive injured proximal tubule cells in CNI nephrotoxicity mouse model vs other proximal tubule cells in CNI nephrotoxicity mouse model ##########

srat <- readRDS("integrated.rds")
srat$celltype.stim <- paste(srat$RNA_snn_res.0.2, srat$stim, sep = "_")
Idents(srat) <- "celltype.stim"

levels(srat) <- c("0_Veh", "1_Veh", "2_Veh", "3_Veh", "4_Veh", "5_Veh", "6_Veh", "7_Veh", "8_Veh", "9_Veh", "10_Veh", "11_Veh", "12_Veh", "13_Veh", "14_Veh", "15_Veh", "16_Veh", "0_CsA", "1_CsA", "2_CsA", "3_CsA", "4_CsA", "5_CsA", "6_CsA", "7_CsA", "8_CsA", "9_CsA", "10_CsA", "11_CsA", "12_CsA", "13_CsA", "14_CsA", "15_CsA", "16_CsA")
new.cluster.ids <- c("0_Veh", "1_Veh", "2_Veh", "3_Veh", "4_Veh", "5_Veh", "6_Veh", "7_Veh", "8_Veh", "9_Veh", "10_Veh", "11_Veh", "12_Veh", "13_Veh", "14_Veh", "15_Veh", "16_Veh", "PTothers_CsA", "1_CsA", "2_CsA", "PTothers_CsA", "4_CsA", "5_CsA", "6_CsA", "7_CsA", "PTothers_CsA", "9_CsA", "10_CsA", "11_CsA", "12_CsA", "13_CsA", "14_CsA", "15_CsA", "16_CsA")

names(new.cluster.ids) <- levels(srat)
srat <- RenameIdents(srat, new.cluster.ids)
srat[["clusters.PTrest"]] <- Idents(srat)

clusterPTothers.markers <- FindMarkers(srat, ident.1 = "4_CsA", ident.2 = "PTothers_CsA")
write.csv(clusterPTothers.markers, "clusterPTothersmarkers.csv")