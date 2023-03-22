```{r}
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
#library(EnsDb.Mmusculus.v75)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(cowplot)
library(patchwork) #latest version is required!
library(TFBSTools)
library(JASPAR2020)
library(gsfisher)
library(EnhancedVolcano)
library(monocle)
setwd("/rds/projects/c/croftap-stia-atac/CM_multiome/STIA_andATAC/")
options(bitmapType='cairo')
library(CellChat)
#load()
```


```{r}
STIA_RNA_counts <- Read10X_h5("/rds/projects/c/croftap-aia-seq-data/STIA_old/filtered_gene_bc_matrices_h5.h5")
STIA_RNA <- CreateSeuratObject(counts = STIA_RNA_counts, assay = "RNA")
rm(STIA_RNA_counts)
samples_ID <- read.csv(file.path("/rds/projects/c/croftap-aia-seq-data/STIA_old/", "LibraryID.csv"))
STIA_RNA$sampleID=(samples_ID[match(rownames(STIA_RNA@meta.data),samples_ID$Barcode),2])
table(samples_ID[,2])
STIA_RNA$orig.ident <- STIA_RNA$sampleID
table(STIA_RNA@meta.data[["orig.ident"]])
STIA_RNA$sampleID <- NULL
rm(samples_ID)
table(STIA_RNA$orig.ident)
```

```{r}
CON_A <- STIA_RNA[,grepl("CON_A", STIA_RNA$orig.ident, ignore.case=TRUE)]
CON_B <- STIA_RNA[,grepl("CON_B", STIA_RNA$orig.ident, ignore.case=TRUE)]
CON_C <- STIA_RNA[,grepl("CON_C", STIA_RNA$orig.ident, ignore.case=TRUE)]

PEAK_A <- STIA_RNA[,grepl("PEAK_A", STIA_RNA$orig.ident, ignore.case=TRUE)]
PEAK_B <- STIA_RNA[,grepl("PEAK_B", STIA_RNA$orig.ident, ignore.case=TRUE)]
PEAK_C <- STIA_RNA[,grepl("PEAK_C", STIA_RNA$orig.ident, ignore.case=TRUE)]

RESED_A <- STIA_RNA[,grepl("RESED_A", STIA_RNA$orig.ident, ignore.case=TRUE)]
RESED_B <- STIA_RNA[,grepl("RESED_B", STIA_RNA$orig.ident, ignore.case=TRUE)]
RESED_C <- STIA_RNA[,grepl("RESED_C", STIA_RNA$orig.ident, ignore.case=TRUE)]

RESING_A <- STIA_RNA[,grepl("RESING_A", STIA_RNA$orig.ident, ignore.case=TRUE)]
RESING_B <- STIA_RNA[,grepl("RESING_B", STIA_RNA$orig.ident, ignore.case=TRUE)]
RESING_C <- STIA_RNA[,grepl("RESING_C", STIA_RNA$orig.ident, ignore.case=TRUE)]

```


```{r}

mito.features1 <- grep(pattern="^mt-", x=rownames(x=CON_A), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = CON_A, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = CON_A, slot = "counts"))
CON_A[["percent.mito"]] <- percent.mito1
VlnPlot(object = CON_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

```
```{r}
CON_Asafe <- CON_A 
CON_A <- subset(x = CON_A, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.1)
CON_A <- NormalizeData(object = CON_A, verbose = F)
CON_A <- FindVariableFeatures(object = CON_A, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCON_A <- rownames(CON_A)
CON_A <- ScaleData(CON_A, features = all.genesCON_A)
VlnPlot(object = CON_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CON_A <- RunPCA(CON_A)
rm(CON_Asafe)
```



```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=CON_B), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = CON_B, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = CON_B, slot = "counts"))
CON_B[["percent.mito"]] <- percent.mito1
VlnPlot(object = CON_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

```






```{r}
CON_Bsafe <- CON_B 
CON_B <- subset(x = CON_B, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.1)
CON_B <- NormalizeData(object = CON_B, verbose = F)
CON_B <- FindVariableFeatures(object = CON_B, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCON_B <- rownames(CON_B)
CON_B <- ScaleData(CON_B, features = all.genesCON_B)
VlnPlot(object = CON_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CON_B <- RunPCA(CON_B)
rm(CON_Bsafe)
```





```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=CON_C), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = CON_C, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = CON_C, slot = "counts"))
CON_C[["percent.mito"]] <- percent.mito1
VlnPlot(object = CON_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

```

```{r}
CON_Csafe <- CON_C 
CON_C <- subset(x = CON_C, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.1)
CON_C <- NormalizeData(object = CON_C, verbose = F)
CON_C <- FindVariableFeatures(object = CON_C, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCON_C <- rownames(CON_C)
CON_C <- ScaleData(CON_C, features = all.genesCON_C)
VlnPlot(object = CON_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CON_C <- RunPCA(CON_C)
rm(CON_Csafe)

```

```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=PEAK_A), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = PEAK_A, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = PEAK_A, slot = "counts"))
PEAK_A[["percent.mito"]] <- percent.mito1
VlnPlot(object = PEAK_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
```
```{r}
PEAK_Asafe <- PEAK_A 
PEAK_A <- subset(x = PEAK_Asafe, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mito < 0.1)
PEAK_A <- NormalizeData(object = PEAK_A, verbose = F)
PEAK_A <- FindVariableFeatures(object = PEAK_A, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesPEAK_A <- rownames(PEAK_A)
PEAK_A <- ScaleData(PEAK_A, features = all.genesPEAK_A)
VlnPlot(object = PEAK_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
PEAK_A <- RunPCA(PEAK_A)
rm(PEAK_Asafe)


```

```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=PEAK_B), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = PEAK_B, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = PEAK_B, slot = "counts"))
PEAK_B[["percent.mito"]] <- percent.mito1
VlnPlot(object = PEAK_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
```
```{r}
PEAK_Bsafe <- PEAK_B 
PEAK_B <- subset(x = PEAK_B, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mito < 0.1)
PEAK_B <- NormalizeData(object = PEAK_B, verbose = F)
PEAK_B <- FindVariableFeatures(object = PEAK_B, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesPEAK_B <- rownames(PEAK_B)
PEAK_B <- ScaleData(PEAK_B, features = all.genesPEAK_B)
VlnPlot(object = PEAK_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
PEAK_B <- RunPCA(PEAK_B)
rm(PEAK_Bsafe)
```
```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=PEAK_C), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = PEAK_C, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = PEAK_C, slot = "counts"))
PEAK_C[["percent.mito"]] <- percent.mito1
VlnPlot(object = PEAK_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
```
```{r}
PEAK_Csafe <- PEAK_C 
PEAK_C <- subset(x = PEAK_C, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mito < 0.1)
PEAK_C <- NormalizeData(object = PEAK_C, verbose = F)
PEAK_C <- FindVariableFeatures(object = PEAK_C, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesPEAK_C <- rownames(PEAK_C)
PEAK_C <- ScaleData(PEAK_C, features = all.genesPEAK_C)
VlnPlot(object = PEAK_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
PEAK_C <- RunPCA(PEAK_C)
rm(PEAK_Csafe)
```
```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=RESED_A), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = RESED_A, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = RESED_A, slot = "counts"))
RESED_A[["percent.mito"]] <- percent.mito1
VlnPlot(object = RESED_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


 
```

```{r}
RESED_Asafe <- RESED_A 
RESED_A <- subset(x = RESED_A, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.1)
RESED_A <- NormalizeData(object = RESED_A, verbose = F)
RESED_A <- FindVariableFeatures(object = RESED_A, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesRESED_A <- rownames(RESED_A)
RESED_A <- ScaleData(RESED_A, features = all.genesRESED_A)
VlnPlot(object = RESED_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
RESED_A <- RunPCA(RESED_A)
rm(RESED_Asafe)
```

```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=RESED_B), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = RESED_B, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = RESED_B, slot = "counts"))
RESED_B[["percent.mito"]] <- percent.mito1
VlnPlot(object = RESED_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
```
```{r}
RESED_Bsafe <- RESED_B 
RESED_B <- subset(x = RESED_B, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.1)
RESED_B <- NormalizeData(object = RESED_B, verbose = F)
RESED_B <- FindVariableFeatures(object = RESED_B, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesRESED_B <- rownames(RESED_B)
RESED_B <- ScaleData(RESED_B, features = all.genesRESED_B)
VlnPlot(object = RESED_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
RESED_B <- RunPCA(RESED_B)
rm(RESED_Bsafe)
```
```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=RESED_C), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = RESED_C, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = RESED_C, slot = "counts"))
RESED_C[["percent.mito"]] <- percent.mito1
VlnPlot(object = RESED_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
```

```{r}
RESED_Csafe <- RESED_C 
RESED_C <- subset(x = RESED_C, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mito < 0.1)
RESED_C <- NormalizeData(object = RESED_C, verbose = F)
RESED_C <- FindVariableFeatures(object = RESED_C, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesRESED_C <- rownames(RESED_C)
RESED_C <- ScaleData(RESED_C, features = all.genesRESED_C)
VlnPlot(object = RESED_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
RESED_C <- RunPCA(RESED_C)
rm(RESED_Csafe)
```
```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=RESING_A), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = RESING_A, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = RESING_A, slot = "counts"))
RESING_A[["percent.mito"]] <- percent.mito1
VlnPlot(object = RESING_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
```


```{r}
RESING_Asafe <- RESING_A 
RESING_A <- subset(x = RESING_A, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mito < 0.1)
RESING_A <- NormalizeData(object = RESING_A, verbose = F)
RESING_A <- FindVariableFeatures(object = RESING_A, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesRESING_A <- rownames(RESING_A)
RESING_A <- ScaleData(RESING_A, features = all.genesRESING_A)
VlnPlot(object = RESING_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
RESING_A <- RunPCA(RESING_A)
rm(RESING_Asafe)
```
```{r}
mito.features1 <- grep(pattern="^mt-", x=rownames(x=RESING_B), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = RESING_B, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = RESING_B, slot = "counts"))
RESING_B[["percent.mito"]] <- percent.mito1
VlnPlot(object = RESING_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
```


```{r}
RESING_Bsafe <- RESING_B 
RESING_B <- subset(x = RESING_B, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mito < 0.1)
RESING_B <- NormalizeData(object = RESING_B, verbose = F)
RESING_B <- FindVariableFeatures(object = RESING_B, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesRESING_B <- rownames(RESING_B)
RESING_B <- ScaleData(RESING_B, features = all.genesRESING_B)
VlnPlot(object = RESING_B, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
RESING_B <- RunPCA(RESING_B)
rm(RESING_Bsafe)
```
```{r}

mito.features1 <- grep(pattern="^mt-", x=rownames(x=RESING_C), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = RESING_C, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = RESING_C, slot = "counts"))
RESING_C[["percent.mito"]] <- percent.mito1
VlnPlot(object = RESING_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


```

```{r}
RESING_Csafe <- RESING_C 
RESING_C <- subset(x = RESING_C, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mito < 0.1)
RESING_C <- NormalizeData(object = RESING_C, verbose = F)
RESING_C <- FindVariableFeatures(object = RESING_C, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesRESING_C <- rownames(RESING_C)
RESING_C <- ScaleData(RESING_C, features = all.genesRESING_C)
VlnPlot(object = RESING_C, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
RESING_C <- RunPCA(RESING_C)
rm(RESING_Csafe)
```


```{r}

reference.list <- c(CON_A, CON_B, CON_C)
anchors_ctrl <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated_ctrl <- IntegrateData(anchorset = anchors_ctrl, dims = 1:30)

reference.list <- c(PEAK_A, PEAK_B, PEAK_C)
anchors_peak <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated_peak <- IntegrateData(anchorset = anchors_peak, dims = 1:30)

reference.list <- c(RESED_A, RESED_B, RESED_C)
anchors_RESED <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated_RESED <- IntegrateData(anchorset = anchors_RESED, dims = 1:30)

reference.list <- c(RESING_A, RESING_B, RESING_C)
anchors_RESING <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated_RESING <- IntegrateData(anchorset = anchors_RESING, dims = 1:30)

reference.list_all <- c(integrated_ctrl,integrated_peak,integrated_RESING,integrated_RESED)
anchors_all <- FindIntegrationAnchors(object.list = reference.list_all, dims = 1:30)
aggr <- IntegrateData(anchorset = anchors_all, dims = 1:30)

```

```{r}
aggr$all<-'1'
Idents(aggr)<-'all'
VlnPlot(object = aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.00)
```




```{r}
DefaultAssay(object=aggr) <- "integrated"
aggr <- ScaleData(object = aggr, verbose=F)
aggr <- RunPCA(object = aggr, verbose=F)
ElbowPlot(object = aggr)
aggr <- SCTransform(aggr)
DefaultAssay(object=aggr) <- "integrated"
aggr <- FindNeighbors(object = aggr, dims = 1:20)
aggr <- FindClusters(object = aggr, dim.use= 1:20, resolution = 0.3, graph.name = "integrated_snn")
aggr <- RunUMAP(object = aggr, reduction = "pca", dims = 1:20)
aggr <- RunTSNE(object = aggr, reduction = "pca", dims = 1:20)
p5 <- DimPlot(object = aggr, reduction = "umap", pt.size=0.5, label = T, label.size = 5, group.by = "orig.ident")
p6 <- DimPlot(object = aggr, reduction = "tsne", pt.size=0.5, label = T, label.size = 5, group.by = "orig.ident")
plot_grid(p5, p6)

DimPlot(object = aggr, reduction = "umap", pt.size=0.5, label = F, label.size = 5, group.by = "orig.ident")


```

```{r}


for (res in c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7)) {
   aggr <- FindClusters(aggr, resolution = res, algorithm = 3, graph.name = "integrated_snn")
                          
}

library(clustree)
clustree(aggr, assay="integrated")


for (res in c(0.01, 0.05, 0.025, 0.075)) {
   aggr <- FindClusters(aggr, resolution = res, algorithm = 3, graph.name = "integrated_snn")
                          
}

library(clustree)
clustree(aggr, assay="integrated")

for (res in c(0.03, 0.04)) {
   aggr <- FindClusters(aggr, resolution = res, algorithm = 3, graph.name = "integrated_snn")
                          
}

library(clustree)
clustree(aggr, assay="integrated")

DimPlot(aggr, group.by = "integrated_snn_res.0.05")

```


```{r}
Idents(aggr)<-'integrated_snn_res.0.05'
all.markers_RNA <- FindAllMarkers(aggr)


```
```{r}
library(viridis)
library(dplyr)
all.markers_top <- all.markers_RNA %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

all.markers_top <- all.markers_top[-c(66), ]
all.markers_top <- all.markers_top[-c(66), ]
all.markers_top <- all.markers_top[-c(80), ]

DotPlot(aggr, features = all.markers_top$gene) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma")  +RotatedAxis()+coord_flip()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))
ggsave("myplot.pdf")
```

```{r}
#pan fibroblast
DefaultAssay(aggr)<-'RNA'

VlnPlot(aggr, features = c("Cd248", "Cd34", "Cdh11", "Col1a1", "Thy1", "Pdgfra", "Pdpn", "Fap"), pt.size = 0, stack = T)


FeaturePlot(aggr, features = "Fap")
FeaturePlot(aggr, features = "Col1a1")
FeaturePlot(aggr, features = "Thy1")
FeaturePlot(aggr, features = "Pdgfra")
FeaturePlot(aggr, features = "Pdpn")
FeaturePlot(aggr, features = "Cd248")


```

```{r}


 

#linning layer
VlnPlot(aggr, features = c("Cd55", "Clic5", "Col22a1", "Hbegf", "Htra4", "Tspan15", "Prg4"),  pt.size = 0, stack = T)
FeaturePlot(aggr, features = c("Cd55", "Clic5", "Col22a1", "Hbegf", "Htra4", "Tspan15"))
```
```{r}
#osteoblasts
VlnPlot(aggr, features = c("Alpl", "Bglap", "Bglap2", "Omd", "Runx2", "Sp7", "Ostn"),  pt.size = 0, stack = T)


```
```{r}
#chondorcytes
VlnPlot(aggr, features = c("Chad",  "Cilp", "Clu", "Sox9"),  pt.size = 0, stack=T)

```
```{r}
#vascular
VlnPlot(aggr, features = c("Cdh5", "Emcn", "Pecam1"),  pt.size = 0, stack=T)
```

```{r}
#pericytes
VlnPlot(aggr, features=c("Acta2", "Des", "Flt1", "Mcam", "Notch3", "Pdgfrb", "Rgs5"), stack = T,  pt.size = 0)
```
```{r}
VlnPlot(aggr, features=c("Actn3", "Aldoa", "Tnnt3"), stack = T,  pt.size = 0)
```
```{r}
#cell cylel
VlnPlot(aggr, features = c("Cdk1", "Cenpa", "Top2a", "Mki67"), stack = T,  pt.size = 0)
```
```{r}
VlnPlot(aggr, features = c("H2-Aa", "H2-Ab1"),  pt.size = 0)
```
```{r}

aggr$cm_clusters <- aggr@active.ident
current.sample.ids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")
new.sample.ids <- c("fibros", "chondrocytes", "Osteoblasts", "Vasc1", "cycling_fibros", "peri_vasc", "Vasc2", "muscle", "contam")
aggr@meta.data[["cm_clusters"]] <- plyr::mapvalues(x = aggr@meta.data[["cm_clusters"]], from = current.sample.ids, to = new.sample.ids)

DimPlot(aggr, group.by = "cm_clusters", label = T)
```
```{r}
Idents(aggr)<-'cm_clusters'
DotPlot(aggr, features = "Hk2")+ RotatedAxis()+ scale_size(range = c(2, 8))
VlnPlot(aggr, features = "Hk2")+ RotatedAxis()+ scale_size(range = c(2, 8))+coord_flip()

```


```{r}
DimPlot(object = aggr, reduction = "umap", pt.size=0.2, label =F, label.size = 5, group.by = "cm_clusters")
```
```{r}

DefaultAssay(aggr)<-'integrated'
all.markers_top_5 <- all.markers_RNA %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

all.markers_top_5 <- all.markers_top_5[-c(42), ]
Idents(aggr)<-"cm_clusters"

Idents(aggr)<-"cm_clusters"
levels(aggr)<-c(  "fibros" ,"chondrocytes", "Osteoblasts" , "Vasc1" , "Vasc2",  "peri_vasc" ,  "muscle" ,             "contam",                   "cycling_fibros")

DotPlot(aggr, features = all.markers_top_5$gene) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma")  +RotatedAxis()+coord_flip()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))
```
```{r}
DefaultAssay(aggr)<-'integrated'
Idents(aggr)<-"cm_clusters"
All_markers_integrated<-FindAllMarkers(aggr, only.pos = T)

All_markers_integrated_top5 <- All_markers_integrated %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

All_markers_integrated_top5 <- All_markers_integrated_top5[-c(42), ]

DotPlot(aggr, features = All_markers_integrated_top5$gene) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma")  +RotatedAxis()+coord_flip()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))
```



```{r}


Idents(aggr)<-"cm_clusters"
levels(aggr)<-c(  "fibros" ,"Osteoblasts" , "chondrocytes", "peri_vasc" ,  "muscle" ,     "Vasc1" , "Vasc2",           "contam",                   "cycling_fibros")

DotPlot(aggr, features = c("Thy1",  "Cd248", "Cdh11", "Pdgfra",   "Col14a1", "Col1a1", "Pdpn", "Col22a1", "F13a1", "Tspan15", "Alpl", "Runx2", "Runx3", "Anxa8", "Ank", "Chd9",  "Notch3", "Acta2", "Rgs5","Actn3", "Aldoa", "Tnnt3", "Flt1", "Cdh5", "Pecam1"), idents=c("fibros" ,"Osteoblasts" , "chondrocytes", "peri_vasc" ,  "muscle" ,     "Vasc1" , "Vasc2",           "contam",                   "cycling_fibros")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```


```{r}
fibro_only <- aggr[,grepl("fibros", aggr$cm_clusters, ignore.case=TRUE)]
fibro_only <- fibro_only[,!grepl("cycling_fibros", fibro_only$cm_clusters, ignore.case=TRUE)]
DefaultAssay(object=fibro_only) <- "integrated"
fibro_only <- FindVariableFeatures(object = fibro_only, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesfibro_only <- rownames(fibro_only)
fibro_only <- ScaleData(fibro_only, features = all.genesfibro_only)
fibro_only <- RunPCA(fibro_only)
fibro_only <- SCTransform(fibro_only)
DefaultAssay(object=fibro_only) <- "integrated"
fibro_only <- FindNeighbors(object = fibro_only, dims = 1:20)
#fibro_only <- FindClusters(object = fibro_only, dim.use= 1:20, resolution = 0.3, graph.name = "integrated_snn")
fibro_only <- RunUMAP(object = fibro_only, reduction = "pca", dims = 1:20)
#fibro_only <- RunTSNE(object = fibro_only, reduction = "pca", dims = 1:20)
DimPlot(object = fibro_only, reduction = "umap", pt.size=0.2, label = F, label.size = 5, group.by = "cm_clusters")

```
