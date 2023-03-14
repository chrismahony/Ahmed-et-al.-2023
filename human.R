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
#library(monocle)
setwd("/rds/projects/c/croftap-stia-atac/CM_multiome/STIA_andATAC/")
options(bitmapType='cairo')
library(CellChat)
#load()
```

```{r}

celseq_matrix_ru1_molecules.tsv.725583 <- read.delim("/rds/projects/c/croftap-stia-atac/all_amp/celseq_matrix_ru1_molecules.tsv.725583.gz", header = T, row.names=1)
celseq_matrix_ru1_molecules.tsv.725583[is.na(celseq_matrix_ru1_molecules.tsv.725583)] <- 0
counts=celseq_matrix_ru1_molecules.tsv.725583
celseq_meta.tsv.725591 <- read.delim("/rds/projects/c/croftap-stia-atac/all_amp/celseq_meta.tsv.725591.gz", row.names=1, header=T)
metadata=celseq_meta.tsv.725591
amp <- CreateSeuratObject(counts = counts, meta.data = metadata)

rm(celseq_matrix_ru1_molecules.tsv.725583)
rm(celseq_meta.tsv.725591)
amp$all <- 'all'
Idents(amp) <- 'all'
VlnPlot(amp, features = c("percent_mt_molecules", "genes_detected", "reads"), ncol = 3)
DefaultAssay(amp) <- 'RNA'
amp <- subset(amp, subset = genes_detected > 200 & genes_detected < 5000 & percent_mt_molecules < 0.25)
amp <- NormalizeData(amp)
amp <- FindVariableFeatures(amp)
all.genes <- rownames(amp)
amp <- ScaleData(amp, features = all.genes)
amp <- RunPCA(amp, features = VariableFeatures(object = amp))
ElbowPlot(amp)
amp <- FindNeighbors(amp, dims = 1:20)
amp <- FindClusters(amp, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
amp <- FindClusters(amp, resolution = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08))
clustree(amp)
amp <- FindClusters(amp, resolution =0.08)
amp <- RunUMAP(amp, dims = 1:20)
DimPlot(amp, reduction = "umap")
library(harmony)
amp <- RunHarmony(amp, group.by.vars = "orig.ident")
amp <- FindNeighbors(amp, dims = 1:20, reduction = "harmony", verbose = FALSE)
amp <- RunUMAP(amp, dims = 1:20, reduction = "harmony")
DimPlot(amp)
amp$type_cm<-amp$type
Idents(amp)<-'type_cm'
p1<-DimPlot(amp, group.by = "type")
amp<-CellSelector(p1, object = amp, ident = "linning")
amp$type_cm<-amp@active.ident
table(amp$type_cm)
saveRDS(amp, "amp.rds")
VlnPlot(amp, features = c("PRG4", "CLIC5", "TSPAN15", "CD55", "COL22A1", "HBEGF", "HTRA4"))
DimPlot(amp, group.by = "type")
```


```{r}
amp_fibrosonly <- amp[,grepl("Fibroblast|linning", amp$type_cm, ignore.case=TRUE)]
amp_fibrosonly <- FindVariableFeatures(object = amp_fibrosonly, selection.method = "vst", nfeatures = 2000, verbose=F)
amp_fibrosonly_genes <- rownames(amp_fibrosonly)
amp_fibrosonly <- ScaleData(amp_fibrosonly, features = amp_fibrosonly_genes)
amp_fibrosonly <- RunPCA(amp_fibrosonly)
amp_fibrosonly <- FindNeighbors(object = amp_fibrosonly, dims = 1:20)
amp_fibrosonly <- RunUMAP(object = amp_fibrosonly, reduction = "pca", dims = 1:20)
DimPlot(object = amp_fibrosonly, pt.size=0.2, label = F, label.size = 5, group.by = "disease")
```
```{r}
library(harmony)
amp<-RunHarmony(amp, group.by.vars = "disease")
amp_fibrosonly <- FindNeighbors(object = amp_fibrosonly, dims = 1:20, reduction = "harmony")
amp_fibrosonly <- RunUMAP(object = amp_fibrosonly, reduction = "harmony", dims = 1:20)
DimPlot(object = amp_fibrosonly, pt.size=0.2, label = F, label.size = 5, group.by = "disease")
DimPlot(object = amp_fibrosonly, pt.size=0.2, label = F, label.size = 5, group.by = "type_cm")



```
```{r}

Idents(amp)<-'type'
amp<-AddModuleScore(amp, features=list("HK2", "GLS"), name="double_pos")

DotPlot(amp, features = c("HK2", "GLS", "double_pos1"))+RotatedAxis()

DotPlot(amp, features = c("HK2"))+RotatedAxis()
DotPlot(amp, features = c("GLS"))+RotatedAxis()


```



```{r}
fibro_only=amp_fibrosonly

EXPR_Hk2 = GetAssayData(object=fibro_only,assay="RNA",slot="data")["HK2",]
EXPR_Hk2_df=data.frame( positive= EXPR_Hk2 > 0, negative = EXPR_Hk2 == 0)
names(EXPR_Hk2_df)<-paste0( c("positive_","negative_"),"HK2")
fibro_only <- AddMetaData(fibro_only,metadata=EXPR_Hk2_df)


EXPR_Gls = GetAssayData(object=fibro_only,assay="RNA",slot="data")["GLS",]
EXPR_Gls_df=data.frame( positive= EXPR_Gls > 0, negative = EXPR_Gls == 0)
names(EXPR_Gls_df)<-paste0( c("positive_","negative_"),"GLS")
fibro_only <- AddMetaData(fibro_only,metadata=EXPR_Gls_df)


Idents(fibro_only)<-'positive_HK2'
fibro_only$Hk2_Glspos <- paste(Idents(fibro_only), fibro_only$positive_GLS, sep = "_")
Idents(fibro_only)<-'Hk2_Glspos'
table(fibro_only$Hk2_Glspos)

fibro_only$Hk2_Glspos_renamed <- fibro_only$Hk2_Glspos
current.sample.ids <- c("FALSE_FALSE" , "FALSE_TRUE" , "TRUE_FALSE" ,  "TRUE_TRUE" )
new.sample.ids <- c("double_neg", "Gls_pos", "Hk2_pos", "double_pos")
fibro_only@meta.data[["Hk2_Glspos_renamed"]] <- plyr::mapvalues(x = fibro_only@meta.data[["Hk2_Glspos_renamed"]], from = current.sample.ids, to = new.sample.ids)
Idents(fibro_only) <- "Hk2_Glspos_renamed"
table(fibro_only$Hk2_Glspos_renamed)

DimPlot(fibro_only, group.by = "Hk2_Glspos_renamed")
VlnPlot(fibro_only, features = c("GLS", "HK2"), stack = T)

```

```{r}
Idents(fibro_only)<-'Hk2_Glspos_renamed'
fibro_only$Hk2_Glspos_renamed_disease <- paste(Idents(fibro_only), fibro_only$disease, sep = "_")
Idents(fibro_only)<-'Hk2_Glspos_renamed_disease'
table(fibro_only$Hk2_Glspos_renamed_disease)

```

```{r}
Idents(fibro_only)<-'Hk2_Glspos_renamed'

DotPlot(fibro_only, features = c("HK2", "GLS", "HAS1", "COL12A1","COL3A1", "CTHRC1","CCL2", "IL6", "MMP3", "TSC22D3", "MYOC", "ANGPTL7", "IGFBP6", "ITM2B"), idents = c("double_pos", "double_neg"))+RotatedAxis()+coord_flip()

levels(fibro_only)<-c("double_pos", "Gls_pos",  "Hk2_pos" , "double_neg")

DotPlot(fibro_only, features = c("HK2", "GLS", "HAS1", "COL12A1","COL3A1", "CTHRC1","CCL2", "IL6", "MMP3", "TSC22D3", "MYOC", "ANGPTL7", "IGFBP6", "ITM2B"), )+RotatedAxis()+coord_flip()



```
