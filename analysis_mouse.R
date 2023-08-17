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
DimPlot(aggr, group.by = "cm_clusters", label=T)
ncol(aggr)
```
```{r}


aggr<-AddModuleScore(aggr, features=list("Hk2", "Gls"), name="double_pos")

DotPlot(aggr, features = c("Hk2", "Gls", "double_pos1"))+RotatedAxis()

```



```{r}
load("~/croftap-stia-atac-path/CM_multiome/STIA_andATAC/Monica_Gaum_analysis/STIA_analysis/2018_analysis.RData")
rm(list=ls()[! ls() %in% c("fibro_only")])


DimPlot(object = fibro_only, reduction = "umap", pt.size=0.2, label = F, label.size = 5, group.by = "orig.ident")

```
```{r}
EXPR_Hk2 = GetAssayData(object=fibro_only,assay="RNA",slot="data")["Hk2",]
EXPR_Hk2_df=data.frame( positive= EXPR_Hk2 > 0, negative = EXPR_Hk2 == 0)
names(EXPR_Hk2_df)<-paste0( c("positive_","negative_"),"Hk2")
fibro_only <- AddMetaData(fibro_only,metadata=EXPR_Hk2_df)


EXPR_Gls = GetAssayData(object=fibro_only,assay="RNA",slot="data")["Gls",]
EXPR_Gls_df=data.frame( positive= EXPR_Gls > 0, negative = EXPR_Gls == 0)
names(EXPR_Gls_df)<-paste0( c("positive_","negative_"),"Gls")
fibro_only <- AddMetaData(fibro_only,metadata=EXPR_Gls_df)


Idents(fibro_only)<-'positive_Hk2'
fibro_only$Hk2_Glspos <- paste(Idents(fibro_only), fibro_only$positive_Gls, sep = "_")
Idents(fibro_only)<-'Hk2_Glspos'
table(fibro_only$Hk2_Glspos)


fibro_only$Hk2_Glspos_renamed <- fibro_only$Hk2_Glspos
current.sample.ids <- c("FALSE_FALSE" , "FALSE_TRUE" , "TRUE_FALSE" ,  "TRUE_TRUE" )
new.sample.ids <- c("double_neg", "Gls_pos", "Hk2_pos", "double_pos")
fibro_only@meta.data[["Hk2_Glspos_renamed"]] <- plyr::mapvalues(x = fibro_only@meta.data[["Hk2_Glspos_renamed"]], from = current.sample.ids, to = new.sample.ids)
Idents(fibro_only) <- "Hk2_Glspos_renamed"
table(fibro_only$Hk2_Glspos_renamed)

DimPlot(fibro_only, group.by = "Hk2_Glspos_renamed")
DimPlot(fibro_only, group.by = "orig.ident")

```
```{r}
Idents(fibro_only) <- "Hk2_Glspos_renamed"
VlnPlot(fibro_only, features = c("Hk2", "Gls"), stack=T)


```
```{r}

table(fibro_only$Hk2_Glspos_renamed)


```
```{r}


DotPlot(fibro_only, features = c("Hk2", "Gls"), idents = c("double_pos", "double_neg"))+RotatedAxis()+coord_flip()

```


```{r}

DefaultAssay(fibro_only)<-"RNA"

markers_double<-FindAllMarkers(fibro_only, only.pos = T)

library(gsfisher)

expressed_genes<-rownames(fibro_only)
annotation_gs <- fetchAnnotation(species="mm", ensembl_version=NULL, ensembl_host=NULL)

index <- match(markers_double$gene, annotation_gs$gene_name)
markers_double$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- expressed_genes
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]

seurat_obj.res <- markers_double
seurat_obj <- fibro_only
seurat_obj.res <- seurat_obj.res[!is.na(seurat_obj.res$ensembl),]
ensemblUni <- na.omit(ensemblUni)



go.results <- runGO.all(results=seurat_obj.res,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="p_val_adj", p_threshold=0.05,
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=3, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)

```
```{r}
markers_double %>%
    group_by(cluster) %>%
    top_n(n = 4, wt = avg_log2FC) -> markers_double_top4

markers_double_top4 <- markers_double_top4[-c(5),] 
markers_double_top4 <- markers_double_top4[-c(8),] 


DotPlot(fibro_only, features = markers_double_top4$gene)+RotatedAxis()

DotPlot(fibro_only, features = c("Hk2", "Gls", "Has1", "Col12a1","Col3a1", "Cthrc1","Ccl2", "Il6", "Mmp3", "Tsc22d3", "Myoc", "Angptl7", "Igfbp6", "Itm2b"), idents = c("double_pos", "double_neg"))+RotatedAxis()+coord_flip()


DotPlot(fibro_only, features = c("Hk2", "Gls", "Has1", "Col12a1","Col3a1", "Cthrc1","Ccl2", "Il6", "Mmp3", "Tsc22d3", "Myoc", "Angptl7", "Igfbp6", "Itm2b"),)+RotatedAxis()+coord_flip()

```




```{r}

DotPlot(fibro_only, features = c("Ccl2", "Il6", "Mmp3"))


```
```{r}


Idents(fibro_only)<-'Hk2_Glspos_renamed'
fibro_only$Hk2_Glspos_renamed_diseased_state <- paste(Idents(fibro_only), fibro_only$diseased_state, sep = "_")
table(fibro_only$Hk2_Glspos_renamed_diseased_state)


Idents(fibro_only)<-'Hk2_Glspos_renamed_diseased_state'

DefaultAssay(fibro_only)<-'RNA'
GlsHk2_pos_peak_vs_rest_RNA<-FindMarkers(fibro_only, ident.1= "double_pos_peak", ident.2="double_pos_resting")
GlsHk2_pos_resing_vs_rest_RNA<-FindMarkers(fibro_only, ident.1= "double_pos_resolving", ident.2="double_pos_resting")
GlsHk2_pos_resolved_vs_rest_RNA<-FindMarkers(fibro_only, ident.1= "double_pos_resolved", ident.2="double_pos_resting")

write.csv(GlsHk2_pos_peak_vs_rest_RNA, "/rds/projects/m/mahonyc-kitwong-runx1/Monica_Guma_paper2/GlsHk2_pos_peak_vs_rest_RNA.csv")
write.csv(GlsHk2_pos_resing_vs_rest_RNA, "/rds/projects/m/mahonyc-kitwong-runx1/Monica_Guma_paper2/GlsHk2_pos_resing_vs_rest_RNA.csv")
write.csv(GlsHk2_pos_resolved_vs_rest_RNA, "/rds/projects/m/mahonyc-kitwong-runx1/Monica_Guma_paper2/GlsHk2_pos_resolved_vs_rest_RNA.csv")

```
```{r}

GlsHk2_pos_peak_vs_rest_RNA$gene<-rownames(GlsHk2_pos_peak_vs_rest_RNA)

EnhancedVolcano(GlsHk2_pos_peak_vs_rest_RNA,
    lab = rownames(GlsHk2_pos_peak_vs_rest_RNA),
    x = 'avg_log2FC',
    y = 'p_val_adj',
        selectLab = c("Ccl2", "Il6", "Mmp3"),
    title = 'GlsHk2_pos_peak_vs_rest_RNA-RNA',
    subtitle = "GEX, red=p_adj<0.05 & FC > 0.5",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5
    )
```
```{r}
EnhancedVolcano(GlsHk2_pos_resing_vs_rest_RNA,
    lab = rownames(GlsHk2_pos_resing_vs_rest_RNA),
    x = 'avg_log2FC',
    y = 'p_val_adj',
        selectLab = c("Ccl2", "Il6", "Mmp3"),
    title = 'GlsHk2_pos_resing_vs_rest_RNA',
    subtitle = "GEX, red=p_adj<0.05 & FC > 0.5",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5
    )
```
```{r}
EnhancedVolcano(GlsHk2_pos_resolved_vs_rest_RNA,
    lab = rownames(GlsHk2_pos_resolved_vs_rest_RNA),
    x = 'avg_log2FC',
    y = 'p_val_adj',
        selectLab = c("Ccl2", "Il6", "Mmp3"),
    title = 'GlsHk2_pos_resolved_vs_rest_RNA',
    subtitle = "GEX, red=p_adj<0.05 & FC > 0.5",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5
    )
```

```{r}
#heatmap on RNA
DefaultAssay(fibro_only)<-'RNA'

markers_double_f= filter(markers_double, p_val_adj < 0.01 )

genes<-markers_double_f$gene
genes<-unique(genes)

Idents(fibro_only)<-'Hk2_Glspos_renamed'
levels(fibro_only)
Dotplpot_RNA<-DotPlot(fibro_only, features=genes)
Dotplpot_data_RNA<-Dotplpot_RNA[["data"]]

library(tidyr)



Dotplpot_data_RNA<-Dotplpot_data_RNA %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 



row.names(Dotplpot_data_RNA) <- Dotplpot_data_RNA$features.plot  
Dotplpot_data_RNA<-Dotplpot_data_RNA %>% select(features.plot, double_neg, Hk2_pos, Gls_pos, double_pos)

library(dplyr)

Dotplpot_data_RNA <- Dotplpot_data_RNA[,-1] %>% as.matrix()

heatmap_RNA<-Heatmap(Dotplpot_data_RNA,cluster_columns =F)
heatmap_RNA

######check cell proportions
pt_orig <- table(fibro_only$Hk2_Glspos_renamed, fibro_only$orig.ident)
pt_orig <- as.data.frame(pt_orig)
pt_orig$Var1 <- as.character(pt_orig$Var1)

library(tidyverse)
pt_double_pos = subset(pt_orig, Var1 == "double_pos")
pt_double_neg = subset(pt_orig, Var1== "double_neg")
pt_Gls = subset(pt_orig, Var1 == "Gls_pos")
pt_Hk2_pos = subset(pt_orig, Var1 == "Hk2_pos")


pt_double_pos <- pt_double_pos %>% select(-one_of('Var1'))
pt_double_pos <- pt_double_pos %>% remove_rownames %>% column_to_rownames(var="Var2")
pt_master = pt_double_pos
colnames(pt_master) <- "double_pos"
#split up for each cluster
pt_master$pt_double_neg <- pt_double_neg$Freq
pt_master$pt_Gls <- pt_Gls$Freq
pt_master$pt_Hk2_pos <- pt_Hk2_pos$Freq
pt_master <- pt_master/rowSums(pt_master)

pt_master <- as.data.frame(pt_master)
pt_master$condition<-c("rest", "rest", "rest", "peak", "peak", "peak", "resolved", "resolved", "resolved", "resolving", "resolving", "resolving")

ggplot(pt_master, aes(x=factor(condition, level=c('rest', 'peak', 'resolving', 'resolved')), y=double_pos)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("double_pos")

res.aov_double_pos <- aov(double_pos ~ condition, data = pt_master)
summary(res.aov_double_pos)
stats_res.aov_double_pos <- TukeyHSD(res.aov_double_pos)
stats_res.aov_double_pos


ggplot(pt_master, aes(x=factor(condition, level=c('rest', 'peak', 'resolving', 'resolved')), y=pt_double_neg)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("double_neg")

res.aov_pt_double_neg <- aov(pt_double_neg ~ condition, data = pt_master)
summary(res.aov_pt_double_neg)
stats_res.aov_pt_double_neg <- TukeyHSD(res.aov_pt_double_neg)
stats_res.aov_pt_double_neg
