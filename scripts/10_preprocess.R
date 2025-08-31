table(grepl("^MT-",rownames(pbmc)))

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(pbmc@meta.data, 5)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
table(pbmc$patient_id,pbmc$Tissue)

table(grepl("^ERCC-",rownames(pbmc)))
pbmc[["percent.ERCC"]] <- PercentageFeatureSet(pbmc, pattern = "^ERCC-")
head(pbmc@meta.data)
summary(pbmc@meta.data)
rownames(pbmc)[grep("^ERCC-",rownames(pbmc))]
sum(pbmc$percent.ERCC< 40)

sum(pbmc$percent.ERCC< 10)   #就只剩下460个cell，明显低于文献中的数量
pctERCC=40
pbmc <- subset(pbmc, subset = percent.ERCC < pctERCC)
dim(pbmc)
# >20047  2142 
dim(a.filt)
# >23460  2343 
pbmc.filt <- subset(pbmc, subset = patient_id=="R1"&sample_type=="PeriTumor")    #BT_S1 2 4 6 Tumor Periphery
dim(pbmc.filt)

plot1 <- FeatureScatter(pbmc.filt, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "patient_id",pt.size=1.5)
plot2 <- FeatureScatter(pbmc.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "patient_id",pt.size=1.5)
plot1 + plot2
############################
selected_f <- rownames(pbmc)[Matrix::rowSums(pbmc@assays$RNA@counts > 0 ) > 3]
pbmc.filt <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 5 & 
                      orig.ident%in%c("X102","X105","X105A","X114","X115","X118","X124","X125","X126","X143"),
                    features = selected_f)
dim(pbmc.filt)
VlnPlot(object = pbmc.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "orig.ident")



pbmc <- NormalizeData(pbmc.filt, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1500)   #变异最大的1500个基因


top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10,repel = TRUE)
plot1+plot2




all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 40)
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(object = pbmc, dims = 1:2, reduction = "pca",nfeatures = 20)  #4个PC 20个基因

DimPlot(pbmc, reduction = "pca",group.by="orig.ident")

DimHeatmap(pbmc, dims = 1:2, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)

pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
JackStrawPlot(object = pbmc, dims = 1:20,reduction = "pca")
ElbowPlot(pbmc,reduction="pca")  #根据ElbowPlot确定PC数量


#devtools::install_github("eddelbuettel/harmony",force = TRUE)

library("harmony")
pbmc=RunHarmony(pbmc,reduction="pca",group.by.vars="orig.ident",reduction.save="harmony")




pbmc <- FindNeighbors(pbmc, reduction = "harmony",dims = 1:15)   
library(clustree)

obj=FindClusters(pbmc, resolution = seq(0.1,1,by=0.1))
clustree(obj)

pbmc <- FindClusters(pbmc, resolution = 0.2)  


pbmc <- RunTSNE(object = pbmc, reduction = "harmony",dims = 1:15)   
DimPlot(pbmc, reduction = "tsne",label = TRUE,pt.size = 1)
pbmc <- RunUMAP(object = pbmc, reduction = "harmony",dims = 1:15)   
DimPlot(pbmc, reduction = "umap",label = TRUE,pt.size = 1)
table(pbmc@meta.data$Tissue,pbmc@meta.data$seurat_clusters)

DimPlot(pbmc, reduction = "umap",group.by = "orig.ident", label = TRUE)


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(pbmc.markers, n = 10)
write.table(pbmc.markers,file = 'pbmc.markers.txt',sep = '\t',row.names=F,quote=F)

library("tidyverse")
sig.markers<- pbmc.markers %>% select(gene,everything()) %>%
  subset(p_val_adj<0.05 & abs(pbmc.markers$avg_log2FC)>1)
dim(sig.markers)
write.table(sig.markers,file="sigmarkers.xls",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

VlnPlot(object = pbmc, features = c("EGFR", "CCL3","AGXT2L1","DCN","GPR17","MOG","SYT1"),group.by = "seurat_clusters",pt.size = 1)
VlnPlot(object = pbmc, features = c("CX3CR1","CD74","TMEM119" ),group.by = "seurat_clusters",pt.size = 1)
saveRDS(pbmc,file="pbmc2.RDS")

pbmc=readRDS("GSE131928tumor.rds")
table(Idents(pbmc))
table(pbmc@meta.data[c('Sample', 'seurat_clusters')])

genes_to_check <- c("EGFR",
                    "PTPRZ1", "MKI67",
                    "TOP2A",
                    "STMN2", "CD24", "CDK4",
                    "PDGFRA","OLIG2",
                    "BCAN",
                    "SOX9",
                    "GFAP",
                    "AQP4",
                    "CHI3L1",
                    "DDIT3")
library(stringr)
library(ggplot2)
genes_to_check = unique(intersect(rownames(pbmc),genes_to_check))
genes_to_check
P13 <- DotPlot(pbmc, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P14 <- VlnPlot(object = pbmc, features =genes_to_check,log =T )
P15 <- FeaturePlot(object = pbmc, features=genes_to_check )
ggsave(filename = "Tumor.pdf",P13/P14/P15,height = 15,width = 10)

genes_to_check <- c("MAG",
                    "MOG", "PLP1"
)

genes_to_check = unique(intersect(rownames(pbmc),genes_to_check))
genes_to_check
P13 <- DotPlot(pbmc, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P14 <- VlnPlot(object = pbmc, features =genes_to_check,log =T )
P15 <- FeaturePlot(object = pbmc, features=genes_to_check )
ggsave(filename = "Normal glial cell.pdf",P13/P14/P15,height = 15,width = 10)

genes_to_check <- c("FCER1A",
                    "CLEC10A", "CD1C", "BACH1", "LYZ", "S100A9", "AIF1", "CD163", "CD68"
                    , "P2RY12", "TMEM119")

genes_to_check = unique(intersect(rownames(pbmc),genes_to_check))
genes_to_check
P13 <- DotPlot(pbmc, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P14 <- VlnPlot(object = pbmc, features =genes_to_check,log =T )
P15 <- FeaturePlot(object = pbmc, features=genes_to_check )
ggsave(filename = "Myeloid cell.pdf",P13/P14/P15,height = 15,width = 10)

genes_to_check <- c("CD3D",
                    "IL7R", "CCR7", "CD8A", "GZMK", "IFNG", "GNLY", "NKG7", "CD19"
                    , "MS4A1", "CD79A", "SDC1")

genes_to_check = unique(intersect(rownames(pbmc),genes_to_check))
genes_to_check
P13 <- DotPlot(pbmc, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P14 <- VlnPlot(object = pbmc, features =genes_to_check,log =T )
P15 <- FeaturePlot(object = pbmc, features=genes_to_check )
ggsave(filename = "Lymphocytes.pdf",P13/P14/P15,height = 15,width = 10)

genes_to_check <- c("CD248",
                    "PDGFRB", "RGS5", "CDH5", "VWF", "PECAM1")

genes_to_check = unique(intersect(rownames(pbmc),genes_to_check))
genes_to_check
P13 <- DotPlot(pbmc, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P14 <- VlnPlot(object = pbmc, features =genes_to_check,log =T )
P15 <- FeaturePlot(object = pbmc, features=genes_to_check )
ggsave(filename = "Vascular cell.pdf",P13/P14/P15,height = 15,width = 10)

levels(pbmc)
celltype <- c("TAM_Microglia","Tumor","Tumor","Tumor","TAM_Machrophage",
                     "Tumor","Normal glial", "Tumor","Lymphocyte")
names(celltype) <- levels(pbmc) 
pbmc <- RenameIdents(pbmc, celltype)       

pbmc$celltype <- Idents(pbmc)
head(pbmc@meta.data)
DimPlot(pbmc,                             
        reduction = 'umap',                
        group.by = 'celltype',               
        pt.size = 1,                       
        #split.by = 'seurat_clusters',       
        label = T)                           



saveRDS(pbmc, file = "GSE.rds")
