
#############Package #################
library(Seurat)
library(dplyr)
library(reticulate)
library(ggplot2)


seurat.outdir <- "./SingleCell"
setwd(seurat.outdir)

# Subclustering of keratinocytes identifies the nail-specific cell --------

Idents(ON.integrated) <-"cell_type"
ON.krt <- subset(x=ON.integrated, idents=c( "Suprabasal KC" ,"Basal KC",
                                            "proliferating KC", "Keratinocyte-1", "Keratinocyte-2", 
                                            "Cornified cells"))

ON.krt <- RunPCA(ON.krt, npcs = 30, verbose = FALSE)
print(ON.krt[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(ON.krt, reduction = "pca")
DimHeatmap(ON.krt, dims = 1:13, cells = 1000, balanced = TRUE)
ElbowPlot(ON.krt, ndims = 30)

ON.krt <- RunUMAP(ON.krt, reduction = "pca", dims = 1:30)

DimPlot(ON.krt, reduction = "umap", label = T, group.by = "orig.ident")
DimPlot(ON.krt, reduction = "umap", 
        label = T)

ON.krt <- RunTSNE(object = ON.krt, dims = 1:30)
ON.krt <- FindNeighbors(ON.krt, dims = 1:20)
ON.krt <- FindClusters(ON.krt, resolution =0.6)


krt.markers <- FindAllMarkers(object = ON.krt, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
seurat.top40.genes <- seurat.markers %>% group_by(cluster) %>% top_n(n = 60, wt = avg_logFC)

ON.nailkc <- subset(x=ON.integrated, idents=c( "Keratinocyte-1" ,"Keratinocyte-2"))
ON.nailkc <- RunPCA(ON.nailkc, npcs = 30, verbose = FALSE)
print(ON.nailkc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(ON.nailkc, reduction = "pca")
DimHeatmap(ON.nailkc, dims = 1:13, cells = 1000, balanced = TRUE)
ElbowPlot(ON.nailkc, ndims = 30)
ON.nailkc <- RunUMAP(ON.nailkc, reduction = "pca", dims = 1:20)

ON.nailkc <- RunTSNE(object = ON.nailkc, dims = 1:30)
ON.nailkc <- FindNeighbors(ON.nailkc, dims = 1:20)
ON.nailkc <- FindClusters(ON.nailkc, resolution =0.6)


# LGR6+ nail matrix cells (basal layer of nail matrix) -------------
# were isolated from the nail keratinocyte-1 (nail matrix)  --------

ON.nailkc@meta.data$nailkc <- ON.nailkc@meta.data$cell_type
ON.nailkc@meta.data$nailkc <- as.character(ON.nailkc@meta.data$nailkc)
LGR6.dt <- data.frame(ON.nailkc@assays$RNA["LGR6",])
LGR6.dt <-data.frame(t(LGR6.dt))
LGR6_ID<-rownames(LGR6.dt)[LGR6.dt>0]
LGR6_ID <- gsub("[.]", "-", LGR6_ID)

ON.nailkc@meta.data$nailkc[rownames(ON.nailkc@meta.data) %in% LGR6_ID] <-"LGR6 KC"

Idents(ON.nailkc) <-"nailkc"
VlnPlot(ON.nailkc, features = "LGR6")

ON.integrated@meta.data$nailkc <- as.character(ON.integrated@meta.data$subcelltype)
ON.integrated@meta.data[rownames(ON.nailkc@meta.data),]$nailkc <- ON.nailkc@meta.data$nailkc
Idents(ON.integrated)<- "nailkc"

# Subclustering of fibroblasts identifies the nail-specific cell population (onychofibroblasts)----------------

Idents(ON.integrated) <-"cell_type"
ON.fibro <- subset(x=ON.integrated, idents=c( "Fibroblast"))
ON.fibro <- RunPCA(ON.fibro, npcs = 30, verbose = FALSE)
ON.fibro <- RunUMAP(ON.fibro, reduction = "pca", dims = 1:20)

ON.fibro <- RunTSNE(object = ON.fibro, dims = 1:20)
ON.fibro <- FindNeighbors(ON.fibro, dims = 1:20)
ON.fibro <- FindClusters(ON.fibro, resolution =0.6)

fibro.markers <- FindAllMarkers(object = ON.fibro, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
current.cluster.ids <- levels(ON.fibro)                     
names(new.cluster.ids) <- levels(ON.fibro)
new.cluster.ids <-c("Fibro-1", "OMD", "OMD", "Fibro-2", "Fibro-1", "Fibro-3", 
                    "Fibro-1", "Fibro-1", "OMD", "Fibro-3", "Fibro-2")

ON.fibro<- RenameIdents(ON.fibro, new.cluster.ids)
ON.fibro@meta.data$subcelltype <- plyr::mapvalues(x = ON.fibro@meta.data$seurat_clusters,
                                                  from = current.cluster.ids, to = new.cluster.ids)

ON.fibro@meta.data$subcelltype <- as.character(ON.fibro@meta.data$subcelltype)
DimPlot(ON.fibro, group.by="subcelltype", label=T)
ON.integrated@meta.data$subcelltype <- as.character(ON.integrated@meta.data$cell_type)
ON.integrated@meta.data[rownames(ON.fibro@meta.data),]$subcelltype <- ON.fibro@meta.data$subcelltype

OMD.markers <- FindAllMarkers(object = ON.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
OMD.top60.genes <- OMD.markers %>% group_by(cluster) %>% slice_min(p_val_adj, n=60) %>% slice_max(avg_logFC, n=60)


