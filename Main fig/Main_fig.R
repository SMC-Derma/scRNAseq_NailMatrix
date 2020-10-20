
#############Package #################
library(Seurat)
library(dplyr)
library(reticulate)
library(ggplot2)


seurat.outdir <- "./SingleCell"
setwd(seurat.outdir)

#######################################

# Set color ---------------------------------------------------------------

use_color <-c( "#e84546", "#c68441", "#c6c640", "#75c047", "#53bcaa", "#419dc4", "#4c5fab", "#6b539f", "#a052a0","#d37fa7",
              "#d9d482", "#a8d37f","#99d5c9", "#aea3d0","#f6b6d3","#8097cd","#d36a6a", "#8dcad7", "#abdda4")


# Main UMAP figure ------------------------------------------------------


DimPlot(ON.integrated, reduction = "umap", group.by = "cell_type", label = TRUE, pt.size = 0.5) + NoLegend()

png(file="Mainseurat.png", width = 600, height = 600)

DimPlot(ON.integrated, group.by = "cell_type", label = F, cols = use_color) +
  NoLegend() + NoAxes()

dev.off()


# Proportion plot ---------------------------------------------------------

Idents(ON.integrated)<- "cell_type"

cell.num <- table(ON.integrated@meta.data$cell_class)
ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

DimPlot(ON.integrated, reduction = "umap",group.by = "cell_type", cols=use_color, label.size = 3, label = T)

prop.cell<- ON.integrated@meta.data %>% dplyr::group_by(orig.ident, cell_type) %>% dplyr::summarise(n=n())
prop.cell <-prop.cell %>% dplyr::group_by(orig.ident) %>% dplyr::mutate(sum = sum(n)) 
prop.cell<- prop.cell %>% mutate(prop = n/sum*100)

prop.cell$orig.ident<- factor(prop.cell$orig.ident, levels =  c("Poly1", "Poly2", "Poly3", "Poly4" ))
levels(prop.cell$orig.ident)


ON.integrated@meta.data$cell_type <- as.factor(ON.integrated@meta.data$cell_type)
ON.integrated@meta.data$cell_type    <- factor(ON.integrated@meta.data$cell_type , 
                                               levels =  c("Basal KC", "Suprabasal KC", "Cornified cells",  "proliferating KC", "Keratinocyte-1",
                                                           "Keratinocyte-2", "Myofibroblast", "Fibroblast", "EC", "LEC",
                                                           "lymphocyte", "Macrophage", "Mast cells", "Eccrine gl", "Eccrine duct", "Chondrocytes",
                                                           "Neural cells", "Melanocyte", "unknown" ))


ggplot(prop.cell, aes(x = orig.ident,  y = prop,  
                      fill = (factor(cell_type)))) + 
  geom_bar(stat = "identity", color="black") + xlab("Pathology") + ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),axis.text.y = element_text(angle = , hjust = 1)) + 
  theme_classic() +
  labs(fill = "Cell type") +scale_fill_manual(values = use_color)

###Volcano plot 
devtools::install_github('kevinblighe/EnhancedVolcano')
BiocManager::install('EnhancedVolcano')

# volcano plot ------------------------------------------------------------------

library(EnhancedVolcano)

OMD.markers<- FindMarkers(ON.integrated, ident.1 =  "RSPO4 OMD", ident.2 = "fibroblast", min.pct = 0.25, logfc.threshold = log(1.3))

Idents(ON.integrated) <- "cell_class"
volcano.markers<- FindMarkers(ON.integrated, ident.1 =  "OMD", ident.2 = "Fibroblast.n", min.pct = 0.25, logfc.threshold = log(1.3))
OMD.markers<- FindMarkers(ON.integrated, ident.1 =  "RSPO4 OMD", min.pct = 0.25, logfc.threshold = log(1.3))


EnhancedVolcano(volcano.markers,
                lab = rownames(volcano.markers),
                x = 'avg_logFC',
                y = 'p_val',
                xlim = c(-2, 2),
                title = 'Fibroblast (1-3) versus Onychofibroblast',
                pCutoff = 0.05,
                FCcutoff = log(1.5),
                pointSize = 1.5,
                labSize = 5.0)

sig.gene<- volcano.markers %>% mutate(Gene = rownames(.)) %>% arrange(desc(avg_logFC))  %>% .[ .$p_val < 0.05,]

#save Genelist

writeClipboard(sig.gene$Gene[sig.gene$avg_logFC>0.4])
#run metascape in https://metascape.org/gp/index.html#/main/step1
#paste a gene list 


