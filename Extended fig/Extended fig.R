


# Feature plot with Mean expresssion value --------------------------------


Idents(ON.integrated)<- "cell_type"

cell.markers <- FindAllMarkers(object = ON.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
cell.markers %>% filter(cluster =="Basal KC")
cell.markers %>% filter(cluster =="Suprabasal KC")
cell.markers %>% filter(cluster =="Cornified cells")
cell.markers %>% filter(cluster =="Keratinocyte-1")
cell.markers %>% filter(cluster =="Keratinocyte-2")
cell.markers %>% filter(cluster =="proliferating KC")
cell.markers %>% filter(cluster =="Fibroblast")
cell.markers %>% filter(cluster =="Myofibroblast")
cell.markers %>% filter(cluster =="EC")
cell.markers %>% filter(cluster =="LEC")
cell.markers %>% filter(cluster =="lymphocyte")
cell.markers %>% filter(cluster =="Macrophage")
cell.markers %>% filter(cluster =="Mast cells")
cell.markers %>% filter(cluster =="Eccrine duct")
cell.markers %>% filter(cluster =="Melanocyte")
cell.markers %>% filter(cluster =="Chondrocytes")


rownames(ON.integrated@assays$integrated@scale.data)

geneset <- c("KRT5", "KRT14", "TP63", "S100A2") # Basal keratinocyte
geneset <- c("KRT1", "KRT10", "DSG1", "SBSN" ) # Suprabasal keratinocyte
geneset <- c("SPRR1B", "ACER1", "CALML5" ) # Cornified cells 
geneset <- c("WNT6", "CXCL14", "SPINK6" ) # Keratinocyte-1
geneset <- c("SPINK6", "KRT6A", "KRT16", "KRT17" ) # Keratinocyte-2
geneset <- c( "MKI67", "HMGB2", "TYMS" ) # Proliferating KC
geneset <- c( "COL1A1", "COL1A2", "LUM", "DCN", "SFRP2" ) # Fibroblast
geneset <- c( "RGS5", "ACTA2", "TAGLN", "MYH11") # Myofibroblast
geneset <- c( "VWF", "CD34", "PECAM1") # EC
geneset <- c( "CCL21", "MMRN1", "LYVE1", "CLDN5") # LEC
geneset <- c( "CD3D", "TNFSF14", "PTPRC", "KLRC1") # Immune cells
geneset <- c( "HLA-DRA", "CTSS", "AIF1", "CD68") # Macrophage
geneset <- c( "TPSB2", "TPSAB1", "HPGDS") # Mast cells
geneset <- c( "SCGB1B2P", "DEFB1", "AQP5") # Eccrine gl/duct

geneset <- c( "TYR", "PMEL", "MLANA") # Eccrine gl/duct
geneset <- c( "DKK1", "COL9A2", "PCP4") # Eccrine gl/duct



# Get mean expression of genes of interest per cell
mean.exp <- colMeans(x = ON.integrated@assays$integrated@scale.data[geneset, ], na.rm = TRUE)

# Add mean expression values in 'object@meta.data$gene.set.score'
if (all(names(x = mean.exp) == rownames(x = ON.integrated@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  ON.integrated@meta.data$gene.set.score <- mean.exp
}

# Plot mean expression using Seurat::FeaturePlot()

FeaturePlot(object = ON.integrated, features = "gene.set.score", cols = c("Midnight Blue","Orange", "red"))+NoAxes()+
  ggtitle("")

png(file="AverageExpseurat.png", width = 600, height = 600)

FeaturePlot(object = ON.integrated, features = "gene.set.score", cols = c("Midnight Blue","Orange", "red"))+NoAxes()+
  ggtitle("")

dev.off()


