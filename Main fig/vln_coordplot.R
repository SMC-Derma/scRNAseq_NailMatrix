
##Vln_coordplot##

Idents(ON.integrated) <-"cell_type"

vln_function <- function(x){
  VlnPlot(ON.integrated, features=x, assay="integrated", pt.size = 0) + coord_flip() + NoLegend() +
    ylab("") +xlab("") + scale_y_continuous(position = "right", breaks = 0:10) +
    scale_x_discrete(limits=positions[18:1]) +
    theme(text=element_text(size=10),
          axis.text.y = element_blank()) + ggtitle(x)
}

vln_list <- list()
use_gene <- c("KRT1","KRT5","SPRR1B", "MKI67",  "WNT6","SPINK6","RGS5","COL1A1",
              "VWF","LYVE1","PTPRC","AIF1", "TPSB2","SCGB1B2P", "DEFB1","COL9A2", "NRXN1",  "PMEL" )

for(i in 1:length(use_gene)){
  vln_list[[i]] <- vln_function(use_gene[i]) 
}

grid.arrange(arrangeGrob(vln_list[[1]], vln_list[[2]], vln_list[[3]],vln_list[[4]], vln_list[[5]], 
                         vln_list[[6]], vln_list[[7]], vln_list[[8]], vln_list[[9]],
                         vln_list[[10]],vln_list[[11]], vln_list[[12]],
                         vln_list[[13]],vln_list[[14]],vln_list[[15]],vln_list[[16]],vln_list[[17]],vln_list[[18]] , ncol=18), nrow=1)


positions <-c("Basal KC", "Suprabasal KC", "Cornified cells",  "proliferating KC", "Keratinocyte-1",
              "Keratinocyte-2", "Myofibroblast", "Fibroblast", "EC", "LEC",
              "lymphocyte", "Macrophage", "Mast cells", "Eccrine gl", "Eccrine duct", "Chondrocytes",
              "Neural cells", "Melanocyte" )


VlnPlot(ON.integrated, features = "RSPO4", cols = use_color, pt.size = 0 ) + coord_flip() + scale_x_discrete(limits=positions[16:1])
