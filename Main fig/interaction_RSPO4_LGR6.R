devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(pals)
library(scales)
library(ineq)
library(tidyverse)
library(circlize)

ON.integrated$nailkc_niche <- ON.integrated$nailkc
ON.integrated@meta.data$nailkc_niche[ON.integrated@meta.data$nailkc=="Fibro-1"] <- "Fibroblast"
ON.integrated@meta.data$nailkc_niche[ON.integrated@meta.data$nailkc=="Fibro-2"] <- "Fibroblast"
ON.integrated@meta.data$nailkc_niche[ON.integrated@meta.data$nailkc=="Fibro-3"] <- "Fibroblast"

#ligand-receptor data 
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

# interactions and their weights in the ligand-receptor + signaling network
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

ON.integrated@active.assay <-"RNA"
Idents(ON.integrated) <- "nailkc_niche"
nailniche.markers <- FindAllMarkers(object = ON.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

# Circosplot --------------------------------------------------------------


expression <- t( ON.integrated@assays$RNA@data)

#set Receptor cluster and ligand cluster 
Rec_cluster <- "LGR6 KC"
ligA_cluster <- "Fibroblast" 
ligB_cluster <- "OMD"

Receptor_ids <-  ON.integrated@meta.data %>% filter(nailkc_niche == Rec_cluster) %>% rownames()
ligA_ids <- ON.integrated@meta.data %>% filter(nailkc_niche == ligA_cluster) %>% rownames()
ligB_ids <- ON.integrated@meta.data %>% filter(nailkc_niche == ligB_cluster) %>% rownames()

expressed_genes_receptor<- expression[Receptor_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()

DE_receiver <-  nailniche.markers %>% filter(cluster==Rec_cluster)
add_geneset <-  DE_receiver %>% filter(p_val_adj <= 0.01 & abs(avg_logFC) >= 0.4) %>% pull(gene)
add_geneset = add_geneset %>% .[. %in% rownames(ligand_target_matrix)]

expressed_genes_receptor <- unique(c(expressed_genes_receptor, add_geneset))

expressed_genes_ligA <-expression[ligA_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_ligB <-expression[ligB_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
## receiver
background_expressed_genes = expressed_genes_receptor %>% .[. %in% rownames(ligand_target_matrix)]

#Perform NicheNet’s ligand activity analysis on the gene set of interest
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands_ligA = intersect(ligands,expressed_genes_ligA)
expressed_ligands_ligB = intersect(ligands,expressed_genes_ligB)
expressed_ligands = union(expressed_ligands_ligA, expressed_ligands_ligB)

receptors = lr_network %>% pull(to) %>% unique()

expressed_receptors = intersect(receptors,expressed_genes_receptor)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)


#Define a gene set of interest
DE_receiver <-  nailniche.markers%>% filter(cluster==Rec_cluster)

geneset_oi <- DE_receiver %>% filter(p_val_adj <= 0.01 & abs(avg_logFC) >= 0.4) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#if u want to check all set of interaction 
#geneset_oi <- colnames(ligand_target_matrix)

head(geneset_oi)

#Perform NicheNet ligand activity analysis
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)


ligand_activities %>% arrange(-pearson)

best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

best_upstream_ligands %>% intersect(expressed_ligands_ligA) 
best_upstream_ligands %>% intersect(expressed_ligands_ligB) 


ligand_expression_tbl = tibble(
  ligand = best_upstream_ligands, 
  ligA = expression[ligA_ids, best_upstream_ligands] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
  ligB = expression[ligB_ids, best_upstream_ligands] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}))


ligA_specific_ligands = ligand_expression_tbl %>% filter(ligA > ligB + 2) %>% pull(ligand)
ligB_specific_ligands = ligand_expression_tbl %>% filter(ligB > ligA + 2) %>% pull(ligand)
general_ligands = setdiff(best_upstream_ligands,c(ligA_specific_ligands,ligB_specific_ligands))

ligand_type_indication_df = tibble(
  ligand_type = c(rep("ligA-specific", times = ligA_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length()),
                  rep("ligB-specific", times = ligB_specific_ligands %>% length())),
  ligand = c(ligA_specific_ligands, general_ligands, ligB_specific_ligands))


#Infer target genes of top-ranked ligands and visualize in a circos plot

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "Rec_cluster" ) %>% inner_join(ligand_type_indication_df) 
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.8)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "#a1d99b",
                   "ligA-specific" = "gold",
                   "ligB-specific" = "#756bb1")
grid_col_target =c( "Rec_cluster" = "#feb24c")


grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

target_order = circos_links$target %>% unique()
ligand_order = c(ligA_specific_ligands,general_ligands,ligB_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)


width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  #width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "ligA-specific") %>% distinct(ligand) %>% nrow() -1)),
  #width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "ligB-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Rec_cluster") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 

circos.clear()

#Render the circos plot 

circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #

# Visualize ligand-receptor interactions of the prioritized ligand --------

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks 
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)
lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "Rec_cluster") %>% inner_join(ligand_type_indication_df)


grid_col_ligand =c("General" = "#a1d99b",
                   "ligA-specific" = "gold",
                   "ligB-specific" = "#756bb1")
grid_col_receptor =c(
  "Rec_cluster" = "#feb24c")




grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 


receptor_order = circos_links$receptor %>% unique()
ligand_order = c(ligA_specific_ligands,general_ligands,ligB_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  # width_ligand_receptor,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "ligA-specific") %>% distinct(ligand) %>% nrow() -1)),
  #width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "ligB-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_receptor,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "Rec_cluster") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
)

circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #

circos.clear()

circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #

circos.clear()

