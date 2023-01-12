#### Single Cell Analysis - Wen Lin & Gabrielle Choonoo - 2023 ####

#### Imports ####
library("Seurat")
library("SeuratWrappers")
library("scRepertoire")
library("scProportionTest")
library("monocle3")
library("ggplot2")
library("dplyr")

#### Package Versions ####
# R version 3.6.3
# Seurat_3.1.4
# SeuratWrappers_0.3.0
# scProportionTest_0.0.0.9000
# monocle3_1.0.0
# ggplot2_3.3.5 
# dplyr_1.0.7

#### Fig 5 - CD4 T cells ####
# Load data
load('CD4.RData')

# Fig 5 A: CD4 clusters
tmp <- DimPlot(seurat.obj_cd4, label.size=4)
LabelClusters(tmp, id="ident", fontface="bold", size=5)

# Fig 5 B: CD4 clusters by tumor
DimPlot(seurat.obj_cd4, pt.size=1, label=F, label.size=20,split.by='bioreplicate_group')+
  theme(
    legend.text=element_text(size=20),
    strip.text = element_text(size=20))

# Fig 5 C: CD4 proportion test
prop_test <- sc_utils(seurat.obj_cd4)
prop_test <- permutation_test(
  prop_test, cluster_identity = "clustname",
  sample_1 = "StRG47-Progressor", sample_2 = "StRG47-Regressor",
  sample_identity = "bioreplicate_group"
)
permutation_plot(prop_test, FDR_threshold = 0.01)

## Fig 5 D: CD4 features
FEATURES <- c('GNLY','GZMB','PRF1','FOXP3','IL2RA','MKI67')
FEATURE_LABELS <- c('GNLY','GZMB','PRF1','FOXP3','IL2RA','MKI67')

sz.fac.wd <- ceiling(sqrt(length(FEATURES)))
sz.fac.ht <- ceiling(length(FEATURES)/sz.fac.wd)

# List of plots, one for each feature
featplot_umap <- FeaturePlot(seurat.obj_cd4, FEATURES, ncol=sz.fac.wd, pt.size=0.5, combine=F)
# For each plot, remove axes and plot margins
for(i in 1:length(featplot_umap)) {
  featplot_umap[[i]] <- featplot_umap[[i]] +
    theme(legend.title=element_blank(),
          plot.title=element_text(size=24),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.margin=unit(c(0, 0, 0, 0), 'pt')) +
    xlab('') + ylab('')+ ggtitle(FEATURE_LABELS[i])
}
# Project grid of all feature plots on UMAP
plot(cowplot::plot_grid(plotlist=featplot_umap, ncol=sz.fac.wd))

## Fig 5 E: CD4 feature violin plots 
# violin plot functions
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),split.by,
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size,split.by=split.by, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size=rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),split.by,
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x,split.by=split.by, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(size=rel(1)), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

stack_v=StackedVlnPlot(obj = seurat.obj_cd4, features = c(
  "GNLY","GZMA","GZMB","PRF1","NKG7", "CCL5","CXCL13", "IFNG","IL7R","TCF7","PDCD1","HAVCR2","IFIT3","KLRB1"),split.by=NULL)
stack_v

## Fig 5 F: CD4 trajectory
# function to set root in trajectory
get_earliest_principal_node <- function(cds, time_bin='C5-Tnaive/Tcm\n559'){
  cell_ids <- which(colData(cds)[, "clustname"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

## Pseudotime UMAP
# filter Tregs
seurat.obj_cd4_subset = subset(seurat.obj_cd4, clustname%in%c('C6-Prolif\n691','C1-PRF1\n1580','C2-IFIT3\n634','C5-Tnaive/Tcm\n559','C4-IL7R\n1764','C3-KLRB1\n308'))

# create trajectory
seurat.obj.cds <- as.cell_data_set(seurat.obj_cd4_subset)
seurat.obj.cds <- cluster_cells(cds = seurat.obj.cds, reduction_method = "UMAP")
seurat.obj.cds <- learn_graph(seurat.obj.cds, use_partition = TRUE)

# order cells
seurat.obj.cds <- order_cells(seurat.obj.cds, reduction_method = "UMAP",root_pr_nodes=get_earliest_principal_node(seurat.obj.cds))

# plot
tmp = plot_cells(
  cds = seurat.obj.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = FALSE,
  label_branch_points=FALSE,
  label_leaves=FALSE,
  graph_label_size=10, cell_size = 0.5)+
  theme(legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.key.size = unit(1, 'cm'),
        plot.title=element_text(size=20),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.margin=unit(c(0, 0, 0, 0), 'pt'))+
  ggtitle('CD4 Tconv')
tmp

## Pseudotime Boxplot
# save pseudotime data from umap 
dat=tmp$data

# rename clusters
dat[,'clustname'] <- gsub('C5-Tnaive/Tcm\n559',"C5-Tnaive/Tcm",dat[,'clustname'])
dat[,'clustname'] <- gsub('C3-KLRB1\n308',"C3-KLRB1",dat[,'clustname'])
dat[,'clustname'] <- gsub('C4-IL7R\n1764',"C4-IL7R",dat[,'clustname'])
dat[,'clustname'] <- gsub('C1-PRF1\n1580',"C1-PRF1",dat[,'clustname'])
dat[,'clustname'] <- gsub('C2-IFIT3\n634',"C2-IFIT3",dat[,'clustname'])
dat[,'clustname'] <- gsub('C6-Prolif\n691',"C6-Prolif",dat[,'clustname'])

# order levels
dat[,'clustname'] <- factor(dat[,'clustname'], levels=c("C5-Tnaive/Tcm","C3-KLRB1","C4-IL7R","C1-PRF1","C2-IFIT3","C6-Prolif"))

# plot
ggplot(dat, aes(reorder(clustname,cell_color,FUN=function(x)median(x,na.rm=T)), cell_color, fill=clustname, color=clustname)) + 
  geom_boxplot(alpha=0.5) + 
  geom_jitter(size=3, width=0.2) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + 
  ylab('Pseudotime') + xlab('CD4 Subsets')


## Fig 5 G: CD4 clone size barplot 
# load data
load('CD4_TCR.RData')

# remove NA
dat_v2[-which(is.na(dat_v2$bin)),] -> dat_v2

# set bins
dat_v2[which(dat_v2$clone.size %in% c(1)),'bin'] <- '1'
dat_v2[which(dat_v2$clone.size %in% c(2:10)),'bin'] <- '2-10'
dat_v2[which(dat_v2$clone.size %in% c(11:30)),'bin'] <- '11-30'
dat_v2[which(dat_v2$clone.size %in% c(31:50)),'bin'] <- '31-50'
dat_v2[which(dat_v2$clone.size %in% c(51:100)),'bin'] <- '51-100'
dat_v2[which(dat_v2$clone.size >100),'bin'] <- '>100'

# summarize data
dat_v3 <- dat_v2 %>% 
  group_by(clustname,bin) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

# order levels
dat_v3$bin <- factor(dat_v3$bin, levels=c('1','2-10','11-30','31-50','51-100','>100'))
dat_v3$clustname <- factor(dat_v3$clustname, levels=c(c('C1-PRF1','C2-IFIT3','C4-IL7R','C3-KLRB1','C6-Prolif',"C5-Tnaive/Tcm",'C7-Treg','C8-Treg-TNFRSF9','C9-Treg-Prolif')))

# plot
ggplot(dat_v3, aes(x=clustname, y=perc*100, fill=bin, group=bin)) +
  geom_bar(position="stack", stat="identity") +
  theme(
    axis.title.x = element_text(face="bold",size=24),
    axis.title.y = element_text(face="bold",size=24),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,face='bold', size=28),
    axis.text.y = element_text(face="bold",size=28)) +
  labs(x = "CD4 subsets", y = "Percentage", fill = "Clone size bin")


## Fig 5 H: CD4 shannon diversity
clonalDiversity(combined_clust, cloneCall = "gene", group = "sample", n.boots = 100) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10, face='bold'))+
  guides(fill=guide_legend(title="Clusters"))



#### Fig S4 - All Cells and CD8 T cells ####

## Load data
load('All.RData')
load('CD8.RData')

## Fig S4 A: hCD45 clusters
tmp <- DimPlot(seurat.obj_all, label.size=4)
LabelClusters(tmp, id="ident", fontface="bold", size=5)

## Fig S4 C: hCD45 proportion test
prop_test <- sc_utils(seurat.obj_all)
prop_test <- permutation_test(
  prop_test, cluster_identity = "clustname",
  sample_1 = "StRG47-Progressor", sample_2 = "StRG47-Regressor",
  sample_identity = "bioreplicate_group"
)
permutation_plot(prop_test, FDR_threshold = 0.01)

## Fig S4 D: CD8 clusters
tmp <- DimPlot(seurat.obj_cd8, label.size=4)
LabelClusters(tmp, id="ident", fontface="bold", size=5)

## Fig S4 E: CD8 features
# set features
FEATURES <- c('GZMA','GZMK','IL7R','GZMB','PRF1','GNLY')
FEATURE_LABELS <- c('GZMA','GZMK','IL7R','GZMB','PRF1','GNLY')
sz.fac.wd <- ceiling(sqrt(length(FEATURES)))
sz.fac.ht <- ceiling(length(FEATURES)/sz.fac.wd)

# loop through features
featplot_umap <- FeaturePlot(seurat.obj_cd8, FEATURES, ncol=sz.fac.wd, pt.size=0.5, combine=F)
for(i in 1:length(featplot_umap)) {
  featplot_umap[[i]] <- featplot_umap[[i]] +
    theme(legend.title=element_blank(),
          plot.title=element_text(size=24),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.margin=unit(c(0, 0, 0, 0), 'pt')) +
    xlab('') + ylab('')+ ggtitle(FEATURE_LABELS[i])
}

# plot
plot(cowplot::plot_grid(plotlist=featplot_umap, ncol=sz.fac.wd))

## Fig S4 F: CD8 proportion test
prop_test <- sc_utils(seurat.obj_cd8)
prop_test <- permutation_test(
  prop_test, cluster_identity = "clustname",
  sample_1 = "StRG47-Progressor", sample_2 = "StRG47-Regressor",
  sample_identity = "bioreplicate_group"
)
permutation_plot(prop_test, FDR_threshold = 0.01)

## Fig S4 G: CD8 clone size barplot
# load data
load('CD8_TCR.RData')

# remove NA
dat_v2[-which(is.na(dat_v2$bin)),] -> dat_v2

# summarize data
dat_v3 <- dat_v2 %>% 
  group_by(clustname,bin) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

# order levels
dat_v3$bin <- factor(dat_v3$bin, levels=c('1','2-10','11-20','21-30','>30'))
dat_v3$clustname <- sapply(strsplit(as.character(dat_v3$clustname), '\n'),"[",1)

# plot
ggplot(dat_v3, aes(x=clustname, y=perc*100, fill=bin, group=bin)) +
  geom_bar(position="stack", stat="identity") +
  theme(
    axis.title.x = element_text(face="bold",size=24),
    axis.title.y = element_text(face="bold",size=24),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,face='bold', size=28),
    axis.text.y = element_text(face="bold",size=28))  +
  labs(x = "CD8 Subsets", y = "Percentage", fill = "Clone size bin")





