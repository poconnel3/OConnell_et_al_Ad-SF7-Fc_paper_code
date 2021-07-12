#BY:Patrick O'Connell
#04302021

#This analysis starts from FCS files and goes from there. W/ downsampling. 

##This is a repeat with alternate gating. This gating include ALL cells even those w/ large FSC v.SSC scatter profiles. 


#Packages
library(CATALYST)
library(tidyverse)
library(flowCore)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(diffcyt)
library(xlsx)


##-------------------------Making metadata file----------------------------------##
#Load the content of the metadata.xlsx file into R (file name must match name of the fcs file you generate for each sample)
setwd("~/B16_TME_profiling_batch2_analys")
metadata_filename <- "TME_metadata.xlsx"
md <- read_excel(metadata_filename)


# Define condition variables as named in metadata
md$response <- factor(md$response, levels = c("yes", "no"))
md$genotype <- factor(md$genotype, levels = c("WT", "SLAMF7_KO"))
md$condition <- factor(md$condition, levels = c("Ad_null", "Ad_SLAMF7_Fc"))
md
##-----------------------------------------------------------##

#Read FCS files in as flowset 
setwd("~/B16_TME_profiling_batch2_analys/clean_fcs_files")
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

##------------------making panel file-----------------------##
setwd("~/B16_TME_profiling_batch2_analys")
panel <- "TME_panel.xlsx"
panel <- read_excel(panel) 

#check
all(panel$fcs_colname %in% colnames(fcs_raw))
##------------------making panel file-----------------------##

#see how many cells per file
fsApply(fcs_raw, nrow)


#downsample to 30,000 cells per sample
# Define a downsampling ceiling
sampling.ceiling <- 30000
# Being reproducible is a plus
set.seed(666)

# BUILD A DOWNSAMPLED FLOWSET
fcs_raw.dsamp <- fsApply(fcs_raw, function(ff) {
  idx <- sample.int(nrow(ff), min(sampling.ceiling, nrow(ff)))
  ff[idx,]  
})

#check it
fcs_raw.dsamp
fsApply(fcs_raw.dsamp, nrow)

#make Catalyst data object and arcsinh t-for all channels by 6000

#First 3 must be present and sample_id must match condition
factors <- list(factors = c("file_name", "condition", "sample_id", "response", "genotype"))

daf_X <- prepData(
  fcs_raw.dsamp,
  features = NULL,
  md = md,
  md_cols = factors,
  panel = panel,
  transform = TRUE,
  cofactor = 6000,
  by_time = FALSE,
  FACS = TRUE
)


# view number of events per sample
table(daf_X$condition)
table(daf_X$sample_id)
table(daf_X$file_name)
table(daf_X$response)

# view non-mass channels
names(int_colData(daf_X))

#view metadata parameters
names(colData(daf_X))

#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")    

#plot global marker expression by condition
p <- plotExprs(daf_X, color_by = "condition")
p$facet$params$ncol <- 6                   
p 

#plot number of cells per condition
plotCounts(daf_X, color_by = "genotype")
table(daf_X$condition)

#Make MDS plot  
pbMDS(daf_X, by = "sample_id", color_by = "condition", shape_by = "genotype", size_by = TRUE, label_by = NULL)

#Heatmap of all markers per sample_id and condition w/ heiracheal clustering.
plotExprHeatmap(daf_X, 
                bin_anno = TRUE, 
                scale = "first",
                row_anno = "condition"
                )

#Perform FlowSOM clustering
daf_X <- cluster(daf_X, 
               features = "type", 
               xdim = 10, 
               ydim = 10, 
               maxK = 22, 
               seed = 1234, 
               verbose = TRUE
               ) 

#plot cluster heatmap
plotClusterHeatmap(daf_X,                               
                   hm2 = NULL, #swap for "state" or "abundances" to get more info. 
                   k = "meta22", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
                   ) 

#plot marker expression by cluster
plotClusterExprs(daf2, k = "meta20")  

# run UMAP                          
set.seed(777)                               
daf_X <- runDR(daf_X, 
             dr = "UMAP", 
             cells = 10000,
             features = "type" 
             )

# run TSNE                          
set.seed(777)                               
daf_X <- runDR(daf_X, 
               dr = "TSNE", 
               cells = 10000,
               features = "type" 
)

#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")

#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")



#UMAP visualizations
pz <- plotDR(daf_X_clean, "UMAP", color_by = "meta22")               
pz + theme(axis.line = element_line(colour = NA),   #use this theme for all plots
    axis.ticks = element_line(colour = NA), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(colour = NA), 
    axis.text = element_text(colour = NA), 
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 14), 
    panel.background = element_rect(fill = NA)) 


## Facet per sample                                          
plotDR(daf, "UMAP", color_by = "meta20", facet = "sample_id")

## Facet per condition                                       
plotDR(daf_alt2, "UMAP", color_by = "meta20", facet = "IV_pos")

#plot FlowSOM codes to visualize similarity of clusters
plotCodes(daf_X, k = "meta22")

#Merge and annotate clusters
merging_tableX <- "merged_clusters.xlsx"
merging_table1 <- read_excel(merging_tableX)
head(data.frame(merging_table1))  

# convert to factor with merged clusters in desired order                
merging_table1$new_cluster <- factor(merging_table1$new_cluster,         
                                     levels = c("pDCs", "DCs", "Unknown myeloid", "Neutrophils", "TAMs", "Pro-B cell derived Macrophages", "Monocytes",            
                                                "CD8+ T cells", "CD4+ T cells", "NK cells", "B cells", "Debris"))        

daf_X <- mergeClusters(daf_X, k = "meta22",                                  
                     table = merging_table1, 
                     id = "merging1",
                     overwrite = TRUE
                     )   

#remove debris cluster and CCR2_high mono (there are none that I can appreciate)
#make a new SCE object when you do this!!
to_remove <- c("Debris", "CCR2_high Monocytes")
`%notin%` <- Negate(`%in%`) #make this opposite of %in% operator to help
daf_X_clean <- filterSCE(daf_X, cluster_id %notin% to_remove, k = "merging1")

#view order of clusters
levels(cluster_ids(daf_X_clean, k = "merging1"))

#change cluster colors   
clust_color <- c("#AA00FF", "#FF6597", "#717171", "#5cd6d6", "#CCCC00", "#006699", "#00CC66", "#FF97FC",
                 "#FFAA00", "#000000", "#990033")

#UMAP of merged clusters
plotDR(daf_X_clean, "UMAP", color_by = "merging1", k_pal = clust_color) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
                                                                axis.ticks = element_line(colour = NA), 
                                                                panel.grid.major = element_line(linetype = "blank"), 
                                                                panel.grid.minor = element_line(linetype = "blank"), 
                                                                axis.title = element_text(colour = NA), 
                                                                axis.text = element_text(colour = NA), 
                                                                legend.text = element_text(size = 10), 
                                                                legend.title = element_text(size = 14), 
                                                                panel.background = element_rect(fill = NA)) 

#heatmap of merged cluster markers 
plotExprHeatmap(daf_X_clean, 
                bin_anno = FALSE, 
                scale = "last",
                k = "merging1",
                by = "cluster_id", 
                k_pal = clust_color, 
                hm_pal = brewer.pal(9, "Greys"),
                features = c("SLAMF7", "CD38", "CD8", "CD90", "CD3", "CD11b", "Ly6C", "CD4", "AF", "F4/80", "Fascin",
                               "B220", "CD19", "CCR2", "NK1.1", "CD11c", "Ly6G", "MHC-II", "IgD", "PD-L1", "CCR7", "CD40", "PD-1")
)

#make UMAP colored by important markers 
plotDR(daf_X_clean, "UMAP", color_by = c("CD38", "NK1.1", "CD8", "CD4", "B220", "CD11c", "Ly6G", "MHC-II", "F4/80", "SLAMF7", "CD19", "CD11b", "CD90", "PD-1"), facet_by = NULL, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)


#save main analysis object
saveRDS(daf_X_clean, file = "./daf_X_clean.rds")

#-------------------------------compare groups now------------------------------##

#plot abundances as barplot
plotAbundances(daf_X_clean, k = "merging1", by = "sample_id", k_pal = clust_color)


#plot abundances as boxplot (WT V. KO)
daf_X_null <- filterSCE(daf_X_clean, condition == "Ad_null")
plotAbundances(daf_X_null, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL,
               )


##To analyze now we can generate a GLMM
ei <- metadata(daf_X_null)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))


da_res1 <- diffcyt(daf_X_null,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  


cluster_stats <- rowData(da_res1$res) %>%
  as.data.frame(.)

#export stat data
write.xlsx(cluster_stats, "~/B16_TME_profiling_batch2_analys/cluster_stats.xlsx")

FDR_cutoff <- 0.05

table(rowData(da_res1$res)$p_adj < FDR_cutoff)


#plot abundances as boxplot (Ad_null V. Ad-SF7-Fc)
daf_X_trt <- filterSCE(daf_X_clean, genotype != "SLAMF7_KO")
plotAbundances(daf_X_trt, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "condition",
               shape = NULL,
)

plotAbundances(daf_X_clean, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "condition",
               shape = "genotype"
)

##To analyze now we can generate a GLMM
ei <- metadata(daf_X_trt)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "condition",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))


da_res5 <- diffcyt(daf_X_trt,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  


cluster_stats5 <- rowData(da_res5$res) %>%
  as.data.frame(.)

#export stat data
write.xlsx(cluster_stats5, "~/B16_TME_profiling_batch2_analys/cluster_stats5.xlsx")


table(rowData(da_res5$res)$p_adj < FDR_cutoff)



#plot expression of markers by cluster (boxplot) (by genotype)
plotPbExprs(daf_X_null,
            k = "merging1",
            features = NULL,
            facet_by = "cluster_id",
            color_by = "genotype",
            group_by = "genotype",
            geom = "both", 
            jitter = TRUE
                   )

#plot expression of markers by cluster (boxplot) (by traetment)
plotPbExprs(daf_X_trt,
            k = "merging1",
            features = NULL,
            facet_by = "cluster_id",
            color_by = "condition",
            group_by = "condition",
            geom = "both", 
            jitter = TRUE
)

##---------now make histogram plots of SLAMF7 expression for each cluster comparing WT and SLAMF7-KO------##
#manually extract df containing merging1 cluster ids matched to SLAMF7 expression

#add merging1 cluster to metadata
daf_X_null$merg_clust <- cluster_ids(daf_X_null, "merging1")
#make df
ggdf <- data.frame(colData(daf_X_null), t(assay(daf_X_null, "exprs")), check.names = FALSE)
#make plot
fg <- ggplot(ggdf, aes_string("SLAMF7", fill = "genotype")) + 
  geom_density(alpha = 0.7) + 
  facet_wrap(ggdf$merg_clust, scales = "free", nrow = 2) +
  scale_fill_manual(values = c("#FF3333", "#636363"))

fg + theme(axis.line = element_line(size = 0.5, 
                                    linetype = "solid"), panel.background = element_rect(fill = NA)) 



##Statistically compare markers on clusters w/ GLMM. 
ei <- metadata(daf_X_null)$experiment_info

ds_formula2 <- createFormula(ei, cols_fixed = "genotype")  #NOTE: cannot have a rondom column variable here or else it fails for some reason. 

contrast <- createContrast(c(0, 1))

FDR_cutoff <- 0.1

ds_res2 <- diffcyt(daf_X_null, 
                   formula = ds_formula2, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)

table(rowData(ds_res2$res)$p_adj < FDR_cutoff)

cluster_stats4 <- rowData(ds_res2$res) %>%
  as.data.frame(.)
#export stat data
write.xlsx(cluster_stats4, "~/B16_TME_profiling_batch2_analys/cluster_stats4.xlsx")

plotDiffHeatmap(daf_X_null, rowData(ds_res2$res), top_n = 20, fdr = FDR_cutoff, lfc = 2, col_anno = c("response", "genotype"), row_anno = TRUE) 
##---------------------------------------------------------------##



##---------------------Now lets subset out DCs and look for mregDCs----------------------##

#load main analysis object back in
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")

#subset DCs cells
daf_DCs <- filterSCE(daf_X_clean, cluster_id == "DCs", k = "merging1")

levels(daf_DCs$condition)
levels(daf_DCs$sample_id)

#Make MDS plot  
pbMDS(daf_DCs, by = "sample_id", color_by = "condition", shape_by = "genotype", size_by = TRUE, label_by = NULL)

#re-cluster just the DCs (i am using just markers  relevant to DCs to cluster here)
daf_DCs <- cluster(daf_DCs, 
                 features = c("CD45", "CD11b", "CD4", "CD8", "Ly6C", "MHC-II", "CD11c", "B220", "CCR2", "CD38", "AF", "CD40", "PD-L1", "Fascin", "CCR7"), 
                 xdim = 10, 
                 ydim = 10, 
                 maxK = 12, 
                 seed = 7689, 
                 verbose = TRUE
) 

#plot cluster heatmap
plotClusterHeatmap(daf_DCs,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta12", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
) 

# run UMAP                          
set.seed(777)                               
daf_DCs <- runDR(daf_DCs, 
               dr = "UMAP", 
               cells = 10000,
               features = c("CD45", "CD11b", "CD4", "CD8", "Ly6C", "MHC-II", "CD11c", "B220", "CCR2", "CD38", "AF", "CD40", "PD-L1", "Fascin", "CCR7"), 
)

#UMAP DCs 
plotDR(daf_DCs, "UMAP", color_by = "meta12") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)


#generating scatter plots of some mregDCs markers on all these cells so I can best identify true mregDC clusters
plotScatter(
  daf_DCs,
  c("CCR7", "Fascin"),
  color_by = "meta12",
  facet_by = NULL,
  bins = 100,
  assay = "exprs",
  label = c("target"),
  zeros = TRUE
) + geom_point(size=0.5) + facet_wrap("genotype")


#Merge and annotate DC clusters
merging_table_DC <- "merged_clusters_DC.xlsx"
merging_table2 <- read_excel(merging_table_DC)
head(data.frame(merging_table2))  

# convert to factor with merged clusters in desired order                
merging_table2new_cluster <- factor(merging_table2$new_cluster,         
                                     levels = c("DC2 B", "DC2 A", "mregDCs", "DC2 C", "DC1 A", "DC1 B", "DC2 D", "DC2 E", "DC2 F"))        

daf_DCs_merg1 <- mergeClusters(daf_DCs, k = "meta12",                                  
                       table = merging_table2, 
                       id = "DC_clusters",
                       overwrite = TRUE
)   


#save objects
saveRDS(daf_DCs_merg1, file = "./daf_DC_merg1.rds")
saveRDS(daf_DCs, file = "./daf_DCs.rds")

#view order of clusters
levels(cluster_ids(daf_DCs_merg1, k = "DC_clusters"))

#change cluster colors
clust_color4 <- c("#00994D", "#3332FE", "#8800CC", "#FF3254", "#CCCC00", "#33BAFE", "#ffad33", "#29a3a3", "#8c8c8c")

#UMAP DC cells
plotDR(daf_DCs_merg1, "UMAP", color_by = "DC_clusters", k_pal = clust_color4) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#heatmap of merged DC cluster markers
plotExprHeatmap(daf_DCs_merg1, 
                bin_anno = FALSE, 
                scale = "last",
                k = "DC_clusters",
                by = "cluster_id", 
                k_pal = clust_color4,
                features = c("SLAMF7", "CD38", "PD-L1", "CD11b", "Ly6C", "AF", "CD4", "F4/80", "MHC-II", "B220", "CCR2", "CCR7", "CD40", "Fascin", "Lyve-1"),
                hm_pal = brewer.pal(9, "Greys")
)

#scatter plot of mregDC markers
plotScatter(
  daf_DCs_merg1,
  c("CCR7", "Fascin"),
  color_by = "DC_clusters",
  facet_by = NULL,
  bins = 100,
  assay = "exprs",
  label = c("target"),
  k_pal = clust_color4,
  zeros = TRUE
) +  geom_point(size=0.5)



#make UMAP colored by important markers on just mregDCs cells
plotDR(daf_DCs_merg1, "UMAP", color_by = c("CD38", "SLAMF7", "Ly6C", "F4/80", "CD40", "PD-L1", "Fascin", "CD11b", "MHC-II", "CCR7", "CCR2", "AF", "CD4", "B220"), facet_by = NULL, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#------------------compare DC cluster frequencies by condition-----------------------------#

#plot abundances as barplot
plotAbundances(daf_DCs_merg1, k = "DC_clusters", by = "sample_id", k_pal = clust_color4)


#plot abundances as boxplot (by gneotype)
daf_DCs_null <- filterSCE(daf_DCs_merg1, condition == "Ad_null")
plotAbundances(daf_DCs_null, 
               k = "DC_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               k_pal = clust_color4,
               shape = NULL
)


##To analyze now we can generate a GLMM
eil <- metadata(daf_DCs_null)$experiment_info 
(da_formula13 <- createFormula(eil,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res13 <- diffcyt(daf_DCs_null,                                            
                   formula = da_formula13, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "DC_clusters", verbose = FALSE)  

#examine output
cluster_stats3 <- rowData(da_res13$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res13$res)$p_adj < FDR_cutoff)

#export stat data
write.xlsx(cluster_stats3, "~/B16_TME_profiling_batch2_analys/cluster_stats_DCs.xlsx")


#plot abundances as boxplot (by treatment)
daf_DCs_trt <- filterSCE(daf_DCs_merg1, genotype != "SLAMF7_KO")
plotAbundances(daf_DCs_trt, 
               k = "DC_clusters", 
               by = "cluster_id",
               group_by = "condition",
               k_pal = clust_color4,
               shape = NULL
)


##To analyze now we can generate a GLMM
eil <- metadata(daf_DCs_trt)$experiment_info 
(da_formula13 <- createFormula(eil,                     
                               cols_fixed = "condition",      
                               cols_random = "sample_id")) 


da_res13 <- diffcyt(daf_DCs_trt,                                            
                    formula = da_formula13, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "DC_clusters", verbose = FALSE)  

#examine output
cluster_stats6 <- rowData(da_res13$res) %>%
  as.data.frame(.)


table(rowData(da_res13$res)$p_adj < FDR_cutoff)

#export stat data
write.xlsx(cluster_stats3, "~/B16_TME_profiling_batch2_analys/cluster_stats_DCs.xlsx")

daf_DCs_merg1 <- readRDS(file = "./daf_DC_merg1.rds")

#plot all
plotAbundances(daf_DCs_merg1, 
               k = "DC_clusters", 
               by = "cluster_id",
               group_by = "condition",
               k_pal = clust_color4,
               shape = "genotype"
)

##--------------------------------------------------------------------------------------##



##------------------Now subset TAMs-----------------------------------------------------##

#load main analysis object back in
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")

#subset TAMS cells
daf_TAM <- filterSCE(daf_X_clean, cluster_id == c("TAMs"), k = "merging1")

levels(daf_TAM$genotype)


#Make MDS plot  
pbMDS(daf_TAM, by = "sample_id", color_by = "condition",  size_by = TRUE, shape_by = "genotype", label_by = NULL)

#re-cluster just the TAMs     
daf_TAM <- cluster(daf_TAM, 
                 features = c("CD45", "CD11b", "CD4", "F4/80", "Ly6C", "MHC-II", "CCR2", "CCR7", "CD40", "CD38", "CD206", "Lyve-1", "SLAMF7", "PD-L1"), 
                 xdim = 10, 
                 ydim = 10, 
                 maxK = 10, 
                 seed = 7089, 
                 verbose = TRUE
) 

#plot cluster heatmap
plotClusterHeatmap(daf_TAM,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta10", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
) 

# run UMAP                          
set.seed(790)                               
daf_TAM <- runDR(daf_TAM, 
               dr = "UMAP", 
               cells = 10000,
               features = c("CD45", "CD11b", "CD4", "F4/80", "Ly6C", "MHC-II", "CCR2", "CCR7", "CD40", "CD38", "CD206", "Lyve-1", "SLAMF7", "PD-L1"), 
)

#UMAP TAMs
plotDR(daf_TAM, "UMAP", color_by = "meta10") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#Merge and annotate TAM clusters
merging_table_TAM <- "merged_clusters_TAM.xlsx"
merging_table4 <- read_excel(merging_table_TAM)
head(data.frame(merging_table4))  

# convert to factor with merged clusters in desired order                
merging_table4$new_cluster <- factor(merging_table4$new_cluster,         
                                     levels = c("TAM A", "TAM B", "TAM C", "TAM D", "TAM E", "CD206_high TAMs", "SLAMF7_high CD38_high TAMs"))        

daf_TAM <- mergeClusters(daf_TAM, k = "meta10",                                  
                       table = merging_table4, 
                       id = "TAM_clusters",
                       overwrite = TRUE
)   

levels(cluster_ids(daf_TAM, k = "TAM_clusters"))

#save objects
saveRDS(daf_TAM, file = "./daf_TAM.rds")

#change cluster colors
clust_color6 <- c("#FFAA00", "#AA00FF", "#555555", "#FF65FD", "#66CBFD", "#00994D", "#FF6666")

#UMAP TAM
plotDR(daf_TAM, "UMAP", color_by = "TAM_clusters", k_pal = clust_color6, facet_by = "genotype", a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#heatmap of merged tam cluster markers
plotExprHeatmap(daf_TAM, 
                bin_anno = FALSE, 
                scale = "last",
                k = "TAM_clusters",
                by = "cluster_id", 
                k_pal = clust_color6,
                features = c("CD11b", "CD4", "Ly6C", "CD38", "AF", "CD206", "Lyve-1", "F4/80", "Ly6C", "MHC-II", "CCR2", "CCR7", "CD206", "SLAMF7", "CD40", "PD-L1"),
                hm_pal = brewer.pal(9, "Greys")
)

#UMAP TAMs (by important markers)
plotDR(daf_TAM, "UMAP", color_by = c("CD38", "Ly6C", "AF", "CD206", "CD40", "CD11b", "F4/80", "MHC-II", "CCR2", "CCR7", "SLAMF7", "PD-L1"), k_pal = clust_color6, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#generating scatter plots of clusters by SLAMF7 and PD-L1
plotScatter(
  daf_TAM,
  c("SLAMF7", "PD-L1"),
  color_by = "TAM_clusters",
  facet_by = NULL,
  k_pal = clust_color6,
  bins = 100,
  assay = "exprs",
  label = c("target"),
  zeros = TRUE
) + geom_point(size=0.5)

##  CALCULATE CORRELATION COEFFICIENT B/W SF7 AND PD-L1 ON ALL TAMS   ##
library(corrplot)
SF7_pdl1_expr <- data.frame(colData(daf_TAM), t(assay(daf_TAM, "exprs")), check.names = FALSE)

#  calculate correlation coefficients 
cor_matrix_1 = cor(SF7_pdl1_expr[ ,c(15,20)], method='pearson',use='pairwise.complete.obs')
# Now produce the plot
all_cor_matrix1 <- corrplot(cor_matrix_1, method='circle', type='upper', order="hclust")
all_cor_matrix1

#------------------compare TAM cluster frequencies by condition-----------------------------#
daf_TAM <- readRDS(file = "./daf_TAM.rds")

#plot abundances as barplot
plotAbundances(daf_TAM, k = "TAM_clusters", by = "sample_id", k_pal = clust_color6)


#plot abundances as boxplot (by genotype)
daf_TAM_trt <- filterSCE(daf_TAM, condition != "Ad_SLAMF7_Fc")
obj6 <- plotAbundances(daf_TAM_trt, 
               k = "TAM_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL)
obj6 + scale_color_manual(values = c("#ff3399", "#009999")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA)


##To analyze now we can generate a GLMM
eiw <- metadata(daf_TAM_trt)$experiment_info 
(da_formula30 <- createFormula(eiw,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res30 <- diffcyt(daf_TAM_trt,                                            
                    formula = da_formula30, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "TAM_clusters", verbose = FALSE)  

#examine output
cluster_stats30 <- rowData(da_res30$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res30$res)$p_adj < FDR_cutoff)


#plot abundances as boxplot (by treatment)
daf_TAM_alt <- filterSCE(daf_TAM, genotype != "SLAMF7_KO")
plotAbundances(daf_TAM_alt, 
               k = "TAM_clusters", 
               by = "cluster_id",
               group_by = "condition",
               k_pal = clust_color6,
               shape = NULL
)


##To analyze now we can generate a GLMM
eiw <- metadata(daf_TAM_alt)$experiment_info 
(da_formula30 <- createFormula(eiw,                     
                               cols_fixed = "condition",      
                               cols_random = "sample_id")) 


da_res30 <- diffcyt(daf_TAM_alt,                                            
                    formula = da_formula30, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "TAM_clusters", verbose = FALSE)  

#examine output
cluster_stats30 <- rowData(da_res30$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res30$res)$p_adj < FDR_cutoff)

#plot all
plotAbundances(daf_TAM, 
               k = "TAM_clusters", 
               by = "cluster_id",
               group_by = "condition",
               k_pal = clust_color6,
               shape = "genotype"
)
##--------------------------------------------------------------------------------------##


