#By:Patrick O'Connell
#04032021

#This analysis starts from FCS files and goes from there. W/ downsampling. 
#This is an analysis of the first batch of B16 tumors using the TME panel


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
setwd("~/B16_TME_profiling_batch1_analys_repeat")
metadata_filename <- "TME_metadata.xlsx"
md <- read_excel(metadata_filename)

# Define condition variables as named in metadata
md$condition <- factor(md$condition, levels = c("Ad_Null", "Ad_SLAMF7_Fc"))
md$response <- factor(md$response, levels = c("yes", "no"))
head(data.frame(md))
##-----------------------------------------------------------##

#Read FCS files in as flowset 
setwd("~/B16_TME_profiling_batch1_analys_repeat/clean_fcs_files")
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

##------------------making panel file-----------------------##
setwd("~/B16_TME_profiling_batch1_analys_repeat")
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
factors <- list(factors = c("file_name", "condition", "sample_id", "response"))

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
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta22", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
                   ) 


# run UMAP                          
set.seed(777)                               
daf_X <- runDR(daf_X, 
             dr = "UMAP", 
             cells = 10000,
             features = "type" 
             )


#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")

#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")


#UMAP visualizations
pz <- plotDR(daf_X, "UMAP", color_by = "meta22")               
pz + theme(axis.line = element_line(colour = NA),   #use this theme for all plots
    axis.ticks = element_line(colour = NA), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(colour = NA), 
    axis.text = element_text(colour = NA), 
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 14), 
    panel.background = element_rect(fill = NA)) 



#Merge and annotate clusters
merging_table1 <- "merged_clusters.xlsx"
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))  

# convert to factor with merged clusters in desired order                
merging_table1$new_cluster <- factor(merging_table1$new_cluster,         
                                     levels = c("pDCs", "DCs", "CCR2_high monocytes", "Neutrophils", "TAMs", "Pro-B cell derived Macrophages", "Monocytes",            
                                                "CD8+ T cells", "CD4+ T cells", "NK cells", "B cells", "Debris"))        

daf_X <- mergeClusters(daf_X, k = "meta22",                                  
                     table = merging_table1, 
                     id = "merging1",
                     overwrite = TRUE
                     )   

#remove debris cluster 
#make a new SCE object when you do this!!
to_remove <- c("Debris")
`%notin%` <- Negate(`%in%`) #make this opposite of %in% operator to help
daf_X_clean <- filterSCE(daf_X, cluster_id %notin% to_remove, k = "merging1")

  
#view order of clusters
levels(cluster_ids(daf_X, k = "merging1"))

#change cluster colors    
clust_color <- c("#996600", "#FF6597", "#717171", "#e6e600", "#006699", "#00CC66", "#FF97FC",
                 "#FFAA00", "#000000", "#b84dff", "#cc0000")

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

#heatmap of merged cluster markers (remove and CD45; misleading)
plotExprHeatmap(daf_X_clean, 
                bin_anno = FALSE, 
                scale = "last",
                k = "merging1",
                by = "cluster_id", 
                k_pal = clust_color, 
                hm_pal = brewer.pal(9, "Greys"),
                features = c("SLAMF7", "CD38", "CD8", "CD90", "CD11b", "Ly6C", "CD4", "AF", "F4/80", "Fascin", "CD3", "CD206",
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


#plot abundances as boxplot
obj1 <- plotAbundances(daf_X_clean, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "condition",
               shape = NULL,
               )
obj1 + scale_color_manual(values = c("#a6a6a6", "#000000")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA)


##To analyze now we can generate a GLMM
ei <- metadata(daf_X_clean)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "condition",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res1 <- diffcyt(daf_X_clean,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  


cluster_stats <- rowData(da_res1$res) %>%
  as.data.frame(.)

#export stat data
write.xlsx(cluster_stats, "~/B16_TME_profiling_batch1_analys/cluster_stats.xlsx")

FDR_cutoff <- 0.05

table(rowData(da_res1$res)$p_adj < FDR_cutoff)


#plot abundances as boxplot by response
plotAbundances(daf_X_clean, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "response",
               shape = NULL,
)

##To analyze now we can generate a GLMM (by resposne)
ei <- metadata(daf_X_clean)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "response",      
                              cols_random = "sample_id")) 


da_res2 <- diffcyt(daf_X_clean,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  


cluster_stats2 <- rowData(da_res2$res) %>%
  as.data.frame(.)

#export stat data
write.xlsx(cluster_stats2, "~/B16_TME_profiling_batch1_analys/cluster_stats2.xlsx")

table(rowData(da_res2$res)$p_adj < FDR_cutoff)

#now separate out just responders from both treatments and compare clusters. 
daf_resp <- filterSCE(daf_X_clean, response == "yes")

#plot abundances as boxplot by treatment
plotAbundances(daf_resp, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "condition",
               shape = NULL,
)

##To analyze now we can generate a GLMM (by condition)
ei <- metadata(daf_resp)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "condition",      
                              cols_random = "sample_id")) 


da_res3 <- diffcyt(daf_resp,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  


cluster_stats3 <- rowData(da_res3$res) %>%
  as.data.frame(.)

#export stat data
write.xlsx(cluster_stats3, "~/B16_TME_profiling_batch1_analys/cluster_stats3.xlsx")

table(rowData(da_res3$res)$p_adj < FDR_cutoff)


#plot expression of markers by cluster (boxplot)
plotPbExprs(daf_X_clean,
            k = "merging1",
            features = NULL,
            facet_by = "cluster_id",
            color_by = "condition",
            group_by = "condition",
            geom = "both", 
            jitter = TRUE
                   )


##Statistically compare markers on clusters w/ GLMM. 
ei <- metadata(daf_X_clean)$experiment_info

ds_formula2 <- createFormula(ei, cols_fixed = "condition")  #NOTE: cannot have a rondom column variable here or else it fails for some reason. 

contrast <- createContrast(c(0, 1))

FDR_cutoff <- 0.5

ds_res2 <- diffcyt(daf_X_clean, 
                   formula = ds_formula2, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)

table(rowData(ds_res2$res)$p_adj < FDR_cutoff)

cluster_stats4 <- rowData(ds_res2$res) %>%
  as.data.frame(.)

tbl_DS <- rowData(ds_res2$res)
plotDiffHeatmap(daf_X_clean, tbl_DS, fdr = 0.05, sort_by = "lfc", col_anno = "condition")

##---------------------------------------------------------------##



##---------------------Now lets subset out DCs and look for mregDCs----------------------##

#load main analysis object back in
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")

#subset DCs cells
daf_DCs <- filterSCE(daf_X_clean, cluster_id == "DCs", k = "merging1")

levels(daf_DCs$condition)
levels(daf_DCs$sample_id)

#Make MDS plot  
pbMDS(daf_DCs, by = "sample_id", color_by = "condition", shape_by = "response", size_by = TRUE)

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
               features = c("CD45", "CD11b", "CD4", "CD8", "Ly6C", "MHC-II", "CD11c", "B220", "CCR2", "CD38", "AF", "CD40", "PD-L1", "Fascin", "CCR7")
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
) + geom_point(size=0.5)


#Merge and annotate DC clusters
merging_table_DC <- "merged_clusters_DC.xlsx"
merging_table2 <- read_excel(merging_table_DC)
head(data.frame(merging_table2))  

# convert to factor with merged clusters in desired order                
merging_table2new_cluster <- factor(merging_table2$new_cluster,         
                                     levels = c("DC2 B", "DC2 A", "mregDCs", "DC2 C", "DC1 A", "DC1 B", "CD40_high DC2"))        

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
clust_color4 <- c("#00994D", "#ff80ff", "#8800CC", "#FF3254", "#CCCC00", "#33BAFE", "#1C1C1C")

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

#heatmap of merged NK cluster markers
plotExprHeatmap(daf_DCs_merg1, 
                bin_anno = FALSE, 
                scale = "last",
                k = "DC_clusters",
                by = "cluster_id", 
                k_pal = clust_color4,
                features = c("SLAMF7", "CD38", "PD-L1", "CD11b", "Ly6C", "AF", "CD4", "F4/80", "MHC-II", "B220", "CCR2", "CCR7", "CD40", "Fascin", "Lyve-1"),
                hm_pal = brewer.pal(9, "Greys")
)
daf_DCs_merg1 <- readRDS(file = "./daf_DC_merg1.rds")

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


#plot abundances as boxplot
obj3 <- plotAbundances(daf_DCs_merg1, 
               k = "DC_clusters", 
               by = "cluster_id",
               group_by = "response",
               shape = NULL
)
obj3 + scale_color_manual(values = c("#0066ff", "#339966")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA)


##To analyze now we can generate a GLMM
eil <- metadata(daf_DCs_merg1)$experiment_info 
(da_formula13 <- createFormula(eil,                     
                              cols_fixed = "response",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res13 <- diffcyt(daf_DCs_merg1,                                            
                   formula = da_formula13, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "DC_clusters", verbose = FALSE)  

#examine output
cluster_stats_DCs_batch1 <- rowData(da_res13$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res13$res)$p_adj < FDR_cutoff)

#export stat data
write.xlsx(cluster_stats_DCs_batch1, "~/B16_TME_profiling_batch1_analys/cluster_stats_DCs.xlsx")

##--------------------------------------------------------------------------------------##






##------------------Now subset TAMs-----------------------------------------------------##

#load main analysis object back in
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")

#subset TAMS cells
daf_TAM <- filterSCE(daf_X_clean, cluster_id == c("TAMs"), k = "merging1")

levels(daf_TAM$condition)


#Make MDS plot  
pbMDS(daf_TAM, by = "sample_id", color_by = "condition", shape_by = "response", size_by = TRUE)

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
               features = c("CD45", "CD11b", "CD4", "F4/80", "Ly6C", "MHC-II", "CCR2", "CCR7", "CD40", "CD38", "CD206", "Lyve-1", "SLAMF7", "PD-L1") 
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
                                     levels = c("TAM A", "TAM B", "TAM C", "TAM D", "TAM E", "TAM F", "TAM G", "CD206_high TAMs","SLAMF7_high CD38_high TAMs"))        

daf_TAM <- mergeClusters(daf_TAM, k = "meta10",                                  
                       table = merging_table4, 
                       id = "TAM_clusters",
                       overwrite = TRUE
)   


#save objects
saveRDS(daf_TAM, file = "./daf_TAM.rds")

#change cluster colors
clust_color6 <- c("#FFAA00", "#AA00FF", "#555555", "#FF65FD", "#66CBFD", "#00994D", "#FF6666", "#3332FE", "#ff3333")

#UMAP TAM
plotDR(daf_TAM, "UMAP", color_by = "TAM_clusters", k_pal = clust_color6, a_pal = rev(hcl.colors(10, "Oranges"))) +
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
  bins = 100,
  k_pal = clust_color6,
  assay = "exprs",
  label = c("target"),
  zeros = TRUE
) + geom_point(size=0.5)


##  CALCULATE CORRELATION COEFFICIENT B/W SF7 AND PD-L1 ON ALL TAMS   ##
library(corrplot)
SF7_pdl1_expr <- data.frame(colData(daf_TAM), t(assay(daf_TAM, "exprs")), check.names = FALSE)

#  calculate correlation coefficients 
cor_matrix_1 = cor(SF7_pdl1_expr[ ,c(14,19)], method='pearson',use='pairwise.complete.obs')
# Now produce the plot
all_cor_matrix1 <- corrplot(cor_matrix_1, method='circle', type='upper', order="hclust")
all_cor_matrix1

#------------------compare TAM cluster frequencies by condition-----------------------------#

#plot abundances as barplot
plotAbundances(daf_TAM, k = "TAM_clusters", by = "sample_id", k_pal = clust_color6)


#plot abundances as boxplot
plotAbundances(daf_TAM, 
               k = "TAM_clusters", 
               by = "cluster_id",
               group_by = "condition",
               k_pal = clust_color6,
               shape = NULL
)


##To analyze now we can generate a GLMM
eiw <- metadata(daf_TAM)$experiment_info 
(da_formula30 <- createFormula(eiw,                     
                               cols_fixed = "condition",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res30 <- diffcyt(daf_TAM,                                            
                    formula = da_formula30, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "TAM_clusters", verbose = FALSE)  

#examine output
cluster_stats30 <- rowData(da_res30$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res30$res)$p_adj < FDR_cutoff)

#export stat data
write.xlsx(cluster_stats30, "~/B16_TME_profiling_batch1_analys/cluster_stats_TAMs.xlsx")

##--------------------------------------------------------------------------------------##

##  CALCULATE CORRELATION COEFFICIENT B/W SF7 AND PD-L1 ON ALL cell subsets   ##
library(ggplot2)
#add cluster names to metadata
daf_X_clean$merg_clust <- cluster_ids(daf_X_clean, "merging1")
#make df
SF7_pdl1_expr <- data.frame(colData(daf_X_clean), t(assay(daf_X_clean, "exprs")), check.names = FALSE)
#remove unnecessary columns
SF7_pdl1_expr_clean <- SF7_pdl1_expr %>%
  select(., "merg_clust", "SLAMF7", "PD-L1")
#rename PD-L1
SF7_pdl1_expr_clean$`PD-L1` -> SF7_pdl1_expr_clean$PDL1
SF7_pdl1_expr_clean <- select(SF7_pdl1_expr_clean, "merg_clust", "PDL1", "SLAMF7")

#  calculate correlation coefficients 
correlate1 <- SF7_pdl1_expr_clean %>%
  group_by(merg_clust) %>% 
  summarise(r = cor(SLAMF7, PDL1, method='spearman', use='pairwise.complete.obs'))

#plot it
cor_plot <- ggplot(data=correlate1, aes(x=reorder(merg_clust, r), y=r)) +
  geom_bar(stat="identity", fill = "black") +
  coord_flip() 
cor_plot + theme(axis.line = element_line(linetype = "solid"), 
    axis.ticks = element_line(colour = "black", 
        size = 0.7), axis.title = element_text(size = 21), 
    axis.text = element_text(size = 14, colour = "black"), 
    plot.title = element_text(size = 21), 
    panel.background = element_rect(fill = NA)) +labs(x = NULL, y = "Spearmans R") 
##--------------------------------------------------------------------------------------##


# comparing cluster frequencies b/w ad_null from batch 1 and 2 ------------

#from batch 1  Gets you number of cells per cluster
batch1_freq <- data.frame(colData(daf_X_clean), check.names = FALSE)
batch1_freq_clean <- batch1_freq %>%
  group_by(sample_id) %>%
  select(merg_clust, condition) %>%
  filter(condition == "Ad_Null") 

batch1_sub_freq <- batch1_freq_clean %>% 
  group_by(sample_id, merg_clust) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame()
  
#now get total number of cells per sample for batch 1
tot_cells <- batch1_sub_freq %>%
group_by(sample_id) %>%
  summarise(tot_cells = sum(n)) %>%
  as.data.frame()

#now calculate frequency from above 2 dfs
merg_df <- left_join(batch1_sub_freq, tot_cells, by = "sample_id")
merg_df <- merg_df %>%
  mutate(freq = (n/tot_cells)*100) %>%
  mutate(exp = 1) %>%
  select("freq", "exp", "merg_clust")


#now batch 2
daf_X_clean_batch2 <- readRDS(file = "./daf_X_clean_batch2.rds")
#add cluster names to metadata
daf_X_clean_batch2$merg_clust <- cluster_ids(daf_X_clean_batch2, "merging1")

#from batch 2  Gets you number of cells per cluster
batch2_freq <- data.frame(colData(daf_X_clean_batch2), check.names = FALSE)
batch2_freq_clean <- batch2_freq %>%
  group_by(sample_id) %>%
  select(merg_clust, condition, genotype) %>%
  filter(condition == "Ad_null") %>%
  filter(genotype == "WT")

batch2_sub_freq <- batch2_freq_clean %>% 
  group_by(sample_id, merg_clust) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame()

#now get total number of cells per sample for batch 2
tot_cells_b2 <- batch2_sub_freq %>%
  group_by(sample_id) %>%
  summarise(tot_cells = sum(n)) %>%
  as.data.frame()

#now calculate frequency from above 2 dfs
merg_df2 <- left_join(batch2_sub_freq, tot_cells_b2, by = "sample_id")
merg_df2 <- merg_df2 %>%
  mutate(freq = (n/tot_cells)*100) %>%
  mutate(exp = 2) %>%
  select("freq", "exp", "merg_clust")


#combine 2 DF's w/ frequencies
#first clean up mismatched names b/w the 2 experiments
merg_df$merg_clust <- gsub("Pro-B cell derived Marcophages","Pro-B cell derived Macrophages", merg_df$merg_clust)

#merge DFs
merg_comb2 <- bind_rows(merg_df, merg_df2)
#remove Unknown Myeloid and CCR2_high mono since they are not in both experiments
merg_comb2 <- merg_comb2 %>%
  filter(merg_clust != "Unknown myeloid") %>%
  filter(merg_clust != "CCR2_high Monocytes")

#rename exp varaible to be more descriptive
merg_comb2$exp <- gsub("1","Control", merg_comb2$exp)
merg_comb2$exp <- gsub("2","Acute Ad-Null Administration", merg_comb2$exp)

#now plot as grouped boxplot
group_comp <- ggplot(merg_comb2, aes(x=merg_clust, y=freq, fill=exp)) +
  geom_boxplot(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize = 0.5) + 
  theme(axis.line = element_line(linetype = "solid"), 
                 axis.ticks = element_line(colour = "black", 
                                           size = 0.7), axis.title = element_text(size = 17), 
                 axis.text = element_text(size = 10, colour = "black"), 
                 plot.title = element_text(size = 14), 
                 panel.background = element_rect(fill = NA), 
        axis.text.x = element_text(vjust = 0.5, 
                                   angle = 45)) +labs(x = NULL, y = "Frequency of immune cells") +
  scale_fill_manual(values=c("#737373", "#ffffff"))

group_comp

group_comp + theme(axis.title = element_text(size = 17), 
    axis.text = element_text(size = 10), 
    axis.text.x = element_text(vjust = 0.5, 
        angle = 45))


##-------For Batch 1---------##
#export df of frequencies of all cell types per sample
#and a second df of median expression of all markers per cluster, per sample

#just edit this to be for Ad-SF7-Fc
# comparing cluster frequencies b/w ad_SF7-Fc from batch 1 and 2 ------------

#from batch 1  Gets you number of cells per cluster
batch1_freq <- data.frame(colData(daf_X_clean), check.names = FALSE)
batch1_freq_clean <- batch1_freq %>%
  dplyr::group_by(sample_id) %>%
  dplyr::select(merg_clust, condition) %>%
  dplyr::filter(condition == "Ad_SLAMF7_Fc") 

batch1_sub_freq <- batch1_freq_clean %>% 
  group_by(sample_id, merg_clust) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame()

#now get total number of cells per sample for batch 1
tot_cells <- batch1_sub_freq %>%
  group_by(sample_id) %>%
  summarise(tot_cells = sum(n)) %>%
  as.data.frame()

#now calculate frequency from above 2 dfs
merg_df <- left_join(batch1_sub_freq, tot_cells, by = "sample_id")
merg_df <- merg_df %>%
  mutate(freq = (n/tot_cells)*100) %>%
  mutate(exp = 1) 

merg_df$merg_clust <- gsub("Pro-B cell derived Marcophages","Pro-B cell derived Macrophages", merg_df$merg_clust)
write.xlsx(merg_df, "~/B16_TME_profiling_batch1_analys_repeat/batch1_clust_freq_iEN.xlsx")

#now batch 2
daf_X_clean_batch2 <- readRDS(file = "./daf_X_clean_batch2.rds")
#add cluster names to metadata
daf_X_clean_batch2$merg_clust <- cluster_ids(daf_X_clean_batch2, "merging1")

#from batch 2  Gets you number of cells per cluster
batch2_freq <- data.frame(colData(daf_X_clean_batch2), check.names = FALSE)
batch2_freq_clean <- batch2_freq %>%
  group_by(sample_id) %>%
  select(merg_clust, condition, genotype) %>%
  dplyr::filter(condition == "Ad_SLAMF7_Fc") %>%
  dplyr::filter(genotype == "WT")

batch2_sub_freq <- batch2_freq_clean %>% 
  group_by(sample_id, merg_clust) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame()

#now get total number of cells per sample for batch 2
tot_cells_b2 <- batch2_sub_freq %>%
  group_by(sample_id) %>%
  summarise(tot_cells = sum(n)) %>%
  as.data.frame()

#now calculate frequency from above 2 dfs
merg_df2 <- left_join(batch2_sub_freq, tot_cells_b2, by = "sample_id")
merg_df2 <- merg_df2 %>%
  mutate(freq = (n/tot_cells)*100) %>%
  mutate(exp = 2) 

write.xlsx(merg_df2, "~/B16_TME_profiling_batch1_analys_repeat/batch2_clust_freq_iEN.xlsx")




##-------now get median marker expression on each cell type for batch 1 and 2
setwd("~/B16_TME_profiling_batch1_analys_repeat")
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")

#add cluster names to metadata
daf_X_clean$merg_clust <- cluster_ids(daf_X_clean, "merging1")
#make df
batch1_freq <- data.frame(colData(daf_X_clean), t(assay(daf_X_clean, "exprs")), check.names = FALSE)
batch1_freq$merg_clust <- gsub("Pro-B cell derived Marcophages","Pro-B cell derived Macrophages", batch1_freq$merg_clust)

#remove unnecessary columns and non-Ad-SF7-Fc
batch1_freq <- batch1_freq %>%
  dplyr::filter(condition == "Ad_SLAMF7_Fc") %>%
  select(., !one_of("file_name", "sample_id.1", "cluster_id", "response", "condition")) 

#get median expression of marker
batch1_freq_med <- batch1_freq %>%
  group_by(sample_id, merg_clust) %>%
  summarise_all(median) %>%
  as.data.frame()
  
write.xlsx(batch1_freq_med, "~/B16_TME_profiling_batch1_analys_repeat/batch1_clust_expr_iEN.xlsx")




#now batch 2
daf_X_clean_batch2 <- readRDS(file = "./daf_X_clean_batch2.rds")
#add cluster names to metadata
daf_X_clean_batch2$merg_clust <- cluster_ids(daf_X_clean_batch2, "merging1")

#make df
batch2_freq <- data.frame(colData(daf_X_clean_batch2), t(assay(daf_X_clean_batch2, "exprs")), check.names = FALSE)

#remove unnecessary columns and non-Ad-SF7-Fc
batch2_freq <- batch2_freq %>%
  dplyr::filter(condition == "Ad_SLAMF7_Fc") %>%
  select(., !one_of("file_name", "sample_id.1", "cluster_id", "response", "condition", "genotype")) 

#get median expression of marker
batch2_freq_med <- batch2_freq %>%
  group_by(sample_id, merg_clust) %>%
  summarise_all(median) %>%
  as.data.frame()

write.xlsx(batch2_freq_med, "~/B16_TME_profiling_batch1_analys_repeat/batch2_clust_expr_iEN.xlsx")

