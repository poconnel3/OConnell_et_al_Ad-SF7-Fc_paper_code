## script for predicting Ad-SF7-Fc activation or inhibition
#3/18/2021


#A higher SLAMF7_score means a cell looks more like a WT cell and less like a SF7-KO cell (since I subtracted SF7-KO from WT)

#-----------------------------

#let's start w/ NK cells

#Packages
library(CATALYST)
library(tidyverse)
library(readxl)
library(ggplot2)
#library(RColorBrewer)
#library(xlsx)


setwd("~/B16_TME_profiling_batch2_analys")

#load the data 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")
daf_X_NK <-  filterSCE(daf_X_clean, cluster_id == "NK cells", k = "merging1")

#---------------Get marker differences b/w KO (Ad-null) and WT (Ad-null)
#Subset to only KO
daf_X_NK_KO <-  filterSCE(daf_X_NK, genotype == "SLAMF7_KO")
levels(daf_X_NK_KO$genotype)

KO_NK_med_mtx <- data.frame(colData(daf_X_NK_KO), t(assay(daf_X_NK_KO, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
KO_NK_med_mtx <- KO_NK_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_NK_KO <- KO_NK_med_mtx %>%
  summarise_all(median)

#now WT
#Subset to only WT
daf_X_NK_WT <-  filterSCE(daf_X_NK, genotype == "WT")
levels(daf_X_NK_WT$genotype)

WT_NK_med_mtx <- data.frame(colData(daf_X_NK_WT), t(assay(daf_X_NK_WT, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
WT_NK_med_mtx <- WT_NK_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_NK_WT <- WT_NK_med_mtx %>%
  summarise_all(median)

#I want the signature of what cell with SLAMF7 signaling looks like so I subtract KO values from WT values
mark_diff <- t(mark_stats_NK_WT-mark_stats_NK_KO) %>%
  as.data.frame() %>%
  arrange(., desc(V1)) %>%
  rownames_to_column(., "markers")

#from a visual inspection it looks like I should be able to use the top 4 and bottom 4 markers as a signature. 

##---------Get NK cell KO marker signature---------------
high <- mark_diff[1:4,] 
low <- mark_diff[22:25,] 

SF7_KO_sig_NK <- bind_rows("up" = high, "down" = low, .id = "sign")




##---------Sum marker expression of NK sig markers in Ad-SF7-Fc NK cells and generate score---------

#subset to only Ad-Sf7-Fc group
daf_X_NK_SF7_Fc <-  filterSCE(daf_X_NK, condition == "Ad_SLAMF7_Fc")
levels(daf_X_NK_SF7_Fc$condition)

Ad_SF7_exprs_NK <- data.frame(colData(daf_X_NK_SF7_Fc), t(assay(daf_X_NK_SF7_Fc, "exprs")), check.names = FALSE)
#clean up
Ad_SF7_exprs_NK <- Ad_SF7_exprs_NK %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#subset out high sig markers and get pan-high marker score
SF7_Fc_high <- Ad_SF7_exprs_NK %>%
  select(one_of("CD11c", "NK1.1", "CD11b", "Ly6C")) %>%
  transmute(sum = rowSums(across(where(is.numeric))))

#subset out low sig markers and get pan-low marker score
SF7_Fc_low <- Ad_SF7_exprs_NK %>%
  select(one_of("CD90", "B220", "CD38", "PD-L1")) %>%
  transmute(sum = rowSums(across(where(is.numeric))))

#subtract low values from high and normalize by total marker number (on a per cell basis)
high_low_dif <- (SF7_Fc_high - SF7_Fc_low)/8 

high_low_dif <- high_low_dif %>%
  rownames_to_column('rn')

#rename score
names(high_low_dif)[2] <- "SLAMF7_score"





##---------Join SLAMF7_score back to main expression df---------

#Add SLAMF7 scores back to main expression df
Ad_SF7_fc_scored_NK <- data.frame(colData(daf_X_NK_SF7_Fc), t(assay(daf_X_NK_SF7_Fc, "exprs")), check.names = FALSE) %>%
  rownames_to_column('rn')

Ad_SF7_fc_scored_NK <- left_join(Ad_SF7_fc_scored_NK, high_low_dif, by = "rn")





##-------------------Now plot it!----------------

#now make scatter plot of SLAMF7_score by SLAMF7 expression
score_plot_NK <- ggplot(Ad_SF7_fc_scored_NK, aes(x=SLAMF7_score, y=SLAMF7)) +
  geom_point(size=0.6, alpha = 0.5) +
  scale_fill_manual(values = c("#b3b3b3"))
score_plot_NK + theme(axis.line = element_line(linetype = "solid"), 
                   axis.ticks = element_line(size = 0.9), 
                   axis.title = element_text(size = 15), 
                   axis.text = element_text(size = 14, colour = "black"), 
                   axis.text.x = element_text(size = 14), 
                   panel.background = element_rect(fill = NA))


#subset out high sig markers and get pan-high marker score
SF7_KO_high <- KO_NK_med_mtx %>%
  select(one_of("CD11c", "NK1.1", "CD11b", "Ly6C")) %>%
  transmute(sum = rowSums(across(where(is.numeric))))

#subset out low sig markers and get pan-low marker score
SF7_KO_low <- KO_NK_med_mtx %>%
  select(one_of("CD90", "B220", "CD38", "PD-L1")) %>%
  transmute(sum = rowSums(across(where(is.numeric))))

#subtract low values from high and normalize by total marker number (on a per cell basis)
high_low_dif_KO <- (SF7_KO_high - SF7_KO_low)/8 

high_low_dif_KO <- high_low_dif_KO %>%
  rownames_to_column('rn')

#rename score
names(high_low_dif_KO)[2] <- "SLAMF7_score"


#Add SLAMF7 scores back to main expression df

KO_NK_med_mtx <- data.frame(colData(daf_X_NK_KO), t(assay(daf_X_NK_KO, "exprs")), check.names = FALSE) %>%
  rownames_to_column('rn')

SF7_KO_scored <- left_join(KO_NK_med_mtx, high_low_dif_KO, by = "rn")

#merge SF7_KO_scored to Ad_SF7_Fc_scored
KO_and_Fc_scored <- bind_rows(SF7_KO_scored, Ad_SF7_fc_scored_NK)
#check it
table(KO_and_Fc_scored$condition, KO_and_Fc_scored$genotype)

#Now plot  
score_plot <- ggplot(KO_and_Fc_scored, aes(x=SLAMF7_score, y=SLAMF7, color = genotype )) +
  geom_point(size=0.5, alpha = 0.5) +
  scale_fill_manual(values = c("#b3b3b3"))
score_plot + theme(axis.line = element_line(linetype = "solid"), 
                   axis.ticks = element_line(size = 0.9), 
                   axis.title = element_text(size = 15), 
                   axis.text = element_text(size = 14, colour = "black"), 
                   axis.text.x = element_text(size = 14), 
                   panel.background = element_rect(fill = NA))




##---------now trying w/ Dr. Carlson's approach-------------##

#subset to pertinant data
diag_calc <- Ad_SF7_fc_scored_NK %>%
  select(one_of("SLAMF7", "SLAMF7_score")) %>%
  as.data.frame(.) 

#get quantiles
SF7q <- quantile(diag_calc$SLAMF7)
SF7_scoreq <- quantile(diag_calc$SLAMF7_score)

#Now plot the data and identify a line that passes through the lower left quantile intersection and the 
#upper right quantile intersection. Then use the slope to identify parallel lines that pass through 
#the upper left and lower right intersections:
plot(diag_calc$SLAMF7~diag_calc$SLAMF7_score, diag_calc, pch=20)
abline(v=SF7_scoreq[2:4], lty=3)
abline(h=SF7q[2:4], lty=3)
diag <- lm(SF7q[c(2, 4)]~SF7_scoreq[c(2, 4)])
points(SF7_scoreq[c(2, 4)], SF7q[c(2, 4)], cex=2, col="red", lwd=2)
abline(diag)
b <- coef(diag)[2]
a1 <- SF7q[4] - b * SF7_scoreq[2]
a2 <- SF7q[2] - b * SF7_scoreq[4]
abline(a1, b)
abline(a2, b)

#Now identify all points above and below the 2 diagonal lines
res1 <- diag_calc$SLAMF7 - (a1 + b * diag_calc$SLAMF7_score)
res2 <- (a2 + b * diag_calc$SLAMF7_score) - diag_calc$SLAMF7
clr <- c("black", "purple", "darkorange")
idx <- ifelse(res1 > 0, 3, ifelse(res2 > 0, 2, 1))
plot(diag_calc$SLAMF7~diag_calc$SLAMF7_score, pch=20, col=clr[idx])
abline(a1, b, col="red")
abline(a2, b, col="red")

#Get identification of outlier points
position <- c("neither", "activated", "blocked")
diag_calc$outlier <- position[idx]
diag_calc <- diag_calc %>% rownames_to_column('rn')
#remove SF7 and SF7_score so no redundancy when I join
diag_calc <- diag_calc %>%
  select(., rn, outlier)
head(diag_calc)

#Join calls of outlier points back to main df
Ad_SF7_fc_scored_NK_called <- left_join(Ad_SF7_fc_scored_NK, diag_calc, by = "rn")

#Now plot  (nicely) (WORKS BEAUTIFULY)
score_plot2 <- ggplot(Ad_SF7_fc_scored_NK_called, aes(x=SLAMF7_score, y=SLAMF7, color = outlier, group = outlier)) +
  geom_point(aes(size = outlier), alpha = 0.9) +
  scale_size_manual(values=c(3,3,1.2)) +
  scale_color_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))
score_plot2 + theme(axis.line = element_line(linetype = "solid"), 
                   axis.ticks = element_line(size = 0.9), 
                   axis.title = element_text(size = 15), 
                   axis.text = element_text(size = 14, colour = "black"), 
                   axis.text.x = element_text(size = 14), 
                   panel.background = element_rect(fill = NA))


#make pie chart of proportion of outlier calls
NK_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_NK_called$outlier,
  condition = Ad_SF7_fc_scored_NK_called$condition) %>%
  as.data.frame(.) 

ggplot(data=NK_pred_pie, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))

NK_pred_pie
sum(NK_pred_pie$Freq)




##-----------------Now repeat for DCs--------------------##

#load the data 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")
daf_X_DC <-  filterSCE(daf_X_clean, cluster_id == "DCs", k = "merging1")

#---------------Get marker differences b/w KO (Ad-null) and WT (Ad-null)
#Subset to only KO
daf_X_DC_KO <-  filterSCE(daf_X_DC, genotype == "SLAMF7_KO")
levels(daf_X_DC_KO$genotype)

KO_DC_med_mtx <- data.frame(colData(daf_X_DC_KO), t(assay(daf_X_DC_KO, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
KO_DC_med_mtx <- KO_DC_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_DC_KO <- KO_DC_med_mtx %>%
  summarise_all(median)

#now WT
#Subset to only WT
daf_X_DC_WT <-  filterSCE(daf_X_DC, genotype == "WT")
levels(daf_X_DC_WT$genotype)

WT_DC_med_mtx <- data.frame(colData(daf_X_DC_WT), t(assay(daf_X_DC_WT, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
WT_DC_med_mtx <- WT_DC_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_DC_WT <- WT_DC_med_mtx %>%
  summarise_all(median)

#I want the signature of what cell with SLAMF7 signaling looks like so I subtract KO values from WT values
mark_diff <- t(mark_stats_DC_WT-mark_stats_DC_KO) %>%
  as.data.frame() %>%
  arrange(., desc(V1)) %>%
  rownames_to_column(., "markers")

#from a visual inspection it looks like I should be able to use the top 4 and bottom 4 markers as a signature. 

##---------Get NK cell KO marker signature---------------
high <- mark_diff[1:4,] 
low <- mark_diff[22:25,] 

SF7_KO_sig_DC <- bind_rows("up" = high, "down" = low, .id = "sign")




##---------Sum marker expression of DC sig markers in Ad-SF7-Fc NK cells and generate score---------

#subset to only Ad-Sf7-Fc group
daf_X_DC_SF7_Fc <-  filterSCE(daf_X_DC, condition == "Ad_SLAMF7_Fc")
levels(daf_X_DC_SF7_Fc$condition)

Ad_SF7_exprs_DC <- data.frame(colData(daf_X_DC_SF7_Fc), t(assay(daf_X_DC_SF7_Fc, "exprs")), check.names = FALSE)
#clean up
Ad_SF7_exprs_DC <- Ad_SF7_exprs_DC %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#subset out high sig markers and get pan-high marker score
SF7_Fc_high <- Ad_SF7_exprs_DC %>%
  select(one_of("CD11c", "PD-L1", "AF", "Ly6C")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subset out low sig markers and get pan-low marker score
SF7_Fc_low <- Ad_SF7_exprs_DC %>%
  select(one_of("Lyve-1", "CD4", "CD38", "Fascin")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subtract low values from high and normalize by total marker number (on a per cell basis)
high_low_dif <- (SF7_Fc_high - SF7_Fc_low)/8 

high_low_dif <- high_low_dif %>%
  rownames_to_column('rn')

#rename score
names(high_low_dif)[2] <- "SLAMF7_score"





##---------Join SLAMF7_score back to main expression df---------

#Add SLAMF7 scores back to main expression df
Ad_SF7_fc_scored_DC <- data.frame(colData(daf_X_DC_SF7_Fc), t(assay(daf_X_DC_SF7_Fc, "exprs")), check.names = FALSE) %>%
  rownames_to_column('rn')

Ad_SF7_fc_scored_DC <- left_join(Ad_SF7_fc_scored_DC, high_low_dif, by = "rn")
summary(Ad_SF7_fc_scored_DC$SLAMF7_score)



##---------now trying w/ Dr. Carlson's approach-------------##

#subset to pertinant data
diag_calc <- Ad_SF7_fc_scored_DC %>%
  select(one_of("SLAMF7", "SLAMF7_score")) %>%
  as.data.frame(.) 

#get quantiles
SF7q <- quantile(diag_calc$SLAMF7)
SF7_scoreq <- quantile(diag_calc$SLAMF7_score)

#Now plot the data and identify a line that passes through the lower left quantile intersection and the 
#upper right quantile intersection. Then use the slope to identify parallel lines that pass through 
#the upper left and lower right intersections:
diag <- lm(SF7q[c(2, 4)]~SF7_scoreq[c(2, 4)])
b <- coef(diag)[2]
a1 <- SF7q[4] - b * SF7_scoreq[2]
a2 <- SF7q[2] - b * SF7_scoreq[4]


#Now identify all points above and below the 2 diagonal lines
res1 <- diag_calc$SLAMF7 - (a1 + b * diag_calc$SLAMF7_score)
res2 <- (a2 + b * diag_calc$SLAMF7_score) - diag_calc$SLAMF7
idx <- ifelse(res1 > 0, 3, ifelse(res2 > 0, 2, 1))


#Get identification of outlier points
position <- c("neither", "activated", "blocked")
diag_calc$outlier <- position[idx]
diag_calc <- diag_calc %>% rownames_to_column('rn')
#remove SF7 and SF7_score so no redundancy when I join
diag_calc <- diag_calc %>%
  select(., rn, outlier)
head(diag_calc)

#Join calls of outlier points back to main df
Ad_SF7_fc_scored_DC_called <- left_join(Ad_SF7_fc_scored_DC, diag_calc, by = "rn")

#Now plot  (nicely) (WORKS BEAUTIFULY)
score_plot2 <- ggplot(Ad_SF7_fc_scored_DC_called, aes(x=SLAMF7_score, y=SLAMF7, color = outlier, group = outlier)) +
  geom_point(aes(size = outlier), alpha = 0.9) +
  scale_size_manual(values=c(1.5,1.5,0.9)) +       #change point sizes based on number of cells
  scale_color_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))
score_plot2 + theme(axis.line = element_line(linetype = "solid"), 
                    axis.ticks = element_line(size = 0.9), 
                    axis.title = element_text(size = 15), 
                    axis.text = element_text(size = 14, colour = "black"), 
                    axis.text.x = element_text(size = 14), 
                    panel.background = element_rect(fill = NA))


#make pie chart of proportion of outlier calls
DC_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_DC_called$outlier,
  condition = Ad_SF7_fc_scored_DC_called$condition) %>%
  as.data.frame(.) 

ggplot(data=DC_pred_pie, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))

DC_pred_pie
sum(DC_pred_pie$Freq)





##-----------------Now repeat for TAMs--------------------##

#load the data 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")
daf_X_TAM <-  filterSCE(daf_X_clean, cluster_id == "TAMs", k = "merging1")

#---------------Get marker differences b/w KO (Ad-null) and WT (Ad-null)
#Subset to only KO
daf_X_tam_KO <-  filterSCE(daf_X_TAM, genotype == "SLAMF7_KO")
levels(daf_X_tam_KO$genotype)

KO_tam_med_mtx <- data.frame(colData(daf_X_tam_KO), t(assay(daf_X_tam_KO, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
KO_tam_med_mtx <- KO_tam_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_tam_KO <- KO_tam_med_mtx %>%
  summarise_all(median)

#now WT
#Subset to only WT
daf_X_tam_WT <-  filterSCE(daf_X_TAM, genotype == "WT")
levels(daf_X_tam_WT$genotype)

WT_tam_med_mtx <- data.frame(colData(daf_X_tam_WT), t(assay(daf_X_tam_WT, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
WT_tam_med_mtx <- WT_tam_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_tam_WT <- WT_tam_med_mtx %>%
  summarise_all(median)

#I want the signature of what cell with SLAMF7 signaling looks like so I subtract KO values from WT values
mark_diff <- t(mark_stats_tam_WT-mark_stats_tam_KO) %>%
  as.data.frame() %>%
  arrange(., desc(V1)) %>%
  rownames_to_column(., "markers")


##---------Get tam cell KO marker signature---------------
high <- mark_diff[1:5,]   #I increased the number of markers here since more were strongly changed
low <- mark_diff[22:25,] 

SF7_KO_sig_tam <- bind_rows("up" = high, "down" = low, .id = "sign")




##---------Sum marker expression of tam sig markers in Ad-SF7-Fc tam cells and generate score---------

#subset to only Ad-Sf7-Fc group
daf_X_tam_SF7_Fc <-  filterSCE(daf_X_TAM, condition == "Ad_SLAMF7_Fc")
levels(daf_X_tam_SF7_Fc$condition)

Ad_SF7_exprs_tam <- data.frame(colData(daf_X_tam_SF7_Fc), t(assay(daf_X_tam_SF7_Fc, "exprs")), check.names = FALSE)
#clean up
Ad_SF7_exprs_tam <- Ad_SF7_exprs_tam %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#subset out high sig markers and get pan-high marker score
SF7_Fc_high <- Ad_SF7_exprs_tam %>%
  select(one_of("CD11c", "CD45", "CD11b", "PD-L1", "MHC-II")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subset out low sig markers and get pan-low marker score
SF7_Fc_low <- Ad_SF7_exprs_tam %>%
  select(one_of("CD38", "F4/80", "AF", "CD206")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subtract low values from high and normalize by total marker number (on a per cell basis)
high_low_dif <- (SF7_Fc_high - SF7_Fc_low)/9     #changed denominator since I added a marker 

high_low_dif <- high_low_dif %>%
  rownames_to_column('rn')

#rename score
names(high_low_dif)[2] <- "SLAMF7_score"


##---------Join SLAMF7_score back to main expression df---------

#Add SLAMF7 scores back to main expression df
Ad_SF7_fc_scored_tam <- data.frame(colData(daf_X_tam_SF7_Fc), t(assay(daf_X_tam_SF7_Fc, "exprs")), check.names = FALSE) %>%
  rownames_to_column('rn')

Ad_SF7_fc_scored_tam <- left_join(Ad_SF7_fc_scored_tam, high_low_dif, by = "rn")
summary(Ad_SF7_fc_scored_tam$SLAMF7_score)



##---------now trying w/ Dr. Carlson's approach-------------##
#subset to pertinant data
diag_calc <- Ad_SF7_fc_scored_tam %>%
  select(one_of("SLAMF7", "SLAMF7_score")) %>%
  as.data.frame(.) 

#get quantiles
SF7q <- quantile(diag_calc$SLAMF7)
SF7_scoreq <- quantile(diag_calc$SLAMF7_score)

#Now plot the data and identify a line that passes through the lower left quantile intersection and the 
#upper right quantile intersection. Then use the slope to identify parallel lines that pass through 
#the upper left and lower right intersections:
diag <- lm(SF7q[c(2, 4)]~SF7_scoreq[c(2, 4)])
b <- coef(diag)[2]
a1 <- SF7q[4] - b * SF7_scoreq[2]
a2 <- SF7q[2] - b * SF7_scoreq[4]

#Now identify all points above and below the 2 diagonal lines
res1 <- diag_calc$SLAMF7 - (a1 + b * diag_calc$SLAMF7_score)
res2 <- (a2 + b * diag_calc$SLAMF7_score) - diag_calc$SLAMF7
idx <- ifelse(res1 > 0, 3, ifelse(res2 > 0, 2, 1))


#Get identification of outlier points
position <- c("neither", "activated", "blocked")
diag_calc$outlier <- position[idx]
diag_calc <- diag_calc %>% rownames_to_column('rn')
#remove SF7 and SF7_score so no redundancy when I join
diag_calc <- diag_calc %>%
  select(., rn, outlier)
head(diag_calc)

#Join calls of outlier points back to main df
Ad_SF7_fc_scored_tam_called <- left_join(Ad_SF7_fc_scored_tam, diag_calc, by = "rn")

#Now plot  (nicely) (WORKS BEAUTIFULY)
score_plot2 <- ggplot(Ad_SF7_fc_scored_tam_called, aes(x=SLAMF7_score, y=SLAMF7, color = outlier, group = outlier)) +
  geom_point(aes(size = outlier), alpha = 0.9) +
  scale_size_manual(values=c(1.5,1.5,0.9)) +       #change point sizes based on number of cells
  scale_color_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))
score_plot2 + theme(axis.line = element_line(linetype = "solid"), 
                    axis.ticks = element_line(size = 0.9), 
                    axis.title = element_text(size = 15), 
                    axis.text = element_text(size = 14, colour = "black"), 
                    axis.text.x = element_text(size = 14), 
                    panel.background = element_rect(fill = NA))


#make pie chart of proportion of outlier calls
tam_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_tam_called$outlier,
  condition = Ad_SF7_fc_scored_tam_called$condition) %>%
  as.data.frame(.) 

ggplot(data=tam_pred_pie, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))

tam_pred_pie
sum(tam_pred_pie$Freq)




##-----------------Now repeat for CD8+ T cells--------------------##

#load the data 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")
daf_X_cd8 <-  filterSCE(daf_X_clean, cluster_id == "CD8+ T cells", k = "merging1")

#---------------Get marker differences b/w KO (Ad-null) and WT (Ad-null)
#Subset to only KO
daf_X_cd8_KO <-  filterSCE(daf_X_cd8, genotype == "SLAMF7_KO")
levels(daf_X_cd8_KO$genotype)

KO_cd8_med_mtx <- data.frame(colData(daf_X_cd8_KO), t(assay(daf_X_cd8_KO, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
KO_cd8_med_mtx <- KO_cd8_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_cd8_KO <- KO_cd8_med_mtx %>%
  summarise_all(median)

#now WT
#Subset to only WT
daf_X_cd8_WT <-  filterSCE(daf_X_cd8, genotype == "WT")
levels(daf_X_cd8_WT$genotype)

WT_cd8_med_mtx <- data.frame(colData(daf_X_cd8_WT), t(assay(daf_X_cd8_WT, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
WT_cd8_med_mtx <- WT_cd8_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_cd8_WT <- WT_cd8_med_mtx %>%
  summarise_all(median)

#I want the signature of what cell with SLAMF7 signaling looks like so I subtract KO values from WT values
mark_diff <- t(mark_stats_cd8_WT-mark_stats_cd8_KO) %>%
  as.data.frame() %>%
  arrange(., desc(V1)) %>%
  rownames_to_column(., "markers")


##---------Get tam cell KO marker signature---------------
high <- mark_diff[1:3,]   #only top 3 since I do not want to include IgD
low <- mark_diff[21:25,]  #last 5

SF7_KO_sig_cd8 <- bind_rows("up" = high, "down" = low, .id = "sign")




##---------Sum marker expression of  sig markers in Ad-SF7-Fc cd8 cells and generate score---------

#subset to only Ad-Sf7-Fc group
daf_X_cd8_SF7_Fc <-  filterSCE(daf_X_cd8, condition == "Ad_SLAMF7_Fc")
levels(daf_X_cd8_SF7_Fc$condition)

Ad_SF7_exprs_cd8 <- data.frame(colData(daf_X_cd8_SF7_Fc), t(assay(daf_X_cd8_SF7_Fc, "exprs")), check.names = FALSE)
#clean up
Ad_SF7_exprs_cd8 <- Ad_SF7_exprs_cd8 %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#subset out high sig markers and get pan-high marker score
SF7_Fc_high <- Ad_SF7_exprs_cd8 %>%
  select(one_of("Ly6C", "AF", "CD90")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subset out low sig markers and get pan-low marker score
SF7_Fc_low <- Ad_SF7_exprs_cd8 %>%
  select(one_of("CD38", "CD4", "Fascin", "NK1.1", "CD11c")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subtract low values from high and normalize by total marker number (on a per cell basis)
high_low_dif <- (SF7_Fc_high - SF7_Fc_low)/8    

high_low_dif <- high_low_dif %>%
  rownames_to_column('rn')

#rename score
names(high_low_dif)[2] <- "SLAMF7_score"


##---------Join SLAMF7_score back to main expression df---------

#Add SLAMF7 scores back to main expression df
Ad_SF7_fc_scored_cd8 <- data.frame(colData(daf_X_cd8_SF7_Fc), t(assay(daf_X_cd8_SF7_Fc, "exprs")), check.names = FALSE) %>%
  rownames_to_column('rn')

Ad_SF7_fc_scored_cd8 <- left_join(Ad_SF7_fc_scored_cd8, high_low_dif, by = "rn")
summary(Ad_SF7_fc_scored_cd8$SLAMF7_score)



##---------now trying w/ Dr. Carlson's approach-------------##
#subset to pertinant data
diag_calc <- Ad_SF7_fc_scored_cd8 %>%
  select(one_of("SLAMF7", "SLAMF7_score")) %>%
  as.data.frame(.) 

#get quantiles
SF7q <- quantile(diag_calc$SLAMF7)
SF7_scoreq <- quantile(diag_calc$SLAMF7_score)

#Now plot the data and identify a line that passes through the lower left quantile intersection and the 
#upper right quantile intersection. Then use the slope to identify parallel lines that pass through 
#the upper left and lower right intersections:
diag <- lm(SF7q[c(2, 4)]~SF7_scoreq[c(2, 4)])
b <- coef(diag)[2]
a1 <- SF7q[4] - b * SF7_scoreq[2]
a2 <- SF7q[2] - b * SF7_scoreq[4]

#Now identify all points above and below the 2 diagonal lines
res1 <- diag_calc$SLAMF7 - (a1 + b * diag_calc$SLAMF7_score)
res2 <- (a2 + b * diag_calc$SLAMF7_score) - diag_calc$SLAMF7
idx <- ifelse(res1 > 0, 3, ifelse(res2 > 0, 2, 1))


#Get identification of outlier points
position <- c("neither", "activated", "blocked")
diag_calc$outlier <- position[idx]
diag_calc <- diag_calc %>% rownames_to_column('rn')
#remove SF7 and SF7_score so no redundancy when I join
diag_calc <- diag_calc %>%
  select(., rn, outlier)
head(diag_calc)

#Join calls of outlier points back to main df
Ad_SF7_fc_scored_cd8_called <- left_join(Ad_SF7_fc_scored_cd8, diag_calc, by = "rn")

#Now plot  (nicely) (WORKS BEAUTIFULY)
score_plot2 <- ggplot(Ad_SF7_fc_scored_cd8_called, aes(x=SLAMF7_score, y=SLAMF7, color = outlier, group = outlier)) +
  geom_point(aes(size = outlier), alpha = 0.9) +
  scale_size_manual(values=c(1.5,1.5,0.9)) +       #change point sizes based on number of cells
  scale_color_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))
score_plot2 + theme(axis.line = element_line(linetype = "solid"), 
                    axis.ticks = element_line(size = 0.9), 
                    axis.title = element_text(size = 15), 
                    axis.text = element_text(size = 14, colour = "black"), 
                    axis.text.x = element_text(size = 14), 
                    panel.background = element_rect(fill = NA))


#make pie chart of proportion of outlier calls
cd8_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_cd8_called$outlier,
  condition = Ad_SF7_fc_scored_cd8_called$condition) %>%
  as.data.frame(.) 

ggplot(data=cd8_pred_pie, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))

cd8_pred_pie
sum(cd8_pred_pie$Freq)



##-----------------Now repeat for pdcs--------------------##

#load the data 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")
daf_X_pdc <-  filterSCE(daf_X_clean, cluster_id == "pDCs", k = "merging1")

#---------------Get marker differences b/w KO (Ad-null) and WT (Ad-null)
#Subset to only KO
daf_X_pdc_KO <-  filterSCE(daf_X_pdc, genotype == "SLAMF7_KO")
levels(daf_X_pdc_KO$genotype)

KO_pdc_med_mtx <- data.frame(colData(daf_X_pdc_KO), t(assay(daf_X_pdc_KO, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
KO_pdc_med_mtx <- KO_pdc_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_pdc_KO <- KO_pdc_med_mtx %>%
  summarise_all(median)

#now WT
#Subset to only WT
daf_X_pdc_WT <-  filterSCE(daf_X_pdc, genotype == "WT")
levels(daf_X_pdc_WT$genotype)

WT_pdc_med_mtx <- data.frame(colData(daf_X_pdc_WT), t(assay(daf_X_pdc_WT, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
WT_pdc_med_mtx <- WT_pdc_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_pdc_WT <- WT_pdc_med_mtx %>%
  summarise_all(median)

#I want the signature of what cell with SLAMF7 signaling looks like so I subtract KO values from WT values
mark_diff <- t(mark_stats_pdc_WT-mark_stats_pdc_KO) %>%
  as.data.frame() %>%
  arrange(., desc(V1)) %>%
  rownames_to_column(., "markers")


##---------Get tam cell KO marker signature---------------
high <- mark_diff[1:4,]   #top 4, but remove IgD
low <- mark_diff[22:25,] 

SF7_KO_sig_pdc <- bind_rows("up" = high, "down" = low, .id = "sign") 
SF7_KO_sig_pdc <- SF7_KO_sig_pdc[-3,]



##---------Sum marker expression of  sig markers in Ad-SF7-Fc pDCs cells and generate score---------

#subset to only Ad-Sf7-Fc group
daf_X_pdc_SF7_Fc <-  filterSCE(daf_X_pdc, condition == "Ad_SLAMF7_Fc")
levels(daf_X_pdc_SF7_Fc$condition)

Ad_SF7_exprs_pdc <- data.frame(colData(daf_X_pdc_SF7_Fc), t(assay(daf_X_pdc_SF7_Fc, "exprs")), check.names = FALSE)
#clean up
Ad_SF7_exprs_pdc <- Ad_SF7_exprs_pdc %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#subset out high sig markers and get pan-high marker score
SF7_Fc_high <- Ad_SF7_exprs_pdc %>%
  select(one_of("CD11c", "Ly6C", "AF")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subset out low sig markers and get pan-low marker score
SF7_Fc_low <- Ad_SF7_exprs_pdc %>%
  select(one_of("CD38", "MHC-II", "CD4", "B220")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subtract low values from high and normalize by total marker number (on a per cell basis)
high_low_dif <- (SF7_Fc_high - SF7_Fc_low)/7     #changed denominator  

high_low_dif <- high_low_dif %>%
  rownames_to_column('rn')

#rename score
names(high_low_dif)[2] <- "SLAMF7_score"


##---------Join SLAMF7_score back to main expression df---------

#Add SLAMF7 scores back to main expression df
Ad_SF7_fc_scored_pdc <- data.frame(colData(daf_X_pdc_SF7_Fc), t(assay(daf_X_pdc_SF7_Fc, "exprs")), check.names = FALSE) %>%
  rownames_to_column('rn')

Ad_SF7_fc_scored_pdc <- left_join(Ad_SF7_fc_scored_pdc, high_low_dif, by = "rn")
summary(Ad_SF7_fc_scored_pdc$SLAMF7_score)



##---------now trying w/ Dr. Carlson's approach-------------##
#subset to pertinant data
diag_calc <- Ad_SF7_fc_scored_pdc %>%
  select(one_of("SLAMF7", "SLAMF7_score")) %>%
  as.data.frame(.) 

#get quantiles
SF7q <- quantile(diag_calc$SLAMF7)
SF7_scoreq <- quantile(diag_calc$SLAMF7_score)

#Now plot the data and identify a line that passes through the lower left quantile intersection and the 
#upper right quantile intersection. Then use the slope to identify parallel lines that pass through 
#the upper left and lower right intersections:
diag <- lm(SF7q[c(2, 4)]~SF7_scoreq[c(2, 4)])
b <- coef(diag)[2]
a1 <- SF7q[4] - b * SF7_scoreq[2]
a2 <- SF7q[2] - b * SF7_scoreq[4]

#Now identify all points above and below the 2 diagonal lines
res1 <- diag_calc$SLAMF7 - (a1 + b * diag_calc$SLAMF7_score)
res2 <- (a2 + b * diag_calc$SLAMF7_score) - diag_calc$SLAMF7
idx <- ifelse(res1 > 0, 3, ifelse(res2 > 0, 2, 1))


#Get identification of outlier points
position <- c("neither", "activated", "blocked")
diag_calc$outlier <- position[idx]
diag_calc <- diag_calc %>% rownames_to_column('rn')
#remove SF7 and SF7_score so no redundancy when I join
diag_calc <- diag_calc %>%
  select(., rn, outlier)
head(diag_calc)

#Join calls of outlier points back to main df
Ad_SF7_fc_scored_pdc_called <- left_join(Ad_SF7_fc_scored_pdc, diag_calc, by = "rn")

#Now plot  (nicely) (WORKS BEAUTIFULY)
score_plot2 <- ggplot(Ad_SF7_fc_scored_pdc_called, aes(x=SLAMF7_score, y=SLAMF7, color = outlier, group = outlier)) +
  geom_point(aes(size = outlier), alpha = 0.9) +
  scale_size_manual(values=c(2,2,1.2)) +       #change point sizes based on number of cells
  scale_color_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))
score_plot2 + theme(axis.line = element_line(linetype = "solid"), 
                    axis.ticks = element_line(size = 0.9), 
                    axis.title = element_text(size = 15), 
                    axis.text = element_text(size = 14, colour = "black"), 
                    axis.text.x = element_text(size = 14), 
                    panel.background = element_rect(fill = NA))


#make pie chart of proportion of outlier calls
pdc_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_pdc_called$outlier,
  condition = Ad_SF7_fc_scored_pdc_called$condition) %>%
  as.data.frame(.) 

ggplot(data=pdc_pred_pie, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))

pdc_pred_pie
sum(pdc_pred_pie$Freq)




##-----------------Now repeat for Monocytes--------------------##

#load the data 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")
daf_X_mono <-  filterSCE(daf_X_clean, cluster_id == "Monocytes", k = "merging1")

#---------------Get marker differences b/w KO (Ad-null) and WT (Ad-null)
#Subset to only KO
daf_X_mono_KO <-  filterSCE(daf_X_mono, genotype == "SLAMF7_KO")
levels(daf_X_mono_KO$genotype)

KO_mono_med_mtx <- data.frame(colData(daf_X_mono_KO), t(assay(daf_X_mono_KO, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
KO_mono_med_mtx <- KO_mono_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_mono_KO <- KO_mono_med_mtx %>%
  summarise_all(median)

#now WT
#Subset to only WT
daf_X_mono_WT <-  filterSCE(daf_X_mono, genotype == "WT")
levels(daf_X_mono_WT$genotype)

WT_mono_med_mtx <- data.frame(colData(daf_X_mono_WT), t(assay(daf_X_mono_WT, "exprs")), check.names = FALSE)
#remove 1st 7 columns of metadata and SF7
WT_mono_med_mtx <- WT_mono_med_mtx %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#Get median expression values for all markers
mark_stats_mono_WT <- WT_mono_med_mtx %>%
  summarise_all(median)

#I want the signature of what cell with SLAMF7 signaling looks like so I subtract KO values from WT values
mark_diff <- t(mark_stats_mono_WT-mark_stats_mono_KO) %>%
  as.data.frame() %>%
  arrange(., desc(V1)) %>%
  rownames_to_column(., "markers")


##---------Get tam cell KO marker signature---------------
high <- mark_diff[1:5,]   #top 5, but remove IgD
low <- mark_diff[21:25,] 

SF7_KO_sig_mono <- bind_rows("up" = high, "down" = low, .id = "sign") 
SF7_KO_sig_mono <- SF7_KO_sig_mono[-3,]



##---------Sum marker expression of  sig markers in Ad-SF7-Fc monos cells and generate score---------

#subset to only Ad-Sf7-Fc group
daf_X_mono_SF7_Fc <-  filterSCE(daf_X_mono, condition == "Ad_SLAMF7_Fc")
levels(daf_X_mono_SF7_Fc$condition)

Ad_SF7_exprs_mono <- data.frame(colData(daf_X_mono_SF7_Fc), t(assay(daf_X_mono_SF7_Fc, "exprs")), check.names = FALSE)
#clean up
Ad_SF7_exprs_mono <- Ad_SF7_exprs_mono %>%
  select( -c("sample_id", "genotype", "response", "condition", "file_name", "sample_id.1", "cluster_id", "SLAMF7"))

#subset out high sig markers and get pan-high marker score
SF7_Fc_high <- Ad_SF7_exprs_mono %>%
  select(one_of("CD11c", "Ly6C", "AF")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subset out low sig markers and get pan-low marker score
SF7_Fc_low <- Ad_SF7_exprs_mono %>%
  select(one_of("CD38", "MHC-II", "CD4", "B220")) %>%   #must edit here each time
  transmute(sum = rowSums(across(where(is.numeric))))

#subtract low values from high and normalize by total marker number (on a per cell basis)
high_low_dif <- (SF7_Fc_high - SF7_Fc_low)/9     #changed denominator  

high_low_dif <- high_low_dif %>%
  rownames_to_column('rn')

#rename score
names(high_low_dif)[2] <- "SLAMF7_score"


##---------Join SLAMF7_score back to main expression df---------

#Add SLAMF7 scores back to main expression df
Ad_SF7_fc_scored_mono <- data.frame(colData(daf_X_mono_SF7_Fc), t(assay(daf_X_mono_SF7_Fc, "exprs")), check.names = FALSE) %>%
  rownames_to_column('rn')

Ad_SF7_fc_scored_mono <- left_join(Ad_SF7_fc_scored_mono, high_low_dif, by = "rn")
summary(Ad_SF7_fc_scored_mono$SLAMF7_score)



##---------now trying w/ Dr. Carlson's approach-------------##
#subset to pertinant data
diag_calc <- Ad_SF7_fc_scored_mono %>%
  select(one_of("SLAMF7", "SLAMF7_score")) %>%
  as.data.frame(.) 

#get quantiles
SF7q <- quantile(diag_calc$SLAMF7)
SF7_scoreq <- quantile(diag_calc$SLAMF7_score)

#Now plot the data and identify a line that passes through the lower left quantile intersection and the 
#upper right quantile intersection. Then use the slope to identify parallel lines that pass through 
#the upper left and lower right intersections:
diag <- lm(SF7q[c(2, 4)]~SF7_scoreq[c(2, 4)])
b <- coef(diag)[2]
a1 <- SF7q[4] - b * SF7_scoreq[2]
a2 <- SF7q[2] - b * SF7_scoreq[4]

#Now identify all points above and below the 2 diagonal lines
res1 <- diag_calc$SLAMF7 - (a1 + b * diag_calc$SLAMF7_score)
res2 <- (a2 + b * diag_calc$SLAMF7_score) - diag_calc$SLAMF7
idx <- ifelse(res1 > 0, 3, ifelse(res2 > 0, 2, 1))


#Get identification of outlier points
position <- c("neither", "activated", "blocked")
diag_calc$outlier <- position[idx]
diag_calc <- diag_calc %>% rownames_to_column('rn')
#remove SF7 and SF7_score so no redundancy when I join
diag_calc <- diag_calc %>%
  select(., rn, outlier)
head(diag_calc)

#Join calls of outlier points back to main df
Ad_SF7_fc_scored_mono_called <- left_join(Ad_SF7_fc_scored_mono, diag_calc, by = "rn")

#Now plot  (nicely) (WORKS BEAUTIFULY)
score_plot2 <- ggplot(Ad_SF7_fc_scored_mono_called, aes(x=SLAMF7_score, y=SLAMF7, color = outlier, group = outlier)) +
  geom_point(aes(size = outlier), alpha = 0.9) +
  scale_size_manual(values=c(1.5,1.5,0.9)) +       #change point sizes based on number of cells
  scale_color_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))
score_plot2 + theme(axis.line = element_line(linetype = "solid"), 
                    axis.ticks = element_line(size = 0.9), 
                    axis.title = element_text(size = 15), 
                    axis.text = element_text(size = 14, colour = "black"), 
                    axis.text.x = element_text(size = 14), 
                    panel.background = element_rect(fill = NA))


#make pie chart of proportion of outlier calls
mono_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_mono_called$outlier,
  condition = Ad_SF7_fc_scored_mono_called$condition) %>%
  as.data.frame(.) 

ggplot(data=mono_pred_pie, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a"))

mono_pred_pie
sum(mono_pred_pie$Freq)


##-------------------activated/blocked ratio calc and graphing--------------------

#just manually putting the numbers into a df
pred_ratios <- data.frame(`Cell_Type` = c('Monocytes', 'DCs', 'pDCs', 'NK cells', 'CD8+ T cells', 'TAMs'), `Activating_Blocking_ratio` = c(1.1479, 1.2286, 0.8350, 0.8309, 1, 1)) %>%
  rownames_to_column()

#w/ log2 normalization
pred_ratios2 <- log2(pred_ratios$Activating_Blocking_ratio) %>%
  as.data.frame() %>%
  rownames_to_column()

log2_inc <- left_join(pred_ratios, pred_ratios2, by = "rowname")
log2_inc$. -> log2_inc$Log2norm
log2_inc <- log2_inc %>%
  select(., c("Log2norm", "Cell_Type", "Activating_Blocking_ratio"))

#plot it
ratio_plot <- ggplot(pred_ratios, aes(x=reorder(Cell_Type, Activating_Blocking_ratio), y=Activating_Blocking_ratio)) +
geom_bar(stat = "identity", fill = "black") +
  coord_flip() 
ratio_plot + theme(axis.line = element_line(linetype = "solid"), 
                   axis.ticks = element_line(size = 0.9), 
                   axis.title = element_text(size = 15), 
                   axis.text = element_text(size = 14, colour = "black"), 
                   axis.text.x = element_text(size = 14), 
                   panel.background = element_rect(fill = NA)) + labs( x="Cell Type", y="Activating/Blocking ratio") + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "Gray") + scale_y_continuous(expand = c(0.009, 0))


#plot it (log2norm)
ratio_plot2 <- ggplot(log2_inc, aes(x=reorder(Cell_Type, Log2norm), y=Log2norm)) +
  geom_bar(stat = "identity", fill = "black") +
  coord_flip() 
ratio_plot2 + theme(axis.line = element_line(linetype = "solid"), 
                   axis.ticks = element_line(size = 0.9), 
                   axis.title = element_text(size = 15), 
                   axis.text = element_text(size = 14, colour = "black"), 
                   axis.text.x = element_text(size = 14), 
                   panel.background = element_rect(fill = NA)) + labs( x="Cell Type", y="Log2(Activating/Blocking ratio)") +  
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "Gray")  

scale_y_continuous(expand = c(0.009, 0))



##-------------------Predictions seperated by responder status--------------------

#start with NK cells
NK_responder_pred_plot <- ggplot(Ad_SF7_fc_scored_NK_called, aes(x=SLAMF7_score, y=SLAMF7, color = outlier, group = outlier)) +
  geom_point(aes(size = outlier), alpha = 0.9) +
  scale_size_manual(values=c(1.5,1.5,0.9)) +       #change point sizes based on number of cells
  scale_color_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) +
  facet_wrap(~ response)
NK_responder_pred_plot + theme(axis.line = element_line(linetype = "solid"), 
                    axis.ticks = element_line(size = 0.9), 
                    axis.title = element_text(size = 15), 
                    axis.text = element_text(size = 14, colour = "black"), 
                    axis.text.x = element_text(size = 14), 
                    panel.background = element_rect(fill = NA))


NK_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_NK_called$outlier,
  condition = Ad_SF7_fc_scored_NK_called$condition, 
  response = Ad_SF7_fc_scored_NK_called$response) %>%
  as.data.frame(.) 

#first responders
NK_pred_pie_resp <- NK_pred_pie %>%
  filter(response == "yes")

ggplot(data=NK_pred_pie_resp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

NK_pred_pie_resp
sum(NK_pred_pie_resp$Freq) -> sum_nk_resp

#add column w/ percent
NK_pred_pie_resp <- NK_pred_pie_resp %>%
  mutate(., percent = (NK_pred_pie_resp$Freq/sum_nk_resp)*100)


#now non-responders
NK_pred_pie_nonresp <- NK_pred_pie %>%
  filter(response == "no")

ggplot(data=NK_pred_pie_nonresp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

NK_pred_pie_nonresp
sum(NK_pred_pie_nonresp$Freq) -> sum_nk_non

#add column w/ percent
NK_pred_pie_nonresp <- NK_pred_pie_nonresp %>%
  mutate(., percent = (NK_pred_pie_nonresp$Freq/sum_nk_non)*100)


##-------------make df w/ indiv. samples in above---------------#
#start w/ responders
NK_pred_ind <- dplyr::filter(Ad_SF7_fc_scored_NK_called, response == "yes") 
NK_pred_ind <- table(
  outlier = NK_pred_ind$outlier,
  sample_id = NK_pred_ind$sample_id) %>% 
  as.data.frame(.)
#remove non-responders and calculate % act. or blocked per sample
NK_pred_ind_resp <- NK_pred_ind[1:9,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#now non-responders
NK_pred_ind_non <- dplyr::filter(Ad_SF7_fc_scored_NK_called, response == "no") 
NK_pred_ind_non <- table(
  outlier = NK_pred_ind_non$outlier,
  sample_id = NK_pred_ind_non$sample_id) %>% 
  as.data.frame(.)
#remove responders and calculate % act. or blocked per sample
NK_pred_ind_non_resp <- NK_pred_ind_non[10:15,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#merge responder and non-responder dfs. 
NK_ind_preds <- bind_rows(NK_pred_ind_resp, NK_pred_ind_non_resp)
##--------------------------------------------------------------------#





#now DCs
dc_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_DC_called$outlier,
  condition = Ad_SF7_fc_scored_DC_called$condition, 
  response = Ad_SF7_fc_scored_DC_called$response) %>%
  as.data.frame(.) 

#first responders
dc_pred_pie_resp <- dc_pred_pie %>%
  filter(response == "yes")

ggplot(data=dc_pred_pie_resp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

dc_pred_pie_resp
sum(dc_pred_pie_resp$Freq) -> sum_dc_resp

#add column w/ percent
dc_pred_pie_resp <- dc_pred_pie_resp %>%
  mutate(., percent = (dc_pred_pie_resp$Freq/sum_dc_resp)*100)


#now non-responders
dc_pred_pie_nonresp <- dc_pred_pie %>%
  filter(response == "no")

ggplot(data=dc_pred_pie_nonresp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

dc_pred_pie_nonresp
sum(dc_pred_pie_nonresp$Freq) -> sum_dc_non

#add column w/ percent
dc_pred_pie_nonresp <- dc_pred_pie_nonresp %>%
  mutate(., percent = (dc_pred_pie_nonresp$Freq/sum_dc_non)*100)


##-------------make df w/ indiv. samples in above---------------#
#start w/ responders
dc_pred_ind <- dplyr::filter(Ad_SF7_fc_scored_DC_called, response == "yes") 
dc_pred_ind <- table(
  outlier = dc_pred_ind$outlier,
  sample_id = dc_pred_ind$sample_id) %>% 
  as.data.frame(.)
#remove non-responders and calculate % act. or blocked per sample
dc_pred_ind_resp <- dc_pred_ind[1:9,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#now non-responders
dc_pred_ind_non <- dplyr::filter(Ad_SF7_fc_scored_DC_called, response == "no") 
dc_pred_ind_non <- table(
  outlier = dc_pred_ind_non$outlier,
  sample_id = dc_pred_ind_non$sample_id) %>% 
  as.data.frame(.)
#remove responders and calculate % act. or blocked per sample
dc_pred_ind_non_resp <- dc_pred_ind_non[10:15,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#merge responder and non-responder dfs. 
dc_ind_preds <- bind_rows(dc_pred_ind_resp, dc_pred_ind_non_resp)
##--------------------------------------------------------------------#





#now TAMs
tam_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_tam_called$outlier,
  condition = Ad_SF7_fc_scored_tam_called$condition, 
  response = Ad_SF7_fc_scored_tam_called$response) %>%
  as.data.frame(.) 

#first responders
tam_pred_pie_resp <- tam_pred_pie %>%
  filter(response == "yes")

ggplot(data=tam_pred_pie_resp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

tam_pred_pie_resp
sum(tam_pred_pie_resp$Freq) -> sum_tam_resp

#add column w/ percent
tam_pred_pie_resp <- tam_pred_pie_resp %>%
  mutate(., percent = (tam_pred_pie_resp$Freq/sum_tam_resp)*100)


#now non-responders
tam_pred_pie_nonresp <- tam_pred_pie %>%
  filter(response == "no")

ggplot(data=tam_pred_pie_nonresp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

tam_pred_pie_nonresp
sum(tam_pred_pie_nonresp$Freq) -> sum_tam_non

#add column w/ percent
tam_pred_pie_nonresp <- tam_pred_pie_nonresp %>%
  mutate(., percent = (tam_pred_pie_nonresp$Freq/sum_tam_non)*100)


##-------------make df w/ indiv. samples in above---------------#
#start w/ responders
tam_pred_ind <- dplyr::filter(Ad_SF7_fc_scored_tam_called, response == "yes") 
tam_pred_ind <- table(
  outlier = tam_pred_ind$outlier,
  sample_id = tam_pred_ind$sample_id) %>% 
  as.data.frame(.)
#remove non-responders and calculate % act. or blocked per sample
tam_pred_ind_resp <- tam_pred_ind[1:9,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#now non-responders
tam_pred_ind_non <- dplyr::filter(Ad_SF7_fc_scored_tam_called, response == "no") 
tam_pred_ind_non <- table(
  outlier = tam_pred_ind_non$outlier,
  sample_id = tam_pred_ind_non$sample_id) %>% 
  as.data.frame(.)
#remove responders and calculate % act. or blocked per sample
tam_pred_ind_non_resp <- tam_pred_ind_non[10:15,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#merge responder and non-responder dfs. 
tam_ind_preds <- bind_rows(tam_pred_ind_resp, tam_pred_ind_non_resp)
##--------------------------------------------------------------------#





#now CD8 T cells
cd8_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_cd8_called$outlier,
  condition = Ad_SF7_fc_scored_cd8_called$condition, 
  response = Ad_SF7_fc_scored_cd8_called$response) %>%
  as.data.frame(.) 

#first responders
cd8_pred_pie_resp <- cd8_pred_pie %>%
  filter(response == "yes")

ggplot(data=cd8_pred_pie_resp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

cd8_pred_pie_resp
sum(cd8_pred_pie_resp$Freq) -> sum_cd8_resp

#add column w/ percent
cd8_pred_pie_resp <- cd8_pred_pie_resp %>%
  mutate(., percent = (cd8_pred_pie_resp$Freq/sum_cd8_resp)*100)


#now non-responders
cd8_pred_pie_nonresp <- cd8_pred_pie %>%
  filter(response == "no")

ggplot(data=cd8_pred_pie_nonresp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

cd8_pred_pie_nonresp
sum(cd8_pred_pie_nonresp$Freq) -> sum_cd8_non

#add column w/ percent
cd8_pred_pie_nonresp <- cd8_pred_pie_nonresp %>%
  mutate(., percent = (cd8_pred_pie_nonresp$Freq/sum_cd8_non)*100)



##-------------make df w/ indiv. samples in above---------------#
#start w/ responders
cd8_pred_ind <- dplyr::filter(Ad_SF7_fc_scored_cd8_called, response == "yes") 
cd8_pred_ind <- table(
  outlier = cd8_pred_ind$outlier,
  sample_id = cd8_pred_ind$sample_id) %>% 
  as.data.frame(.)
#remove non-responders and calculate % act. or blocked per sample
cd8_pred_ind_resp <- cd8_pred_ind[1:9,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#now non-responders
cd8_pred_ind_non <- dplyr::filter(Ad_SF7_fc_scored_cd8_called, response == "no") 
cd8_pred_ind_non <- table(
  outlier = cd8_pred_ind_non$outlier,
  sample_id = cd8_pred_ind_non$sample_id) %>% 
  as.data.frame(.)
#remove responders and calculate % act. or blocked per sample
cd8_pred_ind_non_resp <- cd8_pred_ind_non[10:15,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#merge responder and non-responder dfs. 
cd8_ind_preds <- bind_rows(cd8_pred_ind_resp, cd8_pred_ind_non_resp)
##--------------------------------------------------------------------#





#now pDCs
pdc_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_pdc_called$outlier,
  condition = Ad_SF7_fc_scored_pdc_called$condition, 
  response = Ad_SF7_fc_scored_pdc_called$response) %>%
  as.data.frame(.) 




##-------------make df w/ indiv. samples in above---------------#
#start w/ responders
pdc_pred_ind <- dplyr::filter(Ad_SF7_fc_scored_pdc_called, response == "yes") 
pdc_pred_ind <- table(
  outlier = pdc_pred_ind$outlier,
  sample_id = pdc_pred_ind$sample_id) %>% 
  as.data.frame(.)
#remove non-responders and calculate % act. or blocked per sample
pdc_pred_ind_resp <- pdc_pred_ind[1:9,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#now non-responders
pdc_pred_ind_non <- dplyr::filter(Ad_SF7_fc_scored_pdc_called, response == "no") 
pdc_pred_ind_non <- table(
  outlier = pdc_pred_ind_non$outlier,
  sample_id = pdc_pred_ind_non$sample_id) %>% 
  as.data.frame(.)
#remove responders and calculate % act. or blocked per sample
pdc_pred_ind_non_resp <- pdc_pred_ind_non[10:15,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#merge responder and non-responder dfs. 
pdc_ind_preds <- bind_rows(pdc_pred_ind_resp, pdc_pred_ind_non_resp)
##--------------------------------------------------------------------#
 

#first responders
pdc_pred_pie_resp <- pdc_pred_pie %>%
  filter(response == "yes")

ggplot(data=pdc_pred_pie_resp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

pdc_pred_pie_resp
sum(pdc_pred_pie_resp$Freq) -> sum_pdc_resp

#add column w/ percent
pdc_pred_pie_resp <- pdc_pred_pie_resp %>%
  mutate(., percent = (pdc_pred_pie_resp$Freq/sum_pdc_resp)*100)


#now non-responders
pdc_pred_pie_nonresp <- pdc_pred_pie %>%
  filter(response == "no")

ggplot(data=pdc_pred_pie_nonresp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

pdc_pred_pie_nonresp
sum(pdc_pred_pie_nonresp$Freq)
sum(pdc_pred_pie_nonresp$Freq) -> sum_pdc_non

#add column w/ percent
pdc_pred_pie_nonresp <- pdc_pred_pie_nonresp %>%
  mutate(., percent = (pdc_pred_pie_nonresp$Freq/sum_pdc_non)*100)



#now monocytes
mono_pred_pie <- table(
  outlier = Ad_SF7_fc_scored_mono_called$outlier,
  condition = Ad_SF7_fc_scored_mono_called$condition, 
  response = Ad_SF7_fc_scored_mono_called$response) %>%
  as.data.frame(.) 

#first responders
mono_pred_pie_resp <- mono_pred_pie %>%
  filter(response == "yes")

ggplot(data=mono_pred_pie_resp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

mono_pred_pie_resp
sum(mono_pred_pie_resp$Freq) -> sum_mono_resp

#add column w/ percent
mono_pred_pie_resp <- mono_pred_pie_resp %>%
  mutate(., percent = (mono_pred_pie_resp$Freq/sum_mono_resp)*100)


#now non-responders
mono_pred_pie_nonresp <- mono_pred_pie %>%
  filter(response == "no")

ggplot(data=mono_pred_pie_nonresp, aes(x="", y=Freq, fill=outlier)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#ff9900", "#9933ff", "#1a1a1a")) 

mono_pred_pie_nonresp
sum(mono_pred_pie_nonresp$Freq) -> sum_mono_non

#add column w/ percent
mono_pred_pie_nonresp <- mono_pred_pie_nonresp %>%
  mutate(., percent = (mono_pred_pie_nonresp$Freq/sum_mono_non)*100)


##-------------make df w/ indiv. samples in above---------------#
#start w/ responders
mono_pred_ind <- dplyr::filter(Ad_SF7_fc_scored_mono_called, response == "yes") 
mono_pred_ind <- table(
  outlier = mono_pred_ind$outlier,
  sample_id = mono_pred_ind$sample_id) %>% 
  as.data.frame(.)
#remove non-responders and calculate % act. or blocked per sample
mono_pred_ind_resp <- mono_pred_ind[1:9,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#now non-responders
mono_pred_ind_non <- dplyr::filter(Ad_SF7_fc_scored_mono_called, response == "no") 
mono_pred_ind_non <- table(
  outlier = mono_pred_ind_non$outlier,
  sample_id = mono_pred_ind_non$sample_id) %>% 
  as.data.frame(.)
#remove responders and calculate % act. or blocked per sample
mono_pred_ind_non_resp <- mono_pred_ind_non[10:15,] %>%
  group_by(sample_id) %>%
  dplyr::mutate(., total = sum(Freq)) %>%
  group_by(sample_id, outlier) %>%
  dplyr::mutate(., percent = Freq/total) %>%
  dplyr::filter(outlier != "neither") %>%
  select(outlier, sample_id, percent)

#merge responder and non-responder dfs. 
mono_ind_preds <- bind_rows(mono_pred_ind_resp, mono_pred_ind_non_resp)
##--------------------------------------------------------------------#


##--------make aggregate df of predictions for all cell types per sample
all_pred_ind <- bind_rows(mono_ind_preds, pdc_ind_preds, dc_ind_preds, tam_ind_preds, cd8_ind_preds, NK_ind_preds, .id = "Cell Type") %>%
  as.data.frame()

#export
write.xlsx(all_pred_ind, "~/B16_TME_profiling_batch2_analys/batch2_preds_iEN.xlsx")




##----------now make plot of aggregate changes b/w responders and non-responders

#going with bar plot, set to 0 in the middle and with each bar being the percent change in an 
#outlier group (activated or inhibited) b/w responders and non-responders

#merge all DFs containing precent data and calculate percent change b/w response

#mono
merge1 <- bind_rows(mono_pred_pie_nonresp, mono_pred_pie_resp, .id = NULL) %>%
  add_column(., cell_type = "Monocytes") %>%
  filter(., outlier != "neither") 

#precent change calc 
mono_chg <- merge1 %>%
  mutate(., act_perc_chg = ((merge1[3,5]-merge1[1,5])/merge1[1,5])*100) %>%
  mutate(., block_prec_chg = ((merge1[4,5]-merge1[2,5])/merge1[2,5])*100)



#pdc
merg2 <- bind_rows(pdc_pred_pie_nonresp, pdc_pred_pie_resp, .id = NULL) %>%
  add_column(., cell_type = "pDCs") %>%
  filter(., outlier != "neither") 

#precent change calc 
pdc_chg <- merg2 %>%
  mutate(., act_perc_chg = ((merg2[3,5]-merg2[1,5])/merg2[1,5])*100) %>%
  mutate(., block_prec_chg = ((merg2[4,5]-merg2[2,5])/merg2[2,5])*100)



#cd8
merg3 <- bind_rows(cd8_pred_pie_nonresp, cd8_pred_pie_resp, .id = NULL) %>%
  add_column(., cell_type = "CD8+ T cells") %>%
  filter(., outlier != "neither") 

#precent change calc 
cd8_chg <- merg3 %>%
  mutate(., act_perc_chg = ((merg3[3,5]-merg3[1,5])/merg3[1,5])*100) %>%
  mutate(., block_prec_chg = ((merg3[4,5]-merg3[2,5])/merg3[2,5])*100)




#TAM
merg4 <- bind_rows(tam_pred_pie_nonresp, tam_pred_pie_resp, .id = NULL) %>%
  add_column(., cell_type = "TAMs") %>%
  filter(., outlier != "neither") 

#precent change calc 
tam_chg <- merg4 %>%
  mutate(., act_perc_chg = ((merg4[3,5]-merg4[1,5])/merg4[1,5])*100) %>%
  mutate(., block_prec_chg = ((merg4[4,5]-merg4[2,5])/merg4[2,5])*100)


#DC
merg5 <- bind_rows(dc_pred_pie_nonresp, dc_pred_pie_resp, .id = NULL) %>%
  add_column(., cell_type = "DCs") %>%
  filter(., outlier != "neither") 

#precent change calc
dc_chg <- merg5 %>%
  mutate(., act_perc_chg = ((merg5[3,5]-merg5[1,5])/merg5[1,5])*100) %>%
  mutate(., block_prec_chg = ((merg5[4,5]-merg5[2,5])/merg5[2,5])*100)


#NK
merg6 <- bind_rows(NK_pred_pie_nonresp, NK_pred_pie_resp, .id = NULL) %>%
  add_column(., cell_type = "NK cells") %>%
  filter(., outlier != "neither") 

#precent change calc 
nk_chg <- merg6 %>%
  mutate(., act_perc_chg = ((merg6[3,5]-merg6[1,5])/merg6[1,5])*100) %>%
  mutate(., block_prec_chg = ((merg6[4,5]-merg6[2,5])/merg6[2,5])*100)



#now just manually putting the numbers into a df
chg_ratios <- data.frame(`Cell_type_response` = c('Monocytes_activate', 'Monocytes_block', 'DCs_activate', 'DCs_block', 'pDCs_activate', 'pDCs_block', 'NK cells_activate', 'NK cells_block', 'CD8+ T cells_activate', 'CD8+ T cells_block', 'TAMs_activate', 'TAMs_block'), 
                          `Percent_change` = c(81.44, 38.42, 26.31, -21.72, 234.11, -41.47, 16.31, 85.93, 192.35, 0.37, 28.08, -18.79))


#plot it
chg_plot <- ggplot(chg_ratios, aes(x=reorder(Cell_type_response, Percent_change), y=Percent_change)) +
  geom_bar(stat = "identity", position=position_dodge(), fill = "black") +
  coord_flip() 
chg_plot + theme(axis.line = element_line(linetype = "solid"), 
                   axis.ticks = element_line(size = 0.9), 
                   axis.title = element_text(size = 15), 
                   axis.text = element_text(size = 14, colour = "black"), 
                   axis.text.x = element_text(size = 14), 
                   panel.background = element_rect(fill = NA)) + labs( x="Cell Type", y="Percent change from non-responder")  




##----------Percent change in ratios---------

#now I calcualte the activating/blocking ratio for each cell type b/w responders and non-resp
#and calculate the percent change b/w responses and plot that. 

#now just manually putting the numbers into a df
ratios_by_response <- data.frame(`Cell_type_response` = c('Monocytes_responder', 'Monocytes_non-responder', 'DCs_responder', 'DCs_non-responder', 'pDCs_responder', 'pDCs_non-responder', 'NK cells_responder', 'NK cells_non-responder', 'CD8+ T cells_responder', 'CD8+ T cells_non-responder', 'TAMs_responder', 'TAMs_non-responder'), 
                         `Activating_Blocking_ratio` = c(1.177663, 0.8984379, 1.352941, 0.8385417, 1.328671, 0.2327586, 0.7699119, 1.230769, 1.098511, 0.3771289, 1.095238, 0.6944445))

#calc percent change in ratios

#mono
((ratios_by_response[1,2]-ratios_by_response[2,2])/ratios_by_response[2,2])*100
#mono_ratio = 31.07895  #this means a 31% increase in activation/blocking ratio in responders 
#compared to non-responders

#DCs
((ratios_by_response[3,2]-ratios_by_response[4,2])/ratios_by_response[4,2])*100
# DC = 61.34451

#pDCs
((ratios_by_response[5,2]-ratios_by_response[6,2])/ratios_by_response[6,2])*100
# pDC = 470.8365

#NK
((ratios_by_response[7,2]-ratios_by_response[8,2])/ratios_by_response[8,2])*100
# NK = -37.44465

#cd8
((ratios_by_response[9,2]-ratios_by_response[10,2])/ratios_by_response[10,2])*100
# cd8 = 191.2826

#TAM
((ratios_by_response[11,2]-ratios_by_response[12,2])/ratios_by_response[12,2])*100
# tam = 57.71426

#manually put all the values in a df
ratio_resp_chg <- data.frame(`Cell_Type` = c('Monocytes', 'DCs', 'pDCs', 'NK cells', 'CD8+ T cells', 'TAMs'), 
                                 `Percent_change_in_Activating_Blocking_ratio` = c(31.07895, 61.34451, 470.8365, -37.44465, 191.2826, 57.71426))

#plot it
ratio_resp_chg_plot <- ggplot(ratio_resp_chg, aes(x=reorder(Cell_Type, Percent_change_in_Activating_Blocking_ratio), y=Percent_change_in_Activating_Blocking_ratio)) +
  geom_bar(stat = "identity", position=position_dodge(), fill = "black") +
  coord_flip() 
ratio_resp_chg_plot + theme(axis.line = element_line(linetype = "solid"), 
                 axis.ticks = element_line(size = 0.9), 
                 axis.title = element_text(size = 15), 
                 axis.text = element_text(size = 14, colour = "black"), 
                 axis.text.x = element_text(size = 14), 
                 panel.background = element_rect(fill = NA)) + labs( x="Cell Type", y="Percent change in Activating/Blocking ratio from non-responder") +
  scale_y_continuous(expand = c(0.009, 0))


#plot ratios_by_response data

#w/ log2 normalization
ratios_by_response <- rownames_to_column(ratios_by_response)
ratios_by_response2 <- log2(ratios_by_response$Activating_Blocking_ratio) %>%
  as.data.frame() %>%
  rownames_to_column()

log2_inc2 <- left_join(ratios_by_response, ratios_by_response2, by = "rowname")
log2_inc2$. -> log2_inc2$Log2norm
log2_inc2 <- log2_inc2 %>%
  select(., c("Log2norm", "Cell_type_response", "Activating_Blocking_ratio"))


ratio_resp_plot <- ggplot(log2_inc2, aes(x=reorder(Cell_type_response, Log2norm), y=Log2norm)) +
  geom_bar(stat = "identity", position=position_dodge(), fill = "black") +
  coord_flip() 
ratio_resp_plot + theme(axis.line = element_line(linetype = "solid"), 
                            axis.ticks = element_line(size = 0.9), 
                            axis.title = element_text(size = 15), 
                            axis.text = element_text(size = 14, colour = "black"), 
                            axis.text.x = element_text(size = 14), 
                            panel.background = element_rect(fill = NA)) + labs( x="Cell Type Response", y="Log2(Activating/Blocking ratio)") 

  scale_y_continuous(expand = c(0.009, 0))





