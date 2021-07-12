#by:Patrick O'Connell

### Attempt at predicting parameters which effictivly discriminate responders from non-responders
#in mice who recieved Ad-SF7-Fc injection into B16-F10 tumors. 

#3/27/2021

#I am taking a machine learning-based approach here and using the iEN algorithm. This was developed
#specificially for high dimensional immune data and allows the user to weight based on known immunological facts.

#https://github.com/Teculos/immunological-EN


setwd("~/iEN_modeling")

#Packages
library(iEN)
library(tidyverse)
library(readxl)
library(ggplot2)
library(xlsx)
library(parallel)
library(caret)


###------------------------Analysis on Ad-SF7-Fc B16 data----------------------------###

#First I will collect and pool my data together in one giant df



##List of data included:
# 1) T cell phenotype data
# 2) Frequency of all immune cell types in the TME (only used cell types common to both analyses)
# 3) Median expression of all markers on all cell types in TME
# 4) Tumor lysate cytokine data
# 5) Plasma cytokine data
# 6) SLAMF7 activating/blocking predictions


#load in all data, aggregate, and identify where missing values are

#missing: 
# 1936 and 1937 lysate cytokine data (both batch1 respodners)

setwd("~/iEN_modeling/data_dfs")
#load the data, and remove unnecessary columns and name columns appropiatly (must be unique)
batch2_pred_iEN <- read_excel("batch2_preds_iEN.xlsx") %>%
  pivot_wider(names_from = c("cell_type", "prediction"), values_from = freq)
batch2_pred_iEN$sample_id <- gsub(batch2_pred_iEN$sample_id, pattern = "SLAMF7", replacement = "SF7") 

batch1_pred_iEN <- read_excel("batch1_preds_iEN.xlsx") %>%
  pivot_wider(names_from = c("cell_type", "prediction"), values_from = freq)
batch1_pred_iEN$sample_id <- gsub(batch1_pred_iEN$sample_id, pattern = "SLAMF7", replacement = "SF7") 

#merg_pred data
merg_pred_data <- bind_rows(batch1_pred_iEN, batch2_pred_iEN)

batch2_clust_expr_iEN <- read_excel("batch2_clust_expr_iEN.xlsx") %>%
  select(!one_of("responder")) %>%
  pivot_wider(names_from = merg_clust, values_from = c("CD45",       "CD11b",      "CD40",       "CD4",       
                                                       "CD8",        "NK1.1",      "CD19",       "SLAMF7",     "F4/80",     "Ly6C",       "Ly6G",      
                                                       "MHC-II",     "PD-L1",      "CD11c",     "B220",      "PD-1",       "IgD",        "CCR2",      
                                                       "CD38",       "CD206",      "AF",         "CD3",       "Fascin",     "Lyve-1",     "CD90",      
                                                       "CCR7"))

batch1_clust_expr_iEN <- read_excel("batch1_clust_expr_iEN.xlsx") %>%
  select(!one_of("responder")) %>%
  pivot_wider(names_from = merg_clust, values_from = c("CD45",       "CD11b",      "CD40",       "CD4",       
                                                       "CD8",        "NK1.1",      "CD19",       "SLAMF7",     "F4/80",     "Ly6C",       "Ly6G",      
                                                       "MHC-II",     "PD-L1",      "CD11c",     "B220",      "PD-1",       "IgD",        "CCR2",      
                                                       "CD38",       "CD206",      "AF",         "CD3",       "Fascin",     "Lyve-1",     "CD90",      
                                                       "CCR7"))
batch1_clust_expr_iEN$sample_id <- gsub(batch1_clust_expr_iEN$sample_id, pattern = "SLAMF7", replacement = "SF7")  

#merg expr data
merg_expr_data <- bind_rows(batch1_clust_expr_iEN, batch2_clust_expr_iEN)


batch2_clust_freq_iEN <- read_excel("batch2_clust_freq_iEN.xlsx")  %>%
  select(!one_of("responder")) %>%
  spread(., merg_clust, freq)
batch2_clust_freq_iEN$sample_id <- gsub(batch2_clust_freq_iEN$sample_id, pattern = "SLAMF7", replacement = "SF7")  

batch1_clust_freq_iEN <- read_excel("batch1_clust_freq_iEN.xlsx") %>%
  select(!one_of("responder")) %>%
  spread(., merg_clust, freq)
batch1_clust_freq_iEN$sample_id <- gsub(batch1_clust_freq_iEN$sample_id, pattern = "SLAMF7", replacement = "SF7")  

#merg freq data
merg_freq_data <- bind_rows(batch1_clust_freq_iEN, batch2_clust_freq_iEN)

#batch1_2_lysate_cyto_df_iEN <- read_excel("batch1_2_lysate_cyto_df_iEN.xlsx")

batch1_2_plasma_cyto_df_iEN <- read_excel("batch1_2_plasma_cyto_df_iEN.xlsx")  %>%
  filter(., time_point == "Day_0")   
colnames(batch1_2_plasma_cyto_df_iEN) <- paste(colnames(batch1_2_plasma_cyto_df_iEN), "plasma_day0", sep = "_")
batch1_2_plasma_cyto_df_iEN <- rename(batch1_2_plasma_cyto_df_iEN, sample_id = sample_id_plasma_day0)
batch1_2_plasma_cyto_df_iEN <- select(batch1_2_plasma_cyto_df_iEN, !one_of("location_plasma_day0", "batch_plasma_day0", "responder_plasma_day0", "location_plasma_day0", "time_point_plasma_day0"))

batch1_2_plasma_cyto_df_iEN2 <- read_excel("batch1_2_plasma_cyto_df_iEN.xlsx")  %>%
  filter(., time_point == "Day_1")   
colnames(batch1_2_plasma_cyto_df_iEN2) <- paste(colnames(batch1_2_plasma_cyto_df_iEN2), "plasma_day1", sep = "_")
batch1_2_plasma_cyto_df_iEN2 <- rename(batch1_2_plasma_cyto_df_iEN2, sample_id = sample_id_plasma_day1)
batch1_2_plasma_cyto_df_iEN2 <- select(batch1_2_plasma_cyto_df_iEN2, !one_of("location_plasma_day1", "batch_plasma_day1", "responder_plasma_day1", "location_plasma_day1", "time_point_plasma_day1"))
#merge them
batch1_2_plasma_cyto_df_iEN_clean <- left_join(batch1_2_plasma_cyto_df_iEN2, batch1_2_plasma_cyto_df_iEN, by = "sample_id")


batch2_T_cell_data_iEN <- read_excel("batch2_T_cell_data_iEN.xlsx") %>%
  select(!one_of("responder")) %>%
  spread(., cell_type, freq)

batch1_T_cell_data_iEN <- read_excel("batch1_T_cell_data_iEN.xlsx") %>%
  select(!one_of("responder")) %>%
  spread(., cell_type, freq)

#merg_T cell data
merg_T_cell_data <- bind_rows(batch1_T_cell_data_iEN, batch2_T_cell_data_iEN)


batch1_2_lysate_cyto_df_iEN_preds <- read_excel("batch1_2_lysate_cyto_df_iENP_w_preds.xlsx") %>%
  select(., !one_of("location", "batch", "responder"))   #contains pred values
colnames(batch1_2_lysate_cyto_df_iEN_preds) <- paste(colnames(batch1_2_lysate_cyto_df_iEN_preds), "lysate", sep = "_")
batch1_2_lysate_cyto_df_iEN_preds <- rename(batch1_2_lysate_cyto_df_iEN_preds, sample_id = sample_id_lysate)
#change name of one mislabeled sample
batch1_2_lysate_cyto_df_iEN_preds$sample_id <- gsub(batch1_2_lysate_cyto_df_iEN_preds$sample_id, pattern = "1798", replacement = "0798")  


#for missing tumor lysate samples, I'll take the median expression for each cytokine from responding mice from both batches (since they were responders)
med_pred <- batch1_2_lysate_cyto_df_iEN %>%
  filter(responder == "yes") 
med_pred <- med_pred[,5:27] 
preds <- med_pred %>%
  as.data.frame(.) %>%
  summarise_all(median)
#since I'm lazy I just export this and manually copy it into the oroginal df and reload
write.xlsx(preds, "~/iEN_modeling/data_dfs/lysate_preds.xlsx")


#Now aggregate all data above into single df of proper form

#Check that sample_id names for all dfs match
unique(batch1_clust_expr_iEN$sample_id, batch1_2_plasma_cyto_df_iEN_clean$sample_id)
setequal(batch2_clust_expr_iEN$sample_id, batch2_T_cell_data_iEN$sample_id)

agg_data <- full_join(batch1_2_lysate_cyto_df_iEN_preds, batch1_2_plasma_cyto_df_iEN_clean, by = "sample_id")
agg_data <- full_join(agg_data, merg_T_cell_data, by = "sample_id")
agg_data <- full_join(agg_data, merg_freq_data, by = "sample_id")
agg_data <- full_join(agg_data, merg_expr_data, by = "sample_id")
agg_data <- full_join(agg_data, merg_pred_data, by = "sample_id")

agg_data <- column_to_rownames(agg_data, var = "sample_id")
agg_dat_mtx <- as.matrix(agg_data) #input this into iEN
saveRDS(agg_dat_mtx, "~/iEN_modeling/agg_dat_mtx")
agg_dat_mtx <- readRDS("~/iEN_modeling/agg_dat_mtx")

##-------------generate priors-------------------##

#this is just a named number vector (as.numeric) w/ the rownames matching the colnames of the data_matrix
#and the values for each one the weighting. 

#My situation is a bit different from the example paper, since I have differnt types of readouts
#I will score the expression variables from TME panel and T cell panel just as they did in the paper. 
#I will not score TME panel frequency data since I do not want to bias it. (all have value of 1)
#Day_0 plasma cytokines do not get scored (all have value of 1)
#Day_1 plasma cytokines get scored based on expected changes from day_0 (from an Adenovirus perspective)
#lysate cytokine data does not get scored since I do not want to bias. (all have value of 1)
#SF7 blocked/activated prediction data does not get scored since we have no ground truth data on SF7-Fc ability to block or activate specific cell types (all values set to 1)

prior_temp <- t(agg_dat_mtx)
write.xlsx(prior_temp, "~/iEN_modeling/prior_temp.xlsx")

#read priors back in
SF7_Fc_priors<- read_excel("prior_temp.xlsx") %>%
  column_to_rownames(., var = "...1") 
#get data in right format
numbs <- SF7_Fc_priors[[1]]
rownames(SF7_Fc_priors) -> SF7_names 
names(numbs) <- SF7_names
numbs -> SF7_priors_clean
rm(numbs)

##----------Set up variables for iEN------------##

#Set response variable
resp_var <- c("yes", "no", "yes", "no", "yes", "no", "yes", "yes", "no", "no", "yes", "yes", "yes") %>%
  as.factor()

alphaGrid <- seq(0,1, length.out=11) 
phiGrid <- c(0, exp(seq(log(1),log(100), length.out=10))) 
nlambda <- 10 
lambdas=NULL
ncores <- detectCores()-1 
eval <- "ROCAUC"
family <- "binomial"
intercept <- FALSE
standardize <- TRUE
center <- TRUE

K <- 10
repeats <- 100 
set.seed(1994)
for(i in seq(repeats)){
  foldid <- caret::createFolds(seq(nrow(agg_dat_mtx)),k=K,list=FALSE)
}

#Run the algorithm
ien.SF7_Fc_100=cv_iEN(agg_dat_mtx, resp_var, foldid, alphaGrid, phiGrid, nlambda, lambdas,  SF7_priors_clean, ncores, eval, family, intercept, standardize, center)
#save model
saveRDS(ien.SF7_Fc_100, "~/iEN_modeling/ien.SF7_Fc_100")
#load the model back
ien.SF7_Fc_100 <- readRDS("~/iEN_modeling/ien.SF7_Fc_100")

#Run the algorithm w/o priors
ien.SF7_Fc_100_wo=cv_iEN(agg_dat_mtx, resp_var, foldid, alphaGrid, phiGrid=c(0), nlambda, lambdas,  SF7_priors_clean, ncores, eval, family, intercept, standardize, center)
#save the model
saveRDS(ien.SF7_Fc_100_wo, "~/iEN_modeling/ien.SF7_Fc_100_wo")
#load the model back
ien.SF7_Fc_100_wo <- readRDS("~/iEN_modeling/ien.SF7_Fc_100_wo")

#both run fine and have solid prediciton strength. Adding in priors increases the CV from 0.8 to 0.9 showing that my priors are useful. 

print(ien.SF7_Fc_100)
print(ien.SF7_Fc_100_wo)


#make ROC curves

#my data and plot
#predicitons are from the iEN object and re-code yes no classifier to 1,0
SF7_ROC_df <- data.frame(predictions = (ien.SF7_Fc_100@cv.preds), 
                         labels = c(1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1))

SF7_ROC_df_wo <- data.frame(predictions = (ien.SF7_Fc_100_wo@cv.preds), 
                            labels = c(1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1))

library(pROC)

#now my data
pROC_obj_SF7 <- roc(SF7_ROC_df$labels,SF7_ROC_df$predictions,
                smoothed = T,
                ci=F, ci.alpha=0.9, stratified=FALSE,
                plot=TRUE, auc.polygon=F, max.auc.polygon=F, grid=F,
                print.auc=TRUE, show.thres=TRUE)

pROC_obj_SF7_wo <- roc(SF7_ROC_df_wo$labels,SF7_ROC_df_wo$predictions,
                    smoothed = T,
                    ci=F, ci.alpha=0.9, stratified=FALSE,
                    plot=TRUE, auc.polygon=F, max.auc.polygon=F, grid=F,
                    print.auc=TRUE, show.thres=TRUE)

#export sensitivity and specifity data from pROC_ob_SF7 and make a ggplot from it. 
#and repeat for non-prior analysis and plot them too on top
both <- list(pROC_obj_SF7, pROC_obj_SF7_wo)

plot1 <- ggroc(list(`With Priors (AUC=0.925)` = pROC_obj_SF7,`Without Priors (AUC=0.825)` = pROC_obj_SF7_wo), size = 1.0, aes= c("linetype", "color")) 
plot2 <- plot1 + scale_color_manual(values = c("black", "#808080")) 
plot2 + theme(axis.line = element_line(linetype = "solid", size = 1.0), 
    axis.ticks = element_line(colour = "black", 
        size = 0.9), panel.grid.major = element_line(colour = "khaki4", 
        linetype = "blank"), axis.title = element_text(size = 19), 
    axis.text = element_text(size = 16), 
    legend.text = element_text(size = 16), 
    legend.title = element_text(colour = NA), 
    panel.background = element_rect(fill = NA), 
    plot.background = element_rect(colour = "white", 
        size = 1.1), legend.key = element_rect(fill = NA), 
    legend.background = element_rect(fill = NA))
##--------------------------------------------------------------##

#run wilcoxin signed rank test on predictive values from iEN comparing responders to non-responders
#this gives P value statistic for how good the model can seperate the groups. 
SF7_ROC_df$labels <- as.factor(SF7_ROC_df$labels)
levels(SF7_ROC_df$labels) <- c("Non-Responder", "Responder")

#calc p value
wilcox.test(SF7_ROC_df$predictions~SF7_ROC_df$labels)$p.value

#plot it!
p_val_plot <- ggplot(SF7_ROC_df, aes(x=labels, y=predictions)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.alpha = 0) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(1), dotsize = 0.5) 
p_val_plot + theme(axis.line = element_line(linetype = "solid"), 
    axis.ticks = element_line(colour = "black", 
        size = 0.9), axis.title = element_text(size = 19), 
    axis.text = element_text(size = 16), 
    panel.background = element_rect(fill = NA)) +
  labs(y = "iEN predicted values", x = "")
##--------------------------------------------------------------##

