#Written by: Tony Culos
#Edited by: Patrick O'Connell
#4/27/2021

#Network Generation  Script 
#using this to genreate correlation network figure and pull out most significant parameters seperating responders 
#fron non-responders

setwd("~/iEN_modeling")
#packages
library(iEN)
library(tidyverse)
library(readxl)
library(ggplot2)
library(xlsx)
library(pROC)

#--------------load my data-----------

ien.SF7_Fc_100 <- readRDS("~/iEN_modeling/ien.SF7_Fc_100")
agg_dat_mtx <- readRDS("~/iEN_modeling/agg_dat_mtx")

#-------------------start analysis------------
Data_x = agg_dat_mtx  
library(Hmisc)
library(igraph)
corg<-graph.adjacency(1-(as.matrix((cor(Data_x, method="spearman")))),weighted=TRUE)
mst<-minimum.spanning.tree(corg)
E(mst)$weight=E(mst)$weight*5
mst=as.undirected(mst)

r<-array(dim=c(ncol(Data_x),ncol(Data_x)))
P<-array(dim=c(ncol(Data_x),ncol(Data_x)))

temp=rcorr(Data_x, type="spearman")
for(i in seq(ncol(Data_x))){
  for(j in seq(ncol(Data_x))){
    P[i,j] = cor.test(Data_x[,i],Data_x[,j], method = "spearman")$p.value
    r[i,j] = cor.test(Data_x[,i],Data_x[,j], method = "spearman")$estimate
    
    if(i == j){
      r[i,j] = 1
      P[i,j] = NA      
    }
  }
}

rownames(P) = colnames(Data_x)
colnames(P) = colnames(Data_x)

rownames(r) = colnames(Data_x)
colnames(r) = colnames(Data_x)

temp=list(r=r,P=P)

pvs=temp$P
pvs=pvs+min(pvs[pvs!=0],na.rm=TRUE)
cors=temp$r
pvs[log10((pvs*(ncol(pvs)*(ncol(pvs)-1)/2)))> log10(0.05)]=NA #this is bonferroni correction 
pvs=-log10(pvs)
pvs[pvs==0]=NA
pvs[cors==1]=NA
pvs[is.na(pvs)]=0
g=graph_from_adjacency_matrix(pvs, mode='undirected', weighted=TRUE)


myweights=vector()
myedges=get.edgelist(g)
for (i in seq(nrow(myedges))){
  myweights[i]=1/abs(cors[myedges[i,1], myedges[i,2]])
}
E(g)$weight=myweights


set.seed(13245)
library(Rtsne)
z=Rtsne(pvs, check_duplicates = FALSE, 2)
lo=z$Y

edge.width=0.5


# calculating the AUROC for each feature as a predictor of response. 

##-----calc AUROC for each feature----------
response_codes <- data.frame(labels = c(1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1))
agg_dat_mtx2 <- as.data.frame(agg_dat_mtx)
comb.input <- bind_cols(agg_dat_mtx2, response_codes)

ROCvector<- colnames(agg_dat_mtx)
mypvs <- rep(NA,363)
for (i in seq_along(ROCvector)) {
  a<-ROCvector[i]
  mypvs[i] <- roc_(data=comb.input, "labels", as.character(a))$auc
}

#sizes based off of AUROC for each measurement
sizes=mypvs*4.5 #to globally size nodes as appropriate. 
V(g)$label.cex=0.05 #label size
cols = "#e6e6e6" #node and label color

#visualize points with names
pdf(paste0("~/labaled_net2.pdf")) 
plot(g,edge.width=edge.width,vertex.color=cols, edge.color= "#d3d3d3",layout=lo, edge.arrow.size=0.01, 
     vertex.shape='circle', vertex.size=sizes, vertex.frame.color=NA)
dev.off()


#set colours by input parameter type

#colors:
# 1) lysate cytokines: "#fdac53"
# 2) plasma cytokines: "#9bb7d4"
# 3) T cell phenotype: "#b55a30"
# 4) TME immune cell freq: "#0072b5"
# 5) Marker expression on immune subsets: "#a0daa9"
# 6) SF7 predictions: "#e9897e"

#I manually annotate my data w/ colors and inport in. 
node_cols <- read_excel("node_col_df.xlsx") %>%
  column_to_rownames( var = "...1") %>%
  t(.)

cols = node_cols
frame.cols = node_cols

#visualize points with dataset specific colour
plot(g,edge.width=edge.width,vertex.color=cols,vertex.frame.color=frame.cols, edge.color= "#d3d3d3",layout=lo, edge.arrow.size=0.01, 
     vertex.shape='circle', vertex.size=sizes, vertex.frame.color=NA, vertex.label=NA)

#make legend
plot(NULL)
legend("bottomleft", 
       legend = c("Tumor lysate cytokines", "Plasma cytokines", "T cell phenotype", "TME immune cell freq.", "Marker expr. immune subsets", "SF7-Fc predictions" ), 
       col = c("#fdac53", "#9bb7d4", "#b55a30", "#0072b5", "#a0daa9", "#e9897e"), 
       pch = c(19,19, 19, 19, 19), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1, 
       text.col = "black", 
       horiz = F )


#now color network by parameter indivdual AUROC
ind_aurocs <- as.data.frame(mypvs) %>%
  rownames_to_column(var = "identifier")
names_param <- as.data.frame(ROCvector) %>%
  rownames_to_column(var = "identifier")
named_aurocs <- left_join(ind_aurocs, names_param, by = "identifier") %>%
  select(mypvs, ROCvector)
named_aurocs <-  named_aurocs %>%
     rename(auroc = mypvs, features = ROCvector)
#save
write.csv(named_aurocs, "~/iEN_modeling/named_aurocs.csv")

#check dist of aurocs
histogram(named_aurocs$auroc)

#set colors
cols = vector()
i <- 0
for(i in 1:nrow(named_aurocs)){
  i=i+1
  if(named_aurocs$auroc[i] >= 0.9){
    cols[i] = "#092530"
  }else if(named_aurocs$auroc[i] <0.9 & named_aurocs$auroc[i] >= 0.775){
    cols[i] = "#214D5F"
  }else if(named_aurocs$auroc[i] <0.775 & named_aurocs$auroc[i] >= 0.6782){
    cols[i] =  "#30738D"
  }else if(named_aurocs$auroc[i] <0.6782 & named_aurocs$auroc[i] > 0.575){
    cols[i] =  "#5CA9C7"
  }else if(named_aurocs$auroc[i] <= 0.575){
    cols[i] = "#CCE4ED"
  }
}
cols <- as.data.frame(cols)
cols[1,] = "#5CA9C7" 

cols_clean <- cols %>%
  as.data.frame() %>%
  rownames_to_column(var = "identifier")

named_aurocs <-  named_aurocs %>%
  rownames_to_column(var = "identifier")

colored_aurocs <- left_join(cols_clean, named_aurocs, by = "identifier") %>%
  select(cols, features) %>%
  column_to_rownames(var = "features") %>%
  t(.)

#visualize points with colored by AUROC
cols = colored_aurocs
frame.cols = colored_aurocs

plot(g,edge.width=edge.width,vertex.color=cols,vertex.frame.color=frame.cols, edge.color= "#d3d3d3",layout=lo, edge.arrow.size=0.01, 
     vertex.shape='circle', vertex.size=sizes, vertex.frame.color=NA, vertex.label=NA)

#make legend
plot(NULL)
legend("bottomleft", 
       legend = c("AUROC > 0.9", "0.9 > AUROC > 0.775", "0.775 > AUROC > 0.678", "0.678 > AUROC > 0.575", "AUROC < 0.575" ), 
       col = c("#092530", "#214D5F", "#30738D", "#5CA9C7", "#CCE4ED"), 
       pch = c(19,19, 19, 19, 19), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1, 
       text.col = "black", 
       horiz = F )


#------------getting betas to see direction of change of each feature in model response
test_betas <- as.data.frame(as.matrix(ien.SF7_Fc_100@betas)) %>%
  summarise_all(median) %>%
  t(.) %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature")
write.csv(test_betas, "~/iEN_modeling/test_betas.csv")

summary(test_betas)
histogram(test_betas)
#-------------------------------


##--------now coloring network by iEN betas 

#set colors
#manually do this since it is giving me issues

#I manually annotate my data w/ colors and inport in. 
node_cols2 <- read_csv("test_betas.csv") %>%
  column_to_rownames(var = "feature") %>%
  select(cols) %>%
  t(.)

cols = node_cols2
frame.cols = node_cols2

plot(g,edge.width=edge.width,vertex.color=cols,vertex.frame.color=frame.cols, edge.color= "#d3d3d3",layout=lo, edge.arrow.size=0.01, 
     vertex.shape='circle', vertex.size=sizes, vertex.frame.color=NA, vertex.label=NA)

#make legend
plot(NULL)
legend("bottomleft", 
       legend = c("Increased in responders", "Decreased in responders", "Noncontributory" ), 
       col = c("#FC600A", "#347B98", "#e6e6e6"), 
       pch = c(19,19, 19), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1, 
       text.col = "black", 
       horiz = F )

#make legend for node size scaling
plot(NULL)
legend("bottomleft", 
       legend = c("", "", "") , 
       col = c("#d3d3d3", "#d3d3d3", "#d3d3d3"), 
       pch = c(19,19, 19), 
       bty = "n", 
       pt.cex = c(1,2,3), 
       horiz = F )

