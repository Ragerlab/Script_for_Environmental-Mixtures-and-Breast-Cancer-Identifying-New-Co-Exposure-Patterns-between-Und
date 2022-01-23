rm(list=ls())

#load libraries
library(cluster)
library(vegan)
library(gplots)
library(factoextra)
library(tidyverse)


#read in presence/absence dataframe and name reference from output of experiment_1
chems <- read.csv(paste(getwd(),"Experiment_1/output/presence_absence_binary_df.csv",sep="/"))
dtxsid_name <- read.csv(paste(getwd(),"Experiment_1/output/dtxsid_name_cas_ref.csv", sep="/"))

#prepare dataframe for evaluation
chems<- merge(chems, dtxsid_name, by="DTXSID", all.x = TRUE)

chems$status <- as.logical(chems$status)


names <- chems$true_chemname
ids <- chems$DTXSID
bc_stat <- chems$status

chems_reduced <- chems %>% select(!c(DTXSID,count,status, true_chemname, true_cas))


#generate distance matrix
D <- vegdist(as.matrix(chems_reduced),method = "jaccard")

#run WSS to determine optimal nuimber of clusters
Sys.time()
wss <- fviz_nbclust(chems_reduced, diss = D, method = "wss", FUN=hcut, hc_func="diana", k.max=100)
Sys.time()

#save figure 
wssplot <- wss$data
wssplot <- wssplot %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(title = "Optimal Number of Clusters: WSS", x="number of clusters k", y="total within sum of squares")
ggsave("wss_test.png", height=10, width=16, unit="in",plot=wssplot)



