rm(list=ls())

#load libraries
library(cluster)
library(vegan)
library(gplots)
library(factoextra)
library(tidyverse)


#read in presence/absence dataframe and name reference from output of experiment_1
chems <- read.csv(paste(getwd(),"Experiment_1/output/presence_absence_binary_df.csv",sep="/"))

#prepare dataframe for evaluation
chems_reduced <- chems %>% column_to_rownames("DTXSID") %>% select(!c(count,status))

#generate distance matrix
D <- vegdist(as.matrix(chems_reduced),method = "jaccard")

#determine optimal clusters using silhouette metric
Sys.time()
silho <- fviz_nbclust(chems_reduced, diss = D, method = "silhouette", FUN=hcut, hc_func="diana", k.max=100)
Sys.time()

#save figure 
silhoplot <- silho$data
write_csv(silhoplot, "output/silhouette_data.csv")
silhoplot <- silhoplot %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(title = "Optimal Number of Clusters: Silhouette", x="number of clusters k", y="average silhouette width")
ggsave("figures/optimal_chemical_clusters_silhouette.png", height=10, width=16, unit="in",plot=silhoplot)
