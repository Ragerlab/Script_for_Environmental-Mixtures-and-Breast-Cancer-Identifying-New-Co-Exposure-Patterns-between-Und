rm(list=ls())

library(cluster)
library(vegan)
library(gplots)
library(factoextra)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(janitor)
library(xlsx)
library(gridExtra)


setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/ChemExpoDB_Breast_Cancer/Experiments/sensitivity_analysis")


#read in data
chems<-read.csv(paste(getwd(),"Experiment_1/output/presence_absence_binary_df.csv",sep="/"))

dtxsid_name <- read.csv(paste(getwd(),"Experiment_1/output/dtxsid_name_cas_ref.csv", sep="/"))


#make breast cancer status (bc, nbc, uc) a factor
chems$status <- as.factor(chems$status)


#store the chemical names, CASRNs, and breast cancer status before removing these fields for clustering
chems <- merge(chems, dtxsid_name, by="DTXSID", all.x=TRUE)
names <- chems$true_chemname
cas <- chems$true_cas
stat <- chems$status

chems_reduced <- chems %>% column_to_rownames("DTXSID") %>% select(!c(count,status, true_chemname,true_cas))


#distance matrix to cluster chemicals
D_chem <- vegdist(as.matrix(chems_reduced),method = "jaccard")

#Run clustering algorithm
Sys.time()
cluster_chems <-diana(D_chem, diss=TRUE)
Sys.time()

plot(cluster_chems)


#Select 19 clusters as determined by figures from scripts plotting the wss and silhouette for 1<=n<=100 clusters
ncut <- 19
cluster_assignments_chem <- cutree(cluster_chems, k = ncut)

# number of chems in each cluster
k_chems <- table(cluster_assignments_chem)


#Add back in breast cancer status, chemical name, and CASRN as well as cluster assignment. Arrange the chemicals by cluster.
#Also add an index term so we can find where to add breaks between clusters in the heatmap.
chems_reduced <- chems_reduced %>% rownames_to_column("DTXSID")
chems_reduced$status <- stat
chems_reduced$true_chemname <- names
chems_reduced$true_cas <- cas
chems_reduced$cluster <- cluster_assignments_chem
chems_reduced <- chems_reduced %>% arrange(cluster)
chems_reduced$index <- 1:nrow(chems_reduced)



#transpose initial presence/absence df so we can also run clustering on exposure source categories
cT <- chems %>% column_to_rownames("DTXSID") %>% select(!c(count,status, true_chemname,true_cas))
cT <- t(cT)
cT <- as.data.frame(cT)

#pull list of exposure source categories
escs <- rownames(cT)

#create distance matrix from transposed df to cluster 
D_esc <- vegdist(as.matrix(cT),method = "jaccard")


#look for optimal number of clusters of exposure source categories
wss <- fviz_nbclust(cT, diss = D_esc, method = "wss", FUN=hcut, hc_func="diana", k.max=31)
wssplot <- wss$data
wssplot <- wssplot %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(x="number of clusters k", y="total within sum of squares")
ggsave("Experiment_2/figures/optimal_esc_clusters_wss.png", height=10, width=16, unit="in",plot=wssplot)


silho <- fviz_nbclust(cT, diss = D_esc, method = "silhouette", FUN=hcut, hc_func="diana", k.max=31)
silhoplot <- silho$data
silhoplot <- silhoplot %>% ggplot(aes(x=clusters, y=y))+geom_point()+geom_line(group=1)+labs(x="number of clusters k", y="average silhouette width")
ggsave("Experiment_2/figures/optimal_esc_clusters_silhouette.png", height=10, width=16, unit="in",plot=silhoplot)

#run clustering algorithm
cluster_escs<-diana(D_esc, diss=TRUE)
plot(cluster_escs)

# Select 12 clusters of exposure source categories
ncut <- 12
cluster_assignments_esc <- cutree(cluster_escs, k = ncut)

# number of exposures source categories in each cluster
k_escs <- table(cluster_assignments_esc)

#Add in the cluster assignment for the exposure source categories and then arrange by cluster. Also
#add in an index term so we can find where to add breaks between clusters on the heatmap
cT$cluster <- cluster_assignments_esc
cT <- cT %>% arrange(cluster)
cT$index <- 1:nrow(cT)


#reorganize columns of the arranged chemicals to reflect clustered exposure source categories
new_col_order <- rownames(cT)
new_col_order <- c(c("DTXSID", "true_cas","true_chemname","status","cluster", "index"),new_col_order)
chems_reduced <- chems_reduced[new_col_order]


### Summary of Chemical Clusters by Cancer Status ##########################################################################
# write a csv containing all information and cluster assignments for each chemical in the analysis
write.csv(chems_reduced %>% select(!index),"Experiment_2/output/all_chems_w_cluster.csv", row.names = FALSE)

#Make a table summarizing the number of chemicals in each cluster and the status (bc, nbc, uc) breakdown of the chemicals
cluster_counts <- chems_reduced %>% group_by(cluster) %>% dplyr::count(status, .drop=FALSE)
bc_chems <- cluster_counts %>% filter(status=="bc")
nbc_chems <- cluster_counts %>% filter(status=="nbc")
uc_chems <- cluster_counts %>% filter(status=="uc")

cluster_summary <- chems_reduced %>% group_by(cluster) %>% dplyr::count(cluster, name = "total_chems")
cluster_summary$bc_chems <- bc_chems$n
cluster_summary$nbc_chems <- nbc_chems$n
cluster_summary$uc_chems <- uc_chems$n

cluster_summary <- cluster_summary[c("cluster","bc_chems","nbc_chems","uc_chems","total_chems")]

#write a csv of the summary table
write.csv(cluster_summary, file="Experiment_2/output/cluster_summary.csv", row.names = FALSE)


#Make lists of where breaks should be placed in the heatmap to separate clusters for both chemicals and exposure source categories
sep_chems <- chems_reduced %>% group_by(cluster) %>% slice_max(n=1, order_by=index)
sep_chems <- sep_chems$index

sep_escs <- cT %>% group_by(cluster) %>% slice_max(n=1, order_by=index)
sep_escs <- sep_escs$index


#Recode presence (originally indicated by 1) as 2 for bcs and 3 for nbcs to change colors in the heatmap
l <- c("DTXSID", "true_chemname","true_cas","status", "cluster", "index")
to_mutate <- setdiff(colnames(chems_reduced), l)

chems_reduced <- chems_reduced %>% mutate_at(.vars = all_of(to_mutate), list(~case_when(status=="uc"~.*1,
                                                                                       status=="bc"~.*2,
                                                                                       status=="nbc"~.*3)))

chems_reduced <- chems_reduced %>% column_to_rownames("DTXSID") %>% select(!c(cluster, true_chemname,true_cas,status, index))


#set colors for heatmap (absent=gray, bc=red, nbc=blue, uc=yellow)
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette){
  stopifnot(length(colors) == 6)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  ramp3 <- colorRampPalette(colors[5:6])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2,ramp3))
}


cutoff.distance <- 1  
colors <-makeColorRampPalette(c("gray", "gold1","gold1","indianred3","indianred3","darkblue"), cutoff.distance / max(chems_reduced),100)


#make master heatmap
pheatmap(as.matrix(chems_reduced), color=colors,
         cluster_rows=FALSE, cluster_cols = FALSE, show_rownames = FALSE,
         legend = FALSE, cellheight = 0.135, cellwidth = 18, fontsize_co1 = 5,  angle_col = 45, gaps_row = sep_chems, gaps_col = sep_escs,
         filename = "Experiment_2/figures/clustered_heatmap.png", height = 17, width = 16)
