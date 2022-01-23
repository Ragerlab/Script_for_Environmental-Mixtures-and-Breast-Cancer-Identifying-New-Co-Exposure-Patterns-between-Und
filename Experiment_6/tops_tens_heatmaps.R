rm(list=ls())

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(janitor)


#read in data
all_chems <- read_csv("Experiment_2/output/all_chems_w_cluster.csv")

dtxsid_name <- read_csv("Experiment_1/output/dtxsid_name_cas_ref.csv")

#make breast cancer status (bc, nbc, uc) a factor
all_chems$status <- as.factor(all_chems$status)

#read in top 10 UCs for each cluster of interest
c1 <- read_csv("Experiment_6/output/cluster_1_prioritization.csv")
c1 <- c1 %>% slice_head(n=10)

c4 <- read_csv("Experiment_6/output/cluster_4_prioritization.csv")
c4 <- c4 %>% slice_head(n=10)

c5 <- read_csv("Experiment_6/output/cluster_5_prioritization.csv")
c5 <- c5 %>% slice_head(n=12)

c6 <- read_csv("Experiment_6/output/cluster_6_prioritization.csv")
c6 <- c6 %>% slice_head(n=11)

c7 <- read_csv("Experiment_6/output/cluster_7_prioritization.csv")
c7 <- c7 %>% slice_head(n=10)

c9 <- read_csv("Experiment_6/output/cluster_9_prioritization.csv")
c9 <- c9 %>% slice_head(n=10)


#Recode presence (1)  for bc as 2 and for nbc as 3 to change colors in heatmap by status
l <- c("DTXSID", "true_chemname","true_cas","status","cluster")
to_mutate <- setdiff(colnames(all_chems), l)

all_chems <- all_chems %>% mutate_at(.vars = all_of(to_mutate), list(~case_when(status=="uc"~.*1, status=="bc"~.*2, status=="nbc"~.*3))) 

sep_cols <- c(4,14,16,19,21,23,25,26, 28,29,31)

#set colors for heatmap (absent=gray, bc=red, nbc=blue, uc=yellow)
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette){
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}


cutoff.distance <- 1  

escs <- all_chems %>% select(to_mutate)
colors <-makeColorRampPalette(c("gray", "gold1","gold1","indianred3"), cutoff.distance / max(escs),100)


#cluster heatmaps
chems1 <- all_chems %>% filter(all_chems$DTXSID %in% c1$DTXSID | all_chems$DTXSID %in% c1$bc_DTXSID) %>% select(!c("DTXSID","status","true_cas","cluster")) %>% column_to_rownames("true_chemname")

pheatmap(as.matrix(chems1), main="Cluster 1", color=colors,
         cluster_rows=FALSE, cluster_cols = FALSE, legend = FALSE, cellheight = 12,
         cellwidth = 16, fontsize_co1 = 5,  angle_col = 45, gaps_col = sep_cols,filename = "Experiment_6/figures/cluster_1_top_ten.png",
         height = 10, width = 18)



chems4 <- all_chems %>% filter(all_chems$DTXSID %in% c4$DTXSID | all_chems$DTXSID %in% c4$bc_DTXSID) %>% select(!c("DTXSID","status","true_cas","cluster")) %>% column_to_rownames("true_chemname")

pheatmap(as.matrix(chems4), main="Cluster 4", color=colors,
         cluster_rows=FALSE, cluster_cols = FALSE, legend = FALSE, cellheight = 12,
         cellwidth = 16, fontsize_co1 = 5,  angle_col = 45, gaps_col = sep_cols,filename = "Experiment_6/figures/cluster_4_top_ten.png",
         height = 10, width = 22)


chems5 <- all_chems %>% filter(all_chems$DTXSID %in% c5$DTXSID | all_chems$DTXSID %in% c5$bc_DTXSID) %>% select(!c("DTXSID","status","true_cas","cluster")) %>% column_to_rownames("true_chemname")

pheatmap(as.matrix(chems5), main="Cluster 5", color=colors,
         cluster_rows=FALSE, cluster_cols = FALSE, legend = FALSE, cellheight = 12,
         cellwidth = 16, fontsize_co1 = 5,  angle_col = 45, gaps_col = sep_cols,filename = "Experiment_6/figures/cluster_5_top_ten.png",
         height = 10, width = 16)


chems6 <- all_chems %>% filter(all_chems$DTXSID %in% c6$DTXSID | all_chems$DTXSID %in% c6$bc_DTXSID) %>% select(!c("DTXSID","status","true_cas","cluster")) %>% column_to_rownames("true_chemname")

pheatmap(as.matrix(chems6), main="Cluster 6", color=colors,
         cluster_rows=FALSE, cluster_cols = FALSE, legend = FALSE, cellheight = 12,
         cellwidth = 16, fontsize_co1 = 5,  angle_col = 45, gaps_col = sep_cols,filename = "Experiment_6/figures/cluster_6_top_ten.png",
         height = 10, width = 18)



chems7 <- all_chems %>% filter(all_chems$DTXSID %in% c7$DTXSID | all_chems$DTXSID %in% c7$bc_DTXSID) %>% select(!c("DTXSID","status","true_cas","cluster")) %>% column_to_rownames("true_chemname")

pheatmap(as.matrix(chems7), main="Cluster 7", color=colors,
         cluster_rows=FALSE, cluster_cols = FALSE, legend = FALSE, cellheight = 12,
         cellwidth = 16, fontsize_co1 = 5,  angle_col = 45, gaps_col = sep_cols,filename = "Experiment_6/figures/cluster_7_top_ten.png",
         height = 10, width = 18)


chems9 <- all_chems %>% filter(all_chems$DTXSID %in% c9$DTXSID | all_chems$DTXSID %in% c9$bc_DTXSID) %>% select(!c("DTXSID","status","true_cas","cluster")) %>% column_to_rownames("true_chemname")

pheatmap(as.matrix(chems9), main="Cluster 9", color=colors,
         cluster_rows=FALSE, cluster_cols = FALSE, legend = FALSE, cellheight = 12,
         cellwidth = 16, fontsize_co1 = 5,  angle_col = 45, gaps_col = sep_cols, filename = "Experiment_6/figures/cluster_9_top_ten.png",
         height = 10, width = 16)


