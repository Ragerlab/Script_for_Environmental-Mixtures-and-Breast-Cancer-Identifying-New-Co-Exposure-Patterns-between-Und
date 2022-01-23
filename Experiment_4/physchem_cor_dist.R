library(tidyverse)
library(janitor)

rm(list=ls())
setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/ChemExpoDB_Breast_Cancer/Experiments/sensitivity_analysis")

#Read in chemicals with cluster assignments
all_chems <- read_csv("Experiment_2/output/all_chems_w_cluster.csv")
name_stat <- all_chems %>% select(DTXSID,true_cas,true_chemname,status)

#read in DTXSID/name mapping
dtxsid_name <- read_csv("Experiment_1/output/dtxsid_name_cas_ref.csv")


#Read in OPERA phys chem data pulled from CompTox Dashboard. 
pt1 <- read.csv("Experiment_4/input/CompToxChemicalsDashboard_all_OPERA_032621_pt1.csv")
pt2 <- read.csv("Experiment_4/input/CompToxChemicalsDashboard_all_OPERA_032621_pt2.csv")

opera <- rbind(pt1,pt2)

opera <- opera %>% filter(INPUT!="INPUT") %>% 
  filter(DTXSID!="-")



features <- colnames(opera %>% select(!c(INPUT,DTXSID,PREFERRED_NAME,FOUND_BY)))


pc_tab <- merge(name_stat,opera, by="DTXSID")
pc_tab <- pc_tab %>% select(!c(INPUT,FOUND_BY,PREFERRED_NAME)) %>% mutate_at(features, function(x) as.numeric(x)) %>% filter(complete.cases(.))
write_csv(pc_tab,"Experiment_4/output/final_opera.csv")



opera_names <- opera %>% select(DTXSID, PREFERRED_NAME)

bc <- all_chems %>% filter(status=="bc")
uc <- all_chems %>% filter(status=="uc")



#initialize vectors to store values
clus_num <- c()
uc_list <- c()
bc_list <- c()
sr_cor <- c()
sr_pval <- c()

c <- 1
i <- 3
for(i in unique(all_chems$cluster)){

  #filter for cluster of interest
  coi <- all_chems %>% filter(cluster==i)
  
  opera_coi <- opera %>% 
    filter(DTXSID %in% coi$DTXSID) %>% 
    mutate_at(features, function(x) as.numeric(x))
  
  opera_coi <- opera_coi %>%
    filter(complete.cases(.))
  
  if(nrow(opera_coi %>% filter(DTXSID%in% all_of(bc$DTXSID)))==0 | nrow(opera_coi %>% filter(DTXSID %in% uc$DTXSID))==0){
    print(paste0("No BCs and/or No UCs with OPERA predictions. Skipping cluster ",i))
    next
  }
  else{
    print(paste0("evaluating cluster ",i))
  }
  
  #scale cols
  opera_coi <- opera_coi %>% column_to_rownames("DTXSID") %>% select(all_of(features))
  opera_coi <- as.data.frame(scale(opera_coi))
  
  opera_bc <- opera_coi %>% rownames_to_column("DTXSID") %>% filter(DTXSID %in% all_of(bc$DTXSID))
  opera_uc <- opera_coi %>% rownames_to_column("DTXSID") %>% filter(DTXSID %in% all_of(uc$DTXSID))
  
  print(paste0("number of BCs: ",nrow(opera_bc)))
  print(paste0("number of UCs: ",nrow(opera_uc)))
  cat("\n")
  

  #Run pairwise spearman rank correlations for each possible pair of UCs and BCs
  
  for(j in 1:nrow(opera_uc)){
    
    uc_DTXSID <- opera_uc %>% dplyr::slice(j) %>% select(DTXSID) %>% unlist() %>% unname()
    uc_props <- opera_uc %>% dplyr::slice(j) %>% select(!DTXSID) %>% unlist() %>% unname()
    
    for(k in 1:nrow(opera_bc)){
      bc_DTXSID <- opera_bc %>% dplyr::slice(k) %>% select(DTXSID) %>% unlist() %>% unname()
      bc_props <- opera_bc %>% dplyr::slice(k) %>% select(!DTXSID) %>% unlist() %>% unname()
      sr_cor_test <- cor.test(uc_props,bc_props, alternative = "t", method="spearman")
      clus_num[c] <- i
      uc_list[c] <- uc_DTXSID
      bc_list[c] <- bc_DTXSID
      sr_cor[c] <- sr_cor_test$estimate
      sr_pval[c] <- sr_cor_test$p.value
      c <- c+1
    }
  }
  
}

res <- data.frame(cluster=clus_num, uc_DTXSID=uc_list, bc_DTXSID=bc_list, sr_cor=sr_cor, sr_pval=sr_pval)

write_csv(res, "Experiment_4/output/all_correlations.csv")

g <- ggplot(res, aes(x=sr_cor)) +
  geom_histogram(color="black", fill="grey", binwidth = 0.05)+
  xlab("Correlation Value across UC-BC Pairs")+
  ylab("Number of UCs")+
  theme(axis.title=element_text(size=16))

ggsave("Experiment_4/figures/physchem_cor_dist.png", device="png", width = 6, height=6)

