rm(list=ls())

library(tidyverse)
library(janitor)


#Read in chemicals with cluster assignments
all_chems <- read_csv("Experiment_2/output/all_chems_w_cluster.csv")

dtxsid_name <- read_csv("Experiment_1/output/dtxsid_name_cas_ref.csv")

#read in enriched chemotypes in BCs
enriched_ct_pos <- read_csv("Experiment_3/output/present_enriched_chemotypes.csv")
enriched_ct_pos <- unlist(enriched_ct_pos$chemotype)


#ToxPrint fingerprints from dashboard, split between 2 files due to the dashboard search limit
fp1_dash <- read.csv("Experiment_3/input/CompToxChemicalsDashboard_all_ToxPrint_030721_pt1.csv")
fp2_dash <- read.csv("Experiment_3/input/CompToxChemicalsDashboard_all_ToxPrint_030721_pt2.csv")

#Additional fingerprints acquired from KI
fp_ki <- read.csv("Experiment_3/input/toxprint_V2_vs_DTXSIDs_no_fingerprint_042221.csv")   
fp_ki <- fp_ki %>% rename(DTXSID=dtxsid)

#Join all fingerprints into single dataframe
fps <- rbind(fp1_dash, fp2_dash)
fps <- fps %>% select(!c(INPUT, PREFERRED_NAME))
fps <- rbind(fps,fp_ki)

#drop DTXSIDs that didn't match in the Dashboard (DTXSID returned as "-") or didn't have a fingerprint
fps <- fps %>% 
  filter(DTXSID!="-") %>% 
  mutate_at(., vars(-DTXSID), function(x) as.numeric(x)) %>% 
  filter(complete.cases(.))


#Read in OPERA phys chem data pulled from CompTox Dashboard. 
pt1 <- read.csv("Experiment_4/input/CompToxChemicalsDashboard_all_OPERA_032621_pt1.csv")
pt2 <- read.csv("Experiment_4/input/CompToxChemicalsDashboard_all_OPERA_032621_pt2.csv")

opera <- rbind(pt1,pt2)

opera <- opera %>% filter(INPUT!="INPUT") %>% 
  filter(DTXSID!="-")

features <- colnames(opera %>% select(!c(INPUT,DTXSID,PREFERRED_NAME,FOUND_BY)))

opera_names <- opera %>% select(DTXSID, PREFERRED_NAME)

#Read in clusters that have data and can be characterized
potential_clusters <- read_csv("Experiment_4/output/all_correlations.csv")


clus_num <- c()
bc_num <- c()
nbc_num <- c()
uc_num <- c()
total_chems <- c()
max_ef <- c()
max_cor <- c()

for(i in 1:length(unique(potential_clusters$cluster))){
  clus_num[i] <- unique(potential_clusters$cluster)[i]

  
  #filter for cluster of interest
  coi <- all_chems %>% filter(cluster==unique(potential_clusters$cluster)[i])
  
  bc_num[i] <- nrow(coi %>% filter(status=="bc"))
  nbc_num[i] <- nrow(coi %>% filter(status=="nbc"))
  uc_num[i] <- nrow(coi %>% filter(status=="uc"))
  total_chems[i] <- nrow(coi)
  
  #isolate UC chemicals in cluster of interest
  uc <- coi %>% filter(status=="uc")
  bc <- coi %>% filter(status=="bc")
  
  
  #Select only enriched chemotypes and tally up the number present for each chemical
  uc_fps <- fps %>% filter(DTXSID %in% uc$DTXSID) %>%
    select(c(DTXSID, all_of(enriched_ct_pos))) %>%
    adorn_totals(where = "col")
  
  
  #Add in the chemical name and sort by the number of present features
  uc_fps <- merge(dtxsid_name,uc_fps, on="DTXSID", all.y = TRUE)
  uc_fps <- uc_fps %>% arrange(desc(Total))
  max_tot <- max(uc_fps$Total)
  
  
  max_ef[i] <- max_tot
  
  
  
  # drop DTXSIDs that didn't match in the Dashboard (DTXSID returned as "-")
  opera_coi <- opera %>% 
    filter(DTXSID %in% coi$DTXSID) %>% 
    column_to_rownames("DTXSID") %>% 
    select(!c(INPUT,PREFERRED_NAME, FOUND_BY)) %>% 
    mutate_all(as.numeric) 
  

  opera_coi <- opera_coi %>%
    filter(complete.cases(.))
  
  #scale cols
  opera_coi <- as.data.frame(scale(opera_coi))
  opera_coi <- opera_coi %>% rownames_to_column("DTXSID")
  

  
  opera_bc <- opera_coi %>% filter(DTXSID %in% all_of(bc$DTXSID))
  opera_uc <- opera_coi %>% filter(DTXSID %in% all_of(uc$DTXSID))
  
  
  
  #initialize vectors to store values
  uc_list <- c()
  bc_list <- c()
  sr_cor <- c()
  sr_pval <- c()
  
  #Run pairwise spearman rank correlations for each possible pair of UCs and BCs
  c <- 1
  for(j in 1:nrow(opera_uc)){

    uc_DTXSID <- opera_uc %>% dplyr::slice(j) %>% select(DTXSID) %>% unlist() %>% unname()
    uc_props <- opera_uc %>% dplyr::slice(j) %>% select(!DTXSID) %>% unlist() %>% unname()
    
    for(k in 1:nrow(opera_bc)){
      bc_DTXSID <- opera_bc %>% dplyr::slice(k) %>% select(DTXSID) %>% unlist() %>% unname()
      bc_props <- opera_bc %>% dplyr::slice(k) %>% select(!DTXSID) %>% unlist() %>% unname()
      sr_cor_test <- cor.test(uc_props,bc_props, alternative = "t", method="spearman")
      uc_list[c] <- uc_DTXSID
      bc_list[c] <- bc_DTXSID
      sr_cor[c] <- sr_cor_test$estimate
      sr_pval[c] <- sr_cor_test$p.value
      c <- c+1
    }
  }
  
  ps_res <- data.frame(uc_DTXSID=uc_list, bc_DTXSID=bc_list, sr_cor=sr_cor, sr_pval=sr_pval)
  ps_res <- ps_res %>% group_by(uc_DTXSID) %>% slice_max(sr_cor)
  
  max_sr_cor <- max(ps_res$sr_cor)
  max_cor[i] <- max_sr_cor
  
}

results <- data.frame(cluster=clus_num, bc_chems=bc_num, nbc_chems=nbc_num, uc_chems=uc_num, total_chems=total_chems, max_enriched_feat=max_ef, max_cor=max_cor)

final_coi <- results %>% filter(bc_chems!=0 & max_enriched_feat>=4 & max_cor>=.8)


write_csv(final_coi, "Experiment_5/output/potential_clusters_of_interest.csv")
