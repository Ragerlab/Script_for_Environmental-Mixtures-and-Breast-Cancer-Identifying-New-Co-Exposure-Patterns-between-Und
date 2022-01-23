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


#selected clusters of interest
coi <- c(1,4,5,6,7,9)

i <- 1
#loop through and prioritize each of the clusters of interest
for(i in coi){
  
  print(i)
  
  #filter for cluster of interest
  clust <- all_chems %>% filter(cluster==i)
  
  #isolate UC chemicals in cluster of interest
  uc <- clust %>% filter(status=="uc")

  ############### ToxPrint Prioritization ####################################################################
  
  #Select only enriched chemotypes and tally up the number present for each chemical
  uc_fps <- fps %>% filter(DTXSID %in% uc$DTXSID) %>%
    select(c(DTXSID, all_of(enriched_ct_pos))) %>%
    adorn_totals(where = "col")
  
  #Add in the chemical name and sort by the number of present features
  uc_fps <- merge(dtxsid_name,uc_fps, on="DTXSID", all.y = TRUE)
  uc_fps <- uc_fps %>% arrange(desc(Total))
  
  #Add Structural Similarity score
  max_tot <- max(uc_fps$Total)
  min_tot <- min(uc_fps$Total)
  dif <- max_tot-min_tot
  
  ss_res <- uc_fps %>% mutate(SS_score=(Total-min_tot)/dif) %>% select(DTXSID,true_cas, true_chemname, Total, SS_score)
  
  ############### OPERA Prioritization #######################################################################
  
  #Filter for DTXSIDs in the cluster of interest 
  opera_clust <- opera  %>% 
    filter(DTXSID %in% clust$DTXSID) %>% 
    column_to_rownames("DTXSID") %>% 
    select(!c(INPUT,PREFERRED_NAME, FOUND_BY)) %>% 
    mutate_all(as.numeric) 
  
  #Drop DTXSIDs that do not have OPERA predictions in the Dashboard.
  opera_clust <- opera_clust %>%
    filter(complete.cases(.))
  
  #scale OPERA features
  opera_clust <- as.data.frame(scale(opera_clust))
  opera_clust <- opera_clust %>% rownames_to_column("DTXSID")
  
  #isolate OPERA predictions for BC and UC chemicals in cluster of interest
  bc_chems <- clust %>% filter(status=="bc") %>%
    select(DTXSID) %>%
    unlist() %>%
    unname()
  
  uc_chems <- clust %>% filter(status=="uc") %>%
    select(DTXSID) %>%
    unlist() %>%
    unname()
  
  opera_bc <- opera_clust %>% filter(DTXSID %in% all_of(bc_chems))
  opera_uc <- opera_clust %>% filter(DTXSID %in% all_of(uc_chems))
  
  #initialize vectors to store values
  uc_list <- c()
  bc_list <- c()
  sr_cor <- c()
  sr_pval <- c()
  
  #Run pairwise spearman-rank correlations for each possible pair of UCs and BCs
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
  
  #Organize correlation results and select the highest correlation for each UC
  ps_res <- data.frame(uc_DTXSID=uc_list, bc_DTXSID=bc_list, sr_cor=sr_cor, sr_pval=sr_pval)
  ps_res<- ps_res %>% group_by(uc_DTXSID) %>% slice_max(sr_cor)
  
  #Add in chemical names
  ps_res <- merge(dtxsid_name,ps_res, by.x="DTXSID", by.y="bc_DTXSID", all.y = TRUE)
  ps_res <- ps_res %>% rename(bc_DTXSID=DTXSID) %>% rename(bc_chemname=true_chemname) %>% rename(bc_cas=true_cas)
  
  ps_res <- merge(dtxsid_name,ps_res, by.x="DTXSID", by.y="uc_DTXSID", all.y = TRUE)
  ps_res <- ps_res %>% rename(uc_DTXSID=DTXSID) %>% rename(uc_chemname=true_chemname) %>% rename(uc_cas=true_cas)
  

  #Add Physicochemical similarity score
  max_cor <- max(ps_res$sr_cor)
  min_cor <- min(ps_res$sr_cor)
  dif <- max_cor-min_cor
  
  
  ps_res <- ps_res %>% mutate(PS_score=(sr_cor-min_cor)/dif)
  ps_res <- ps_res %>% arrange(desc(PS_score))
  
  
  ############### Final Prioritization #######################################################################

  # make sure all UCs have a chemical name, combine SS and PS prioritizations, and sum SS_score and PS_score to get final prioritization score
  results <- merge(ss_res,ps_res,by.x="DTXSID", by.y = "uc_DTXSID", all=TRUE)
  results <- results %>%
    rename(ss_chemname=true_chemname) %>% 
    rename(ss_cas=true_cas) %>% 
    mutate(prioritization_score=PS_score+SS_score) %>%
    rename(total_enrich_feat=Total) %>% 
    arrange(desc(prioritization_score))
  
  #export results as a csv
  write_csv(results, paste0("Experiment_6/output/cluster_",i,'_prioritization.csv'))
}
