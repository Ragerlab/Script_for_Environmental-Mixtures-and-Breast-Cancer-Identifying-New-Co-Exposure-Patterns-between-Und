rm(list=ls())

library(janitor)
library(tidyverse)

setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/ChemExpoDB_Breast_Cancer/Experiments/sensitivity_analysis")

#ToxPrint fingerprints from dashboard, split between 2 files due to the dashboard search limit
fp1_dash <- read.csv("Experiment_3/input/CompToxChemicalsDashboard_all_ToxPrint_030721_pt1.csv")
fp2_dash <- read.csv("Experiment_3/input/CompToxChemicalsDashboard_all_ToxPrint_030721_pt2.csv")

#Additional fingerprints acquired from KI
fp_ki <- read.csv("Experiment_3/input/toxprint_V2_vs_DTXSIDs_no_fingerprint_042221.csv")   
fp_ki <- fp_ki %>% dplyr::rename(DTXSID=dtxsid)

#Join all fingerprints into single dataframe
fps <- rbind(fp1_dash, fp2_dash)
fps <- fps %>% select(!c(INPUT, PREFERRED_NAME))
fps <- rbind(fps,fp_ki)

#drop DTXSIDs that didn't match in the Dashboard (DTXSID returned as "-") or didn't have a fingerprint
fps <- fps %>% 
  filter(DTXSID!="-") %>% 
  mutate_at(., vars(-DTXSID), function(x) as.numeric(x)) %>% 
  filter(complete.cases(.))

#make lists of breast cancer positive (bc) and negative (nbc) chemicals
chems <- read_csv("Experiment_2/output/all_chems_w_cluster.csv")
names <- chems %>% select(DTXSID,true_cas, true_chemname, status)

bc_chems <- chems %>% filter(status=="bc") %>% select(DTXSID) %>% unlist()
nbc_chems <- chems %>% filter(status=="nbc") %>% select(DTXSID) %>% unlist()


fps_tab <- merge(names,fps, by="DTXSID")
fps_tab <- fps_tab %>% filter(DTXSID %in% all_of(c(bc_chems,nbc_chems)))
write_csv(fps_tab,"Experiment_3/output/bc_nbc_fps.csv")

chemotypes <- colnames(fps)[2:length(fps)]

#make dataframe to store results
results_pos <- data.frame(chemotype=character(), odds_ratio=numeric(), fisher_pval=numeric(), TP=numeric())
results_neg <- data.frame(chemotype=character(), odds_ratio=numeric(), fisher_pval=numeric(), TP=numeric())

#loop through each feature and make a confusion matrix  where TP are bc chems where the feature is present (1), FP are bc chems
#where the feature is absent (0), TN are nbc chems where the feature is present, and FN are nbc chems where the feature is absent.
#From the confusion matrix, the odds ratio and one tailed fisher's exact p-value are calculated. This approach is described here:
#https://www.sciencedirect.com/science/article/pii/S0160412018321196?via%3Dihub#bb0065. Features are considered enriched if the
#OR>=3, pval<=0.05, and TP>=3.


for(i in 1:length(chemotypes)){
  ct <- chemotypes[i]

  data <- fps %>% select(DTXSID,ct)
  
  bc <- data %>% filter(DTXSID %in% bc_chems) %>% select(all_of(ct)) %>% unlist() %>% unname()
  nbc <- data %>% filter(DTXSID %in% nbc_chems) %>% select(all_of(ct)) %>% unlist() %>% unname()
  
  tab <- rbind(c(sum(bc),length(bc)-sum(bc)),c(sum(nbc),length(nbc)-sum(nbc)))
  colnames(tab) <- c("1","0")
  rownames(tab) <- c("bc","nbc")

  stats_pos <- fisher.test(tab, alternative="g")
  stats_neg <- fisher.test(tab, alternative="l")

  ct_res_pos <- unname(list(ct, stats_pos$estimate, stats_pos$p.value, tab[1,1]))
  results_pos <- rbind(results_pos, ct_res_pos)
  
  ct_res_neg <- unname(list(ct, stats_neg$estimate, stats_neg$p.value, tab[1,2]))
  results_neg <- rbind(results_neg, ct_res_neg)
}

colnames(results_pos) <- c("chemotype", "odds_ratio","fisher_pval","TP")
rownames(results_pos) <- seq(1,nrow(results_pos),1)

colnames(results_neg) <- c("chemotype", "odds_ratio","fisher_pval","FP")
rownames(results_neg) <- seq(1,nrow(results_neg),1)


#enriched features 
enriched_features_pos <- results_pos %>% filter(odds_ratio>=3 & fisher_pval<0.05 & TP>=3)
enriched_features_neg <- results_neg %>% filter(odds_ratio<=0.3 & fisher_pval<0.05 & FP>=3)


write_csv(enriched_features_pos, file="Experiment_3/output/present_enriched_chemotypes.csv")
write_csv(enriched_features_neg, file="Experiment_3/output/absent_enriched_chemotypes.csv")

