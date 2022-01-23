library(tidyverse)
library(janitor)

rm(list=ls())
setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/ChemExpoDB_Breast_Cancer/Experiments/sensitivity_analysis")

#Read in chemicals with cluster assignments
all_chems <- read_csv("Experiment_2/output/all_chems_w_cluster.csv")

dtxsid_name <- read_csv("Experiment_1/output/dtxsid_name_cas_ref.csv")

#read in enriched chemotypes in BCs
enriched_ct_pos <- read_csv("Experiment_3/output/present_enriched_chemotypes.csv")
enriched_ct_pos <- unlist(enriched_ct_pos$chemotype)


#ToxPrint fingerprints from dashboard, split between 2 files due to the dashboard search limit
fp1_dash <- read.csv("Experiment_3/input/CompToxChemicalsDashboard_all_ToxPrint_030721_pt1.csv")
fp2_dash <- read.csv("Experiment_3/input/")

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

#isolate UC chemicals 
uc <- all_chems %>% filter(status=="uc")

#Select only present enriched chemotypes and tally up the number present for each chemical
uc_fps_pres <- fps %>% filter(DTXSID %in% uc$DTXSID) %>%
  select(c(DTXSID, all_of(enriched_ct_pos))) %>%
  adorn_totals(where = "col")


#Select only absent enriched chemotypes and tally up the number absent for each chemical
# uc_fps_abs <- fps %>% filter(DTXSID %in% uc$DTXSID) %>%
#   select(c(DTXSID, all_of(enriched_ct_neg))) %>%
#   adorn_totals(where = "col") %>% 
#   mutate(missing=length(enriched_ct_neg)-Total)

#Add in the chemical name and sort by the number of present features
uc_fps_pres <- merge(dtxsid_name,uc_fps_pres, on="DTXSID", all.y = TRUE)
uc_fps_pres <- uc_fps_pres %>% arrange(desc(Total))
uc_fps_pres$Total <- as.factor(uc_fps_pres$Total)


#log10 y axis
g <- ggplot(uc_fps_pres, aes(x=Total)) +
  geom_bar() +
  scale_y_log10() +
  xlab("Number of enriched features")+
  ylab("Number of UCs")+
  theme(axis.title=element_text(size=16))

ggsave("Experiment_4/figures/present_enriched_chemotypes_dist.png", device="png", width = 6, height = 6)
