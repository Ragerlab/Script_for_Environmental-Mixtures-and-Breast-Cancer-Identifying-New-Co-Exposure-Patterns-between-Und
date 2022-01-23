**prioritization.R**- This script reads in

- The organized, clustered presence/absence dataframe **all_chems_w_cluster.csv** (*Experiment_2/output*)
- The name/CASRN/DTXSID mapping **dtxsid_name_ref.csv** (*Experiment_1/output*)
- The presently enriched structural features **present_enriched_chemotypes.csv** (*Experiment_3/output*)
- ToxPrint fingerprints **CompToxChemicalsDashboard_all_ToxPrint_030721_pt1.csv**, **CompToxChemicalsDashboard_all_ToxPrint_030721_pt2.csv**, and **toxprint_V2_vs_DTXSIDs_no_fingerprint_042221.csv** (*Experiment_3/output*)
-  OPERA predicitions **CompToxChemicalsDashboard_all_OPERA_032621_pt1.csv** and **CompToxChemicalsDashboard_all_OPERA_032621_pt2.csv** (*Experiment_4/input*)

The script then loops through the clusters identified of being of particular interest: 1, 4, 5, 6, 7, and 9. A table is produced for each cluster that contains each UC in the cluster as a row and it's total number of presently enriched features and the highest spearman-rank correlation for that UC amongst all possible BC correlations in the cluster. The p-value of the correlation is included, as is the identifying information of the corresponding highest correlated BC. A structural similarity score and physicochemical similarity score, which are values calculated based on the distributions of enriched features and physicochemical correlations across the cluster, were also included. A final prioritization score is also calculated as the sum of the structural similarity score and the physicochemical similarity scores. The results for each cluster are written to the *output* folder as 

-  cluster_1_prioritization.csv
-  cluster_4_prioritization.csv
-  cluster_5_prioritization.csv
-  cluster_6_prioritization.csv
-  cluster_7_prioritization.csv
-  cluster_9_prioritization.csv
