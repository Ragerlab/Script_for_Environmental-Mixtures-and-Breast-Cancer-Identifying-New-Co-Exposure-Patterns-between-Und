**prioritization.R**- This script reads in

- The organized, clustered presence/absence dataframe **all_chems_w_cluster.csv** (*Experiment_2/output*)
- The name/CASRN/DTXSID mapping **dtxsid_name_ref.csv** (*Experiment_1/output*)
- The presently enriched structural features **present_enriched_chemotypes.csv** (*Experiment_3/output*)
- ToxPrint fingerprints **CompToxChemicalsDashboard_all_ToxPrint_030721_pt1.csv**, **CompToxChemicalsDashboard_all_ToxPrint_030721_pt2.csv**, and **toxprint_V2_vs_DTXSIDs_no_fingerprint_042221.csv** (*Experiment_3/output*)
-  OPERA predicitions **CompToxChemicalsDashboard_all_OPERA_032621_pt1.csv** and **CompToxChemicalsDashboard_all_OPERA_032621_pt2.csv** (*Experiment_4/input*)

The script then loops through the clusters identified of being of particular interest: 1, 4, 5, 6, and 9. A table is produced for each cluster that contains each UC in the cluster as a row and it's total number of presently enriched features and the highest spearman-rank correlation for that UC amongst all possible BC correlations in the cluster. The p-value of the correlation is included, as is the identifying information of the corresponding highest correlated BC. A structural similarity score and physicochemical similarity score, which are values calculated based on the distributions of enriched features and physicochemical correlations across the cluster, were also included. A final overall prioritization score is also calculated as the sum of the structural similarity score and the physicochemical similarity scores. The results for each cluster are arranged in desceding order of the overall prioritization score and are written to the *output* folder as 

-  cluster_1_prioritization.csv
-  cluster_4_prioritization.csv
-  cluster_5_prioritization.csv
-  cluster_6_prioritization.csv
-  cluster_9_prioritization.csv

**top_tens_heatmaps.R**- This script reads in the above-listed prioritization results, and the clustered, organized presence/absence data from the *output* of *Experiment_2*, as well as the name/CASRN/DTXSID mapping from the *output* of *Experiment_1*. For each of the clusters, 1,4,5,6, and 9, a heatmap is made from the presence/absence data of the 10 highest ranking UCs in the cluster by the overall score found in the prioritization results. The relevant, highest correlated BCs are also shown in the heatmap. The heatmaps are saved in teh *figures* folder as

- cluster_1_top_ten.png
- cluster_4_top_ten.png
- cluster_5_top_ten.png
- cluster_6_top_ten.png
- cluster_9_top_ten.png
