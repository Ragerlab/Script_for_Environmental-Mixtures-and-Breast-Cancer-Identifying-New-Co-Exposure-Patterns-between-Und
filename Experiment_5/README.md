**cluster_of_interest_identification**- This script reads in

- The organized, clustered presence/absence dataframe **all_chems_w_cluster.csv** (*Experiment_2/output*)
- The name/CASRN/DTXSID mapping **dtxsid_name_ref.csv** (*Experiment_1/output*)
- The presently enriched structural features **present_enriched_chemotypes.csv** (*Experiment_3/output*)
- ToxPrint fingerprints **CompToxChemicalsDashboard_all_ToxPrint_030721_pt1.csv**, **CompToxChemicalsDashboard_all_ToxPrint_030721_pt2.csv**, and **toxprint_V2_vs_DTXSIDs_no_fingerprint_042221.csv** (*Experiment_3/output*)
-  OPERA predicitions **CompToxChemicalsDashboard_all_OPERA_032621_pt1.csv** and **CompToxChemicalsDashboard_all_OPERA_032621_pt2.csv** (*Experiment_4/input*)
-  The correlations of every UC-BC pair in a cluster **all_correlations.csv** (*Experiment_4/output*)

The script loops through each of the 19 clusters and makes a summary table that includes the total number of chemicals in the cluster, the number of BC, NBCs, and UCs, the maximum number of presently enriched structural features in a UC in the cluster, and the maximum spearman-rank correlation of any UC-BC pair in the cluster. The table is exported to the *output* folder as **potential_clusters_of_interest.csv**.
