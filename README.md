# Script_for_Environmental-Mixtures-and-Breast-Cancer-Identifying-New-Co-Exposure-Patterns-between-Und
Script for 'Environmental Mixtures and Breast Cancer: Identifying New Co-Exposure Patterns between Understudied vs Breast Cancer-Associated Chemicals using Chemical Inventory Informatics'.

This study aimed to to identify environmental chemicals in use inventories that co-occur and share properties with chemicals that have association with breast cancer, highlighting exposure combinations that may alter disease risk. The analysis was broken down into 6 distinct *in silico* experiments that sequentially build off each other. An overview of each experiment is supplied below, with further detail on each experiment provided in the corresponding folder above.


**Expermient_1: Data preparation and organization**- This analysis organizes the chemical dataset then bins the chemicals into three categories related to their association with breast cancer (breast cancer associated, not breast cancer associated, and untested in relation to breast cancer). Chemicals are then mapped to exposure data and the dataset is prepped for further analysis.

**Experiment_2: Clustering Analysis**- This analysis groups the chemical data into clusters by use patterns derived from the exposure data. Hierearchical clustering is employed to generate these groupings, and a heatmap is produced to show the clusters and global landscape of the chemicals in the analysis.

**Experiment_3: Structural Feature Enrichment**- This analysis identifies chemical structural features from ToxPrint fingerprints that are more commonly present, followed by more commonly absent, in chemicals associated with breast cancer, when compared to chemicals with no known association to breast cancer. These were deemed enriched features.

**Experiment_4: Distributions of Enriched Features and Correlations**- This analysis first examines chemicals that have been untested in relation to breast cancer for enriched features that were commonly present in breast cancer associated chemicals. Then a correlation test is run on the physicochemical properties, in the form of OPERA predictions, for every possible pair of untested chemicals and breast cancer associated chemicals within an individual cluster. Plots are generated showing the distribution of the number of enriched features in untested chemicals and the correlation values between the untested chemicals and breast cancer associated chemicals.

**Experiment_5: Cluster of Interest Identification**- This analysis examines each cluster of chemicals for the number of untested chemicals, the number of enriched features in untested chemicals, and the highest correlation value for physicochemical properties between an untested chemical and breast cancer associated chemical in the cluster. A summary table is produced, that when in conjunction with visual inspection of the heatmap from **Experiment_2**, informed which clusters of chemicals should be further evaluated. 

**Experiment_6: Prioritization and Visualization of Clusters of Interest**- This analysis ranks the untested chemicals within each of cluster 1, 4, 5, 6, and 9. This ranking is based on a structural similarity score and physicochemical similiarity score which are derived for each untested chemical by comparing the number of enriched features and correlation values, respectively, to all other untested chemicals in the corresponding cluster. Heatmaps are made for the top ten ranking untested chemicals in each cluster showing use patterns, along with co-occuring breast cancer associated chemicals.
<br>


An overview of exposure data curation:
![Study_Overview_Fig](https://user-images.githubusercontent.com/72747901/146388561-7cd49881-55c1-42b1-8065-e3d02796379a.png)


