import pandas as pd
import os

os.chdir("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/ChemExpoDB_Breast_Cancer/Experiments/sensitivity_analysis/Experiment_1")

#read in CPDat Chemical List Presence dataset
cedb=pd.read_csv("input/all_factotum_DTXSID_listPresence_2020-09-23.csv")

#read in file of missing chemical name/DTXSID pairs acquired from CompTox Dashboard
dashboard=pd.read_csv("input/missing_chemnames_dashboard_020521.csv", usecols=[2,3,4], names=["DTXSID","true_chemname","true_cas"], header=0)

#create dataframe from records from Chemical List Presence dataset that have a chemical name
cedb_drop_none=cedb.loc[cedb.true_chemname!="None"]

#create dataframe from records from the Chemical List Presence dataset that don't have a chemical name
none=cedb.loc[cedb.true_chemname=="None"]
none=none[["DTXSID","docID","raw_cas","raw_chem_name","listSets"]]

#merge in chemical names from Dashboard to records missing chemical names by DTXSID
none_fix=none.merge(dashboard, how="left", on="DTXSID")
none_fix=none_fix[cedb_drop_none.columns]

#rejoin all records into a single dataframe
fixed=pd.concat([cedb_drop_none,none_fix], axis=0)

#Check that the updated dataframe is the same size as the original Chemical List Presence dataset then write it to a csv
if fixed.shape==cedb.shape:
    fixed.to_csv("output/factotum_listPresence_092320_updated_chemnames_020621.csv")
