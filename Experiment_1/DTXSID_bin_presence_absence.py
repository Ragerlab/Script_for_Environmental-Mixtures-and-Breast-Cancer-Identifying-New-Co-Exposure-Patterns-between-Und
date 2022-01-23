import numpy as np
import os
import pandas as pd

os.chdir("Experiment_1")

#read in DTXSID/keyword sets
kw_sets=pd.read_csv("output/factotum_listPresence_092320_updated_chemnames_020621.csv", dtype=str)

#make reference table of DTXSID/true chemical name/ true CASRN
name_dtxsid=kw_sets[["DTXSID","true_chemname","true_cas"]]
name_dtxsid.drop_duplicates(inplace=True)
name_dtxsid.to_csv("output/dtxsid_name_cas_ref.csv", index=False)


#read in Silent Spring cancer associated chemicals with CASRN-DTXSID mapping from CompTox Dashboard
ssi = pd.read_csv("input/SSI_BC_dashboard.csv")
ssi = ssi[["DTXSID"]]
ssi=ssi.loc[ssi.DTXSID.str.contains("-")==False]

#read in ToxRefDB cancer associated chemicals
toxref_pos=pd.read_csv("input/ToxRefDB_043021_cancer_related.csv")
toxref_pos=toxref_pos.rename(columns={"dsstox_substance_id":"DTXSID"})
toxref_pos=toxref_pos[["DTXSID"]]
toxref_pos.drop_duplicates(inplace=True)

#make full list of breast cancer associated DTXSIDs
bc=pd.concat([toxref_pos,ssi])
bc.drop_duplicates(inplace=True)
bc.reset_index(inplace=True,drop=True)
bc=list(bc.DTXSID)

#Read in non-breast cancer associated chemicals and make full list.
toxref_negs=pd.read_csv("input/ToxRefDB_043021_unique_tested_chemicals.csv")
toxref_negs=toxref_negs[["dsstox_substance_id","dsstox_matchwith_unique_cancer_chemicals"]]
toxref_negs=toxref_negs.loc[toxref_negs.dsstox_matchwith_unique_cancer_chemicals=="no"]
toxref_negs=toxref_negs.loc[toxref_negs.dsstox_substance_id.isin(bc)==False]

nbc=list(toxref_negs.dsstox_substance_id.unique())


#read in keyword to exposure source category (esc) mapping
esc=pd.read_csv("input/keyword_esc.csv", usecols=[0,1,4])


#list of keywords that, when present in a set, warrant the chemical record being dropped from the DTXSID/keyword dataset
kw_1=list(esc.loc[esc.Status==1].Keyword)
kw_1_search=list(x.split()[0] for x in kw_1)

#list of keywords that, when present in a set, warrant removal from the set, but the other keywords remain om the DTXSID/keyword dataset
kw_2=list(esc.loc[esc.Status==2].Keyword)

#dictionary of keywords to consider and the exposure source category they map to for analysis
kw2esc=esc.loc[esc.Status==0]
kw2esc=dict(zip(kw2esc.Keyword,kw2esc.Exposure_Source_Cat))



#KW Filter 1: Remove Chemical List Presence records that don't have a kw set
kw_sets=kw_sets.loc[pd.isnull(kw_sets.listSets)==False]
kw_sets=kw_sets.loc[pd.isnull(kw_sets.true_chemname)==False]
kw_sets=kw_sets[["DTXSID","listSets"]]
kw_sets.reset_index(inplace=True, drop=True)


#KW Filter 2: Remove Chemical List Presence records that contain keywords in kw_1.
kw_sets=kw_sets.loc[kw_sets.listSets.str.contains("|".join(kw_1_search))==False]


#KW Filter 3: Remove instances of keywords in kw_2
kw_sets.loc[kw_sets.listSets.str.contains(";"),["listSets"]]=kw_sets.listSets.str.split(";")
kw_sets=kw_sets.rename(columns={"listSets":"keyword"}).explode("keyword")
kw_sets.keyword=kw_sets.keyword.str.strip()
kw_sets=kw_sets.loc[kw_sets.keyword.isin(kw_2)==False]



#Map keywords to exposure source categories.
chem_esc=kw_sets.copy()
chem_esc["esc"]=chem_esc.keyword.replace(kw2esc)

chem_esc=chem_esc[["DTXSID","esc"]]
chem_esc.drop_duplicates(inplace=True)


#Create binary df of DTXSIDs (rows) and exposure source categories (cols) based on whether the DTXSID is associated with the esc (1) or not (0)
pres_abs=pd.pivot_table(chem_esc,index="DTXSID",columns="esc",aggfunc=np.count_nonzero, fill_value=0).astype(bool)*1
pres_abs.reset_index(inplace=True)
pres_abs.columns.name=None
pres_abs.columns=[x.strip(" ") for x in pres_abs.columns]


#add count column with number of associated escs for each DTXSID and select DTXSIDs with at least 2 associated escs
pres_abs["count"]=np.sum(pres_abs,axis=1)
pres_abs= pres_abs.loc[pres_abs["count"]>=2]


#set breast cancer status as breast-cancer associated (bc), non-breast cancer associated (nbc), or untested (uc)
pres_abs["status"]="uc"
pres_abs.loc[pres_abs.DTXSID.isin(bc), ["status"]]="bc"
pres_abs.loc[pres_abs.DTXSID.isin(nbc), ["status"]]="nbc"


# write presence/absence dataframe to csv
pres_abs.to_csv("output/presence_absence_binary_df.csv", index=False)


#Table for SI
si_tab=chem_esc.loc[chem_esc.DTXSID.isin(pres_abs.DTXSID)]
si_tab=si_tab.groupby("DTXSID")["esc"].apply("; ".join).reset_index()
si_tab=si_tab.merge(name_dtxsid, on="DTXSID", how="left")

si_tab=si_tab[["DTXSID","true_chemname","true_cas","esc"]]

si_tab.to_csv("output/chems_for_analysis_esc_map.csv", index=False)
