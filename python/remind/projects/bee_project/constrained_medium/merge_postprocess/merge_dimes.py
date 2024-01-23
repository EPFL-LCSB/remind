import networkx as nx
from remind.utils.postprocessing import *
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from remind.utils.merge_files import *
import os
from pytfa.thermo.tmodel_struct import ThermoModelStructure
import os
from os.path import join
import cobra.test
from cobra.io import read_sbml_model
import time
from pytfa.optim.utils import symbol_sum
import pandas as pd
import numpy as np
from sys import argv
from itertools import combinations
import os

plt.rcParams.update({'font.size':25,
                      'font.family':'Calibri',
                     'legend.fontsize':16,
                   'xtick.labelsize' : 18,
                   'ytick.labelsize' : 18
  })
plt.rcParams['svg.fonttype'] = 'none'




# frame_removed=pd.concat(list_removed, ignore_index=True)

#find essential_alternates


# def get_essential_mets(frame,groupby=True):
#     a=frame.metabolites.value_counts()
#     if groupby:
#         max_alt_number=frame.groupby(['yield_perc']).alternative.nunique().sum()
#     else:
#         max_alt_number=frame.alternative.nunique()
#     essential_mets=a[a==max_alt_number].index.to_list()
#     return essential_mets
#
# frame_subs=frame[frame.BU>0.5]
# essential_subs=frame_subs.groupby(['model']).apply(lambda gr: get_essential_mets(gr,groupby=False) ).reset_index(name='essential_subs')
# essential_subs['len']=[len(k) for k in essential_subs.essential_subs]
# # frame.to_hdf('/remind/projects/bee_project/stats/combined_dimes_180522_bee.h5',key='s')
# biomass_rxn='Growth'
# essential_subs.to_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_subs_170822.h5',key='s')
#
# list_es=[]
# for k in essential_subs.essential_subs:
#     list_es+=k
#
#
# frame_prods=frame[frame.FU>0.5]
# essential_prods=frame_prods.groupby(['model','yield_perc']).apply(lambda gr: get_essential_mets(gr) ).reset_index(name='essential_subs')
# essential_prods['len']=[len(k) for k in essential_prods.essential_subs]
# # frame.to_hdf('/remind/projects/bee_project/stats/combined_dimes_180522_bee.h5',key='s')
# biomass_rxn='Growth'
#
# l_ess=[]
# for k in essential_subs.essential_subs:
#     l_ess+=k


# path="/remind/projects/bee_project/PRAY_LATEST_111023/"
path="/remind/projects/bee_project/data/DiMES_bee"
allFolders = os.listdir(path)
'get a removed list'

list_removed=[]
list_=[]
# allFolders=allFolders[1:]
for subfolder in allFolders:
    list_organism=[]
    path_new=path+'/'+subfolder+'/alternatives'
    allFiles=glob.glob(path_new + "/*"+"h5")
    for file_ in allFiles:
        df = pd.read_hdf(file_)
        df['model']=[subfolder]*df.shape[0]
        list_organism.append(df)
    frame = pd.concat(list_organism, ignore_index=True)
    frame_removed=remove_redundant_sols(frame,groups=['yield_perc','alternative'])
    list_removed.append(frame_removed)
    list_.append(frame)
    #list_removed.append(d)
frame = pd.concat(list_, ignore_index=True)
frame_removed=pd.concat(list_removed, ignore_index=True)


frame_removed.to_hdf('/remind/projects/bee_project/stats/combined_dimes_unique_111023_bee.h5',key='s')
frame.to_hdf('/remind/projects/bee_project/stats/combined_dimes_not_unique_111023_bee.h5',key='s')



#add_size(frame,group=['alternative','model','yield_perc'],name='alt_size')
#checkfor unique alternatives

frame=frame_removed
for k in frame.model.unique():
    f=get_subs_prods_per_group(frame.groupby('model').get_group(k),gr=['alternative','yield_perc'])
    print('For model {} unique alternatives found are {} total alternatives are {} '.format(k,f.groupby(['subs_list','prods_list']).ngroups,f.groupby('yield_perc').alternative.nunique().sum()))


frame_removed=pd.read_hdf('/remind/projects/bee_project/stats/combined_dimes_unique_180522_bee.h5')

f = frame_s.groupby(['yield_perc', 'alternative']).apply(
    lambda gr: tuple(gr[gr.flux < -1.001].metabolites.unique())).reset_index(name='c_source')
f['len'] = [len(k) for k in f.c_source]

import re
def calculate_carbon_influx_outflux(df_alternative,frame_all):

    df_alternative['formula'] = [frame_all[frame_all.mets == k.replace("EX_", "")].formula.iloc[0] for k in
                                     df_alternative.metabolites]

    df_alternative['formula_weight'] = [frame_all[frame_all.mets == k.replace("EX_", "")].formula_weight.iloc[0] for
                                            k in df_alternative.metabolites]
    expr=[]
    pattern = r"[C]\d{1,3}|[C][A-Z]"
    for m in df_alternative.formula:
        result = re.findall(pattern, m)
        # returns a list access the string by [0]
        # if the next  number is digit remove C and take the float numver
        if result[0][1].isdigit():
            # print('yes')
            C_number = result[0].replace('C', '')
            expr+=[float(C_number)]
        # 'this means it has H or other molecule after so C number is 1'
        else:
            expr += [1]
    df_alternative['c_number']=expr
    df_alternative['cflux']=df_alternative.c_number*df_alternative.flux
    df_cflux=df_alternative.groupby(['model','yield_perc','alternative']).cflux.sum().reset_index()

    return df_cflux




d=pd.read_hdf("/remind/projects/bee_project/PRAY_3_081023_010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_1346.0_pars_c_uptakes_only_C_moles_tol_pars_11.0_tol_lb_10.0.h5")


frame=frame_dime_all_b[frame_dime_all_b.yield_perc.round(1)==1.0]

frame=frame_b[frame_b.yield_perc.round(1)==0.2]

"""after giving all lb ubs correct for inorganics et al checking if the dimes produces growth"""
#trying to fix the dime

for alt in (frame.alternative.unique()):

    fr=frame[frame.alternative==alt]
    s,p=get_subs_prods(fr)

    for k in list_carbon_sources:
        if k in s:
            model.reactions.get_by_id(k).lower_bound=fr[fr.metabolites==k].flux.iloc[0]*1.01
            model.reactions.get_by_id(k).upper_bound = 0
        else:

            model.reactions.get_by_id(k).lower_bound=0

        if k in p:
            model.reactions.get_by_id(k).upper_bound = fr[fr.metabolites == k].flux.iloc[0]*1.01
            model.reactions.get_by_id(k).lower_bound=0.0

        else:

            model.reactions.get_by_id(k).upper_bound = 0
    try:
        sol=model.optimize()
        print("Solution for alternativt {} is {}".format(alt,sol.objective_value))

    except:
        print("This alternative is not feasible {}".format(alt))

