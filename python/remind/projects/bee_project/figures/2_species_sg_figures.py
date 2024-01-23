import pandas as pd
import numpy as np
from figure_functions import *



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: ReMIND Team
"""
from remind.io.json import load_json_model
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse
from remind.optim.alternatives import find_alternative_solutions

from pytfa.optim.utils import symbol_sum
#
import os
import glob
from sys import argv
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

comm_name = '7_member_bee'
comm_name = '2_member_bee'
 #for 2speciespos int data #031022

model_of_int_list=[['Gilliamella_apicola_wkB1GF','Lactobacillus_mellifer_Bin4'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_mellifer_Bin4','Snodgrassella_alvi_wkB2GF'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_mellis_Hon2'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_mellis_Hon2','Snodgrassella_alvi_wkB2GF'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_apis_Hma11GF'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_apis_Hma11GF','Snodgrassella_alvi_wkB2GF'],
            ['Gilliamella_apicola_wkB1GF', 'Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'],
            ['Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'],
            ['Gilliamella_apicola_wkB1GF', 'Bifidobacterium_asteroides_PRL2011GF'],
            ['Gilliamella_apicola_wkB1GF', 'Snodgrassella_alvi_wkB2GF'],
            ]

model_of_int_list=[
            ['Gilliamella_apicola_wkB1GF', 'Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'],


            ]
def add_label(df, group=[""],name='abiotic', overshoot=0):
    df['label_{}'.format(name)] = df.groupby(by=group).grouper.group_info[0] + overshoot

output_file_figures="/remind/projects/bee_project/constrained_medium/analysis/figures_updated"
if not os.path.exists(output_file_figures):
    os.makedirs(output_file_figures)


data_path = '/remind/projects/bee_project/ilp_solutions_031022_{}'.format(comm_name)



data_path = '/remind/projects/bee_project/241022_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/101122_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/101122_2_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_2_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_3_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_new_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_newest_2_data_integ_ilp_solutions_2_member_bee'.format(comm_name)

data_path='/remind/projects/bee_project/060623_ilp_solutions_nplusone_2_member_bee/obj_num_2_modelofint_10_alternatives_interaction_positive_directional_general_yield/'

data_path="/remind/projects/bee_project/ilp_solutions_280823_nplusone_2_member_bee"



data_path="/remind/projects/bee_project/111023_2member_ilp_solutions_031022_2_member_bee"




data_path="/remind/projects/bee_project/111023_2member_ilp_solutions_031022_2_member_bee/"

allDirs = os.listdir(data_path)

obj_numbers = {1: 'min competition',
               2: 'max cooperation',
               3: 'max competition',
               4: 'min cooperation',
               5: "weighted",
               6: "min abiotic sources"}

obj_num_sep=True
list_ = []
for k in allDirs:
    # k.split('obj_num_')
    if obj_num_sep:
        obj_num = int(k.split('obj_num_')[1][0])
    if "alternatives" in k:
        alternation = k.split('alternatives_')[1]

    elif "alternation" in k:
        alternation = k.split('alternation_')[1]

    path_new = data_path + '/' + k
    allFiles = glob.glob(path_new + "/*" + "h5")
    for file in allFiles:
        if 'sols_active' in file:
            d = pd.read_hdf(file)
            if obj_num_sep:
                d['obj_num'] = [obj_num] * d.shape[0]
            d['alternation'] = [alternation] * d.shape[0]
            if "modelofint" in file:
                d['modelofint_no'] = [int(file.split('modelofint_')[1][0])] * d.shape[0]

                d['modelofint'] = tuple([model_of_int_list[int(file.split('modelofint_')[1][0])]]) * d.shape[0]
            list_.append(d)

# contains merged all groups
df = pd.concat(list_, ignore_index=True)

df['obj_num'] = [obj_numbers[k] for k in df.obj_num]
df.groupby(['obj_num', 'alternation']).objective.unique()
df.groupby(['obj_num', 'alternation']).alternative.nunique()
gb = df.groupby(['obj_num', 'alternation'])
keys = [k for k in gb.groups.keys()]
keys_int = [k for k in keys]
# keys_int=[k for k in keys if 'interaction' in k]

# newgrouping for n+1 studies
if "modelofint" in df.columns:
    gb = df.groupby(['modelofint'])
    keys = [k for k in gb.groups.keys()]
    keys_int = [k for k in keys]

    # groupofints = [keys_int[3], keys_int[4]]
    #
    # groupofints = [keys_int[5], keys_int[6]]
    # d = pd.concat([gb.get_group(gr) for gr in groupofints])
    grouping = ['alternation', 'modelofint','obj_num']
    gb = df.groupby(grouping)
    keys = [k for k in gb.groups.keys()]
    keys_int = [k for k in keys]

d_int=gb.get_group(keys_int[0])
d_int_direct=gb.get_group(keys_int[1])
d_int_yield=gb.get_group(keys_int[2])

models_dict={
             "model1":"Gapicola",
             "model2":"Salvi"}
frame_int=get_species_data(d_int,models_dict=models_dict,len_models=2)
frame_int=get_species_data(d_int_direct,models_dict=models_dict,len_models=2)
frame_int_yield=get_species_data(d_int_yield,models_dict=models_dict,len_models=2)
#in this version we had helper
#sort the positive interactions

frame_int['pos_int']=[tuple(sorted(list(k))) for k in frame_int.pos_int]
frame_int_yield['pos_int']=[tuple(sorted(list(k))) for k in frame_int_yield.pos_int]


#for data  integration objective function
# frame_int=get_species_data(df,models_dict=models_dict,len_models=2)
# frame_int_unique=frame_int.groupby('abiotic').first().reset_index()
# frame_int_unique['len_abiot']=[len(k)for k in frame_int_unique.abiotic.unique()]

# list_int=[]
# for k in keys_int:
#     f=gb_all.get_group(k)
#     list_int.append(f)

label=True
"only for positive_interaction"
add_label(frame_int,['pos_int'],name='interaction')
add_label(frame_int_yield,['pos_int'],name='interaction')

label=True
# positive
l = []
for k in frame_int.pos_int.unique():
    l += list(k)
# negative
l_ = []
for k in frame_int.neg_int.unique():
    l_ += list(k)

s_ = pd.DataFrame(l).value_counts()
s__ = pd.DataFrame(l_).value_counts()

cols = [k[0] for k in s_.index]
if label:
    cols+=['label_interaction']
dummyarray = np.empty((frame_int.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols)
for k in range(frame_int.alternative.nunique()):
    for met in cols:
        if met in frame_int.iloc[k].pos_int:
            if met in frame_int.Salvi_secretions.iloc[k]:
                if met in frame_int.Gapicola_uptake.iloc[k]:
                    data_coop.iloc[k][met] = 0.5
            elif met in frame_int.Gapicola_secretions.iloc[k]:
                if met in frame_int.Salvi_uptake.iloc[k]:
                    data_coop.iloc[k][met] = 0.0
    if 'label_interaction' in frame_int.columns:
        data_coop['label_interaction'].iloc[k] = frame_int.iloc[k].label_interaction

helper = data_coop.groupby('label_interaction').apply(lambda gr: gr.nunique())

"used for posiitve interaction and directional"
"usage of  helper"
cols.remove('label_interaction')

frame_int_unique = frame_int.groupby('label_interaction').first().reset_index()
dummyarray = np.empty((frame_int_unique.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols)
for k in range(frame_int_unique.alternative.nunique()):
    label_interaction = frame_int_unique.iloc[k].label_interaction
    for met in cols:
        # now helper is added
        if helper.iloc[label_interaction][met] <= 1.0:
            if met in frame_int_unique.iloc[k].pos_int:
                if met in frame_int_unique.Salvi_secretions.iloc[k]:
                    if met in frame_int_unique.Gapicola_uptake.iloc[k]:
                        data_coop.iloc[k][met] = 0.5
                elif met in frame_int_unique.Gapicola_secretions.iloc[k]:
                    if met in frame_int_unique.Salvi_uptake.iloc[k]:
                        data_coop.iloc[k][met] = 1.0
        else:
            #this indicates it can be both ways
            if met in frame_int_unique.iloc[k].pos_int:
                data_coop.iloc[k][met] = 0.0


cols_pos_int=[k[0] for k in s_.index]
data_coop['label_interaction']=frame_int_unique.label_interaction
size = data_coop.shape[0]
data_coop_ap = data_coop.append(data_coop.count(), ignore_index=True)
# data_coop_ap = data_coop.append(data_coop.where(data_coop.notna()).sum(), ignore_index=True)
data_coop_ap = data_coop_ap.sort_values(by=size, axis=1, ascending=False)
data_coop = data_coop_ap.drop(index=size, axis=1)
cols = [k for k in data_coop.columns]
data_coop['len_pos_int']=data_coop[cols_pos_int].count(axis=1)
data_coop['direct']=data_coop[cols_pos_int].nunique(axis=1)
data_coop['direct']=data_coop[cols_pos_int].nunique(axis=1)-data_coop[cols_pos_int].sum(axis=1)/3

data_coop_sorted=data_coop.sort_values(['len_pos_int','direct'],ascending=[True,False]) #before false
data_coop_sorted=data_coop_sorted.drop(['len_pos_int','direct'],axis=1)


##new
data_coop_sorted['direct_new']=data_coop_sorted[cols_pos_int].nunique(axis=1)
data_coop_sorted['direct_val_new']=data_coop_sorted[cols_pos_int].min(axis=1)



fr=pd.DataFrame()
gr='label_interaction'
fr['Gapic_yield_min']=frame_int_yield.groupby(gr).Gapicola_yield.min().values
fr['Gapic_yield_max']=frame_int_yield.groupby(gr).Gapicola_yield.max().values

fr['Salvi_yield_min']=frame_int_yield.groupby(gr).Salvi_yield.min().values
fr['Salvi_yield_max']=frame_int_yield.groupby(gr).Salvi_yield.max().values

#here complicated
# frame_int_yield_unique=frame_int_yield.groupby('label_interaction').first().reset_index()
#do reindexing
data_coop_sorted.set_index('label_interaction', inplace=True)
cols = [k for k in data_coop_sorted.columns]
fr=fr.reindex(data_coop_sorted.index)

#plot data_coop_sorted and fr



##for Data Integration
d_int_abiotic_yield=gb.get_group(keys_int[4])
frame_int=get_species_data(d_int_abiotic_yield,models_dict=models_dict,len_models=2)


frame_int_all = frame_int
add_label(frame_int_all, "abiotic")

# add label


label = True
f_ = []
for k in frame_int.abiotic.unique():
    f_ += list(k)

f__ = pd.DataFrame(f_).value_counts()

cols = [k[0] for k in f__.index]
if label:
    cols_all = cols + ['label_abiotic']
else:
    cols_all = cols

"""copying"""
dummyarray = np.empty((frame_int.alternative.nunique(), len(cols_all)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols_all)
for k in range(frame_int_all.alternative.nunique()):
    for met in cols:
        if (met in frame_int_all.Salvi_uptake.iloc[k]):
            if not (met in frame_int_all.pos_int.iloc[k]):
                data_coop.iloc[k][met] = 0.5
        elif met in frame_int_all.Gapicola_uptake.iloc[k]:
            if not (met in frame_int_all.pos_int.iloc[k]):
                data_coop.iloc[k][met] = 0.0
        if (met in frame_int_all.Gapicola_uptake.iloc[k]) and (met in frame_int_all.Salvi_uptake.iloc[k]):
            data_coop.iloc[k][met] = 1.0

    if 'label_abiotic' in cols_all:
        data_coop['label_abiotic'].iloc[k] = frame_int_all.iloc[k].label_abiotic


helper = data_coop.groupby('label_abiotic').apply(lambda gr: gr.nunique())
cols_all = cols
frame_int = frame_int_all.groupby('abiotic').first()#.reset_index()
dummyarray = np.empty((frame_int.alternative.nunique(), len(cols_all)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols_all)
for k in range(frame_int.alternative.nunique()):
    label_abiot = frame_int.iloc[k].label_abiotic
    for met in cols:
        # now helper is added
        if helper.iloc[label_abiot][met] <= 1.0:
            if (met in frame_int.Salvi_uptake.iloc[k]):
                if not (met in frame_int.pos_int.iloc[k]):
                    data_coop.iloc[k][met] = 0.5
            elif met in frame_int.Gapicola_uptake.iloc[k]:
                if not (met in frame_int.pos_int.iloc[k]):
                    data_coop.iloc[k][met] = 0.3
            if (met in frame_int.Gapicola_uptake.iloc[k]) and (met in frame_int.Salvi_uptake.iloc[k]):
                data_coop.iloc[k][met] = 1.0
        else:
            data_coop.iloc[k][met] = 0.0

size = data_coop.shape[0]
data_coop_ap = data_coop.append(data_coop.sum()*data_coop.count(),ignore_index=True)
data_coop_ap = data_coop_ap.sort_values(by=size, axis=1, ascending=False)
data_coop = data_coop_ap.drop(index=size, axis=1)
cols = [k for k in data_coop.columns]

data_coop['size']=data_coop.count(axis=1)
data_coop= data_coop.sort_values('size')
data_coop=data_coop.drop('size',axis=1)
data_coop=data_coop.reindex([3,4,5,6,7,8,0,1,2])
cols = [k for k in data_coop.columns]
cols_new =[k for k in data_coop.columns]
# cols_new=['thm_e',
#  'arg__L_e',
#  'ade_e',
#  'fol_e',
#  'cit_e',
#  'pydx_e',
#  'nmn_e',
#  'his__L_e',
#  'glc__D_e',
#  'acmana_e',
#  '26dap__M_e',
#  'akg_e',
# 'glu__L_e',
#  'uri_e',
#  'cmp_e',
# 'dha_e',
#       ]
cols_new=['thm_e',
 'arg__L_e',
 'ade_e',
 'fol_e',
 'pydx_e',
 'nmn_e',
 'his__L_e',
 'acmana_e',
 '26dap__M_e',
 'uri_e',
 'cmp_e',
 'glu__L_e',
    'akg_e',
 'sucr_e',
 'glc__D_e',
 'fru_e',
]


data_coop_new=data_coop[cols_new]
cols = [k for k in data_coop_new.columns]

fr=pd.DataFrame()
frame_int_yield=frame_int_all
gr='label_abiotic'
fr['Gapic_yield_min']=frame_int_yield.groupby(gr).Gapicola_yield.min()
fr['Gapic_yield_max']=frame_int_yield.groupby(gr).Gapicola_yield.max().values

fr['Salvi_yield_min']=frame_int_yield.groupby(gr).Salvi_yield.min().values
fr['Salvi_yield_max']=frame_int_yield.groupby(gr).Salvi_yield.max().values


fr=fr.reindex([3,4,5,6,7,8,0,1,2])
"HERE IT IS TO GENERATE THE INTERACTION MAPS"


label = True
mets_ = []
for k in interested.Salvi_dimes.unique():
    mets_ += list(k)
for k in interested.Gapicola_dimes.unique():
    mets_ += list(k)

mets_ = list(set(mets_))
no_mets = len(mets_)
pos = get_coordinates_in_circle(no_mets)
# create the dictionary after create it for all unique mets
positions = dict(zip(mets_, pos))

data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'

data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_111023_bee.h5'
frame_dime = pd.read_hdf(data_path)

plt.rcParams.update({
    'font.family': 'Avenir'})
# chosen min alternative
interested_label = 0#0 #before 4
interested_alternative = 1
interested = frame_int_all[frame_int_all.label_abiotic == interested_label]
intr = interested.iloc[interested_alternative]
frame_s = frame_dime[(frame_dime.model == 'Snodgrassella_alvi_wkB2GF') & (frame_dime.yield_perc == intr.Salvi_yield) & (
            frame_dime.alternative == intr.Salvi_alt)]
frame_g = frame_dime[
    (frame_dime.model == 'Gilliamella_apicola_wkB1GF') & (frame_dime.yield_perc == intr.Gapicola_yield) & (
                frame_dime.alternative == intr.Gapicola_alt)]
frame_fig = frame_s.append(frame_g)
fig = get_interaction_map(frame_fig, positions, intr)
plt.savefig(
    output_file_figures + '/051123_070823_min_abiot_LABEL_{}_ALTERNATIVE_{}interaction_1_map_graph_function_bee_all_above_all_yield.svg'.format(
        interested_label, interested_alternative))
plt.close()

"""from here on 2 species abiotic plot with yields"""

"""SUBPLOT TRIAL for 2 species copied from 7 species script
"""
all_exchanges = pd.read_csv("/remind/projects/bee_project/stats/all_exchanges_stats.csv")
all_exchanges = all_exchanges.groupby('mets').first().reset_index()
met_name = [all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols]
# now heatmap


heatmap_abiotic = "tab20b"
# heatmap_abiotic = "Accent"

# heatmap_abiotic = "Paired"
# heatmap_abiotic = "Set2"
heatmap_yield = "rainbow"
heatmap_yield = "coolwarm"

# df_pivoted = data_coop_sorted
df_pivoted = data_coop_new

# df_pivoted = data_distance_fornow
# fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(30, 20)) #20 was10 before
fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(50, 50)) #20 was10 before

# fig, ax = plt.subplots(figsize=(30, 10))
no_categ = 4
# no_categ = 5
heatmap = cm.get_cmap(heatmap_abiotic, lut=no_categ)
# heatmap = cm.Greens
alpha_chosen = 0.9
current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
cax2 = ax.imshow(df_pivoted, cmap=current_cmap, vmin=df_pivoted.min().min(), vmax=df_pivoted.max().max(),
                 alpha=alpha_chosen)
# for the yield only
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=0.0,vmax=1.0,alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax.xaxis.tick_top()
ax.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                          df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax.set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
                          df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
# ax.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=12,style='italic')
# ax.set_xticklabels([str(k) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=12,style='italic')
# use this normally
ax.set_xticklabels([all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols], rotation=90, fontsize=16,
                   style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.index.to_list()],rotation=0,fontsize=8,style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in yield_df.met_name.unique()],rotation=0,fontsize=8,style='italic')
ax.set_yticklabels([k + 1 for k in df_pivoted.index.to_list()], rotation=0, fontsize=16)

# for annotation
for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):
        text = (df_pivoted.iloc[y, x])
        if not math.isnan(text):
            text = str(np.round(float(text), 1))
            # plt.text(x , y ,  text,
            #      horizontalalignment='center',
            #      verticalalignment='center',
            #     fontsize=15 )

# ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.yaxis.set_minor_locator(AutoMinorLocator())

#
ax.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax)
# cbar.set_ticks(np.arange(-.5, 6.5, 1))
# depends on the value of the lables

# to put ticks and labels if required
ticks = [k + 0.25 / 2 for k in np.arange(0.0, 1.0, 0.25)]
cbar.set_ticks(ticks)
labels = [k for k in np.unique(df_pivoted) if not math.isnan(k)]
cbar.set_ticklabels(labels)
# end puttticks and labels
# cbar.set_ticklabels(np.linspace(1, 7, 7))
plt.tight_layout()
# plt.savefig(output_file_figures+'/bsg_min_abiot.svg')


# plt.savefig(output_file_figures+'/nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
# # plt.savefig(output_file_figures+'/SALVI_yield_nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
# plt.close('all')



# second figure
# frame_int=frame_int_unique
frame_int=fr
df_pivoted = frame_int[['Gapic_yield_min', 'Gapic_yield_max']]

no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.tab20b
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

alpha_chosen = 0.6
current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax2.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax2.xaxis.tick_top()
ax2.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax2.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax2.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
                    style='italic')
# ax.set_xticklabels([str(k) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=12,style='italic')
# use this normally
# ax.set_xticklabels([all_exchanges[all_exchanges.mets==k].name_met.values[0] for k in cols],rotation=90,fontsize=16,style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.index.to_list()],rotation=0,fontsize=8,style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in yield_df.met_name.unique()],rotation=0,fontsize=8,style='italic')
# ax2.set_yticklabels([k+1 for k in df_pivoted.index.to_list()],rotation=0,fontsize=16)

# for annotation
for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):
        text = (df_pivoted.iloc[y, x])
        if not math.isnan(text):
            text = str(np.round(float(text), 1))
            # ax2.text(x , y ,  text,
            #      horizontalalignment='center',
            #      verticalalignment='center',
            #     fontsize=15 )

#
ax2.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax2.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax2.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
cbar = fig.colorbar(cax2, ax=ax2)
# plt.tight_layout()
# plt.savefig(output_file_figures+'/subplot_nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
# # plt.savefig(output_file_figures+'/SALVI_yield_nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
# plt.close('all')


# third figure
df_pivoted = frame_int[['Salvi_yield_min', 'Salvi_yield_max']]
no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.summer
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax3.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax3.xaxis.tick_top()
ax3.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax3.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax3.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
                    style='italic')

# for annotation
for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):
        text = (df_pivoted.iloc[y, x])
        if not math.isnan(text):
            text = str(np.round(float(text), 1))
            # ax3.text(x , y ,  text,
            #      horizontalalignment='center',
            #      verticalalignment='center',
            #     fontsize=15 )

ax3.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax3.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax3.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endppwdoint=True)
cbar = fig.colorbar(cax2, ax=ax3)
plt.tight_layout()
plt.savefig(output_file_figures + '/051123_abiot_2321023_2_pos_int_070823_NEW_2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
fig.tight_layout()
plt.close('all')


""" from here on only to show interaction heatmap"""
fr_int=frame_int_all[frame_int_all.label_abiotic==10]
fr_int=frame_int_all[frame_int_all.label_abiotic==4]

label=True
# positive
l = []
for k in fr_int.pos_int.unique():
    l += list(k)
# negative
l_ = []
for k in fr_int.neg_int.unique():
    l_ += list(k)

cols=set(l+l_)
# s_ = pd.DataFrame(l).value_counts()
# s__ = pd.DataFrame(l_).value_counts()
cols=[k for k in cols]
# cols = [k[0] for k in s_.index]

frame_int_unique=fr_int
dummyarray = np.empty((frame_int_unique.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols)
for k in range(frame_int_unique.alternative.nunique()):
    for met in cols:
        if met in frame_int_unique.iloc[k].pos_int:
            if met in frame_int_unique.Salvi_secretions.iloc[k]:
                if met in frame_int_unique.Gapicola_uptake.iloc[k]:
                    data_coop.iloc[k][met] = 0.5
            elif met in frame_int_unique.Gapicola_secretions.iloc[k]:
                if met in frame_int_unique.Salvi_uptake.iloc[k]:
                    data_coop.iloc[k][met] = 0.0
        elif met in frame_int_unique.iloc[k].neg_int:
            if (met in frame_int_unique.Gapicola_uptake.iloc[k]) and (met in frame_int_unique.Salvi_uptake.iloc[k]):
                data_coop.iloc[k][met] = 1.0

        else:


            if met in frame_int_unique.Salvi_uptake.iloc[k]:
                 if met not in frame_int_unique.iloc[k].pos_int:
                     data_coop.iloc[k][met] = 0.2

            if met in frame_int_unique.Gapicola_uptake.iloc[k]:
                if met not in frame_int_unique.iloc[k].pos_int:
                    data_coop.iloc[k][met] = 0.6

cols_new=['mal__L_e',
          'akg_e',
 'dha_e',
 'lac__D_e',
          'glu__L_e',
    'ade_e',
 'thm_e',
'arg__L_e',

]


cols_new=[
    'ade_e',
 'thm_e',
'arg__L_e',
'glu__L_e',
'akg_e',
'lac__D_e'
]
data_coop=data_coop[cols_new]
cols=cols_new

fr=pd.DataFrame()
fr['Gapicola_yield']=frame_int_unique.Gapicola_yield
fr['Salvi_yield']=frame_int_unique.Salvi_yield



#here plot this
heatmap_abiotic='Set3'
df_pivoted=data_coop
# df_pivoted = data_distance_fornow
fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(30, 20)) #20 was10 before

# fig, ax = plt.subplots(figsize=(30, 10))
no_categ = 5
# no_categ = 5
heatmap = cm.get_cmap(heatmap_abiotic, lut=no_categ)
# heatmap = cm.Greens
alpha_chosen = 1.0
current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
cax2 = ax.imshow(df_pivoted, cmap=current_cmap, vmin=df_pivoted.min().min(), vmax=df_pivoted.max().max(),
                 alpha=alpha_chosen)
# for the yield only
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=0.0,vmax=1.0,alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax.xaxis.tick_top()
ax.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                          df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax.set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
                          df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
# ax.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=12,style='italic')
# ax.set_xticklabels([str(k) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=12,style='italic')
# use this normally
ax.set_xticklabels([all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols], rotation=90, fontsize=16,
                   style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.index.to_list()],rotation=0,fontsize=8,style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in yield_df.met_name.unique()],rotation=0,fontsize=8,style='italic')
ax.set_yticklabels([k + 1 for k in df_pivoted.index.to_list()], rotation=0, fontsize=16)

# for annotation
for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):
        text = (df_pivoted.iloc[y, x])
        if not math.isnan(text):
            text = str(np.round(float(text), 1))
            # plt.text(x , y ,  text,
            #      horizontalalignment='center',
            #      verticalalignment='center',
            #     fontsize=15 )

# ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.yaxis.set_minor_locator(AutoMinorLocator())

#
ax.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax)
# cbar.set_ticks(np.arange(-.5, 6.5, 1))
# depends on the value of the lables

# to put ticks and labels if required
ticks = [k + 0.25 / 2 for k in np.arange(0.0, 1.0, 0.25)]
cbar.set_ticks(ticks)
labels = [k for k in np.unique(df_pivoted) if not math.isnan(k)]
cbar.set_ticklabels(labels)
# end puttticks and labels
# cbar.set_ticklabels(np.linspace(1, 7, 7))
plt.tight_layout()
# plt.savefig(output_file_figures+'/bsg_min_abiot.svg')


# plt.savefig(output_file_figures+'/nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
# # plt.savefig(output_file_figures+'/SALVI_yield_nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
# plt.close('all')



# second figure
# frame_int=frame_int_unique
frame_int=fr
df_pivoted = frame_int[['Gapicola_yield']]

no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.tab20b
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

alpha_chosen = 0.6
current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax2.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax2.xaxis.tick_top()
ax2.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax2.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax2.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
                    style='italic')
# ax.set_xticklabels([str(k) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=12,style='italic')
# use this normally
# ax.set_xticklabels([all_exchanges[all_exchanges.mets==k].name_met.values[0] for k in cols],rotation=90,fontsize=16,style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.index.to_list()],rotation=0,fontsize=8,style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in yield_df.met_name.unique()],rotation=0,fontsize=8,style='italic')
# ax2.set_yticklabels([k+1 for k in df_pivoted.index.to_list()],rotation=0,fontsize=16)

# for annotation
for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):
        text = (df_pivoted.iloc[y, x])
        if not math.isnan(text):
            text = str(np.round(float(text), 1))
            # ax2.text(x , y ,  text,
            #      horizontalalignment='center',
            #      verticalalignment='center',
            #     fontsize=15 )

#
ax2.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax2.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax2.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
cbar = fig.colorbar(cax2, ax=ax2)
# plt.tight_layout()
# plt.savefig(output_file_figures+'/subplot_nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
# # plt.savefig(output_file_figures+'/SALVI_yield_nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
# plt.close('all')


# third figure
df_pivoted = frame_int[['Salvi_yield']]
no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.summer
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax3.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax3.xaxis.tick_top()
ax3.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax3.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax3.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
                    style='italic')

# for annotation
for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):
        text = (df_pivoted.iloc[y, x])
        if not math.isnan(text):
            text = str(np.round(float(text), 1))
            # ax3.text(x , y ,  text,
            #      horizontalalignment='center',
            #      verticalalignment='center',
            #     fontsize=15 )

ax3.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax3.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax3.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax3)
plt.tight_layout()
plt.savefig(output_file_figures + '/Interested_interaction_070823_NEW_2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
fig.tight_layout()
plt.close('all')

