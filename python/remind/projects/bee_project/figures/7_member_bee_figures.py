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
from remind.utils.postprocessing import *
comm_name = '7_member_bee'


output_file_figures="/remind/projects/bee_project/constrained_medium/analysis/figures_trial"

output_file_figures="/remind/projects/bee_project/constrained_medium/analysis/figures_pray"

if not os.path.exists(output_file_figures):
    os.makedirs(output_file_figures)


# model = load_json_model('../projects/community_models/{}_community.json'.format(comm_name))

alternation_opts = ['reaction', 'interaction']


data_path='/remind/projects/bee_project/ilp_solutions_051022_correct_allpos_7_member_bee'
data_path='/remind/projects/bee_project/ilp_solutions_051022_correct_7_member_bee'


data_path='/remind/projects/bee_project/ilp_solutions_290623_7_member_bee'

data_path='/remind/projects/bee_project/121023_ilp_solutions_290623_7_member_bee'

# data_path="/remind/projects/bee_project/ilp_solutions_290623_7_member_bee/
allDirs = os.listdir(data_path)

obj_numbers = {1: 'min competition',
               2: 'max cooperation',
               3: 'max competition',
               4: 'min cooperation',
               5: "weighted",
               6: "min abiotic sources"}


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

list_ = []
for k in allDirs:
    k.split('obj_num_')
    obj_num = int(k.split('obj_num_')[1][0])
    alternation = k.split('alternatives_')[1]
    path_new = data_path + '/' + k
    allFiles = glob.glob(path_new + "/*" + "h5")
    for file in allFiles:
        if 'sols_active' in file:
            d = pd.read_hdf(file)
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

# d = gb.get_group(keys_int[1])
d=gb.get_group(keys_int[1])
d2=gb.get_group(keys_int[2])

from figure_functions import *
models_dict = {"model1": "Bifido",
               "model2": "Gapicola",
               "model3": "Lapis",
               "model4": "Lkulla",
               "model5": "Lmellifer",
               "model6": "Lmellis",
               "model7": "Salvi"}

frame_int=get_species_data(d,models_dict=models_dict,len_models=7)
frame_abiot=get_species_data(d2,models_dict=models_dict,len_models=7)

#to get the length of abiotics
abiotic_len=[]
for kk in range(frame_int.shape[0]):
    abiotic_len+=[len([k for k in frame_int.upt_counts.iloc[kk].keys() if k not in frame_int.pos_int.iloc[kk]])]

frame_int['len_abiot']=abiotic_len
frame_abiot['len_pos_int']=[len(k) for k in frame_abiot.pos_int]
frame_abiot['len_neg_int']=[len(k) for k in frame_abiot.neg_int]


frame_int_unique=frame_int.groupby('pos_int').first().reset_index()
for key, value in models_dict.items():
    frame_int_unique['{}_yield_min'.format(value)]=frame_int.groupby('pos_int')['{}_yield'.format(value)].min().values
    frame_int_unique['{}_yield_max'.format(value)]=frame_int.groupby('pos_int')['{}_yield'.format(value)].max().values

""" to create data """
# list_upt = []
# list_store = []
# for k in range(frame_int.shape[0]):
#     dat = frame_int[col_of_int].iloc[k]
#     dat_pos_int = frame_int['pos_int'].iloc[k]
#     print(k)
#     list_new = []
#     for d in dat:
#         list_upt += list(d)
#         list_new += list(d)
#     list_new = list(set(list_new))
#     list_new = [k for k in list_new if k not in list(dat_pos_int)]
#     # list_upt+=list(dat[0])
#     list_store.append(list_new)
# frame_int['uptakes'] = list_store
# frame_int['len_abiot'] = [len(k) for k in frame_int.uptakes]
# s = pd.DataFrame(list_upt).value_counts()


# add_label(frame_int,['pos_int','neg_int'],name='interaction')

"only for positive_interaction"
add_label(frame_int,['pos_int'],name='interaction')

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
cols_neg = [k[0] for k in s__.index]
# cols += cols_neg
if label:
    cols+=['label_interaction']

"""new 100823"""

frame_int_all=frame_int

col_of_int="upt_counts"
col_of_int="secretion_counts"
cols_all = cols
dummyarray = np.empty((frame_int_all.alternative.nunique(), len(cols_all)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols_all)
for k in range(frame_int_all.alternative.nunique()):
    # label_abiot = frame_int.iloc[k].label_abiotic
    for met in cols:
        if 'label' not in met:
            if met in frame_int_all.iloc[k].pos_int:
                data_coop.iloc[k][met] = frame_int_all.iloc[k][col_of_int][met]
        else:
            data_coop.iloc[k][met] = frame_int_all.iloc[k]['label_interaction']



dummyarray2 = np.empty((data_coop.label_interaction.nunique(), len(cols_all)))
dummyarray2[:] = np.nan
data_coop_abiot_7 = pd.DataFrame(dummyarray2, columns=cols_all)
cols_count=[]
for label in data_coop.label_interaction.unique():
    d=data_coop[data_coop.label_interaction==label]
    for col in d:
        col_count=len(d[col].unique())
        cols_count+=([col_count])
        label_interaction = d['label_interaction'].unique()[0]
        print(col_count,d[col].unique())
        data_coop_abiot_7.iloc[int(label_interaction)]['label_interaction'] = label_interaction
        if col != 'label_interaction':
            data_coop_abiot_7.iloc[int(label_interaction)][col] = max(d[col].unique())
            # data_coop_abiot_3.iloc[k]['met']=d[col].unique()[0]
        # print(len(d[col].unique()))
    # print(col,d[col].unique())

np.unique(data_coop_abiot_7[cols_all].values)


data_coop_uptakes=data_coop_abiot_7
data_coop_secretions=data_coop_abiot_7

"plot"

cols_new=[
    'mal__L_e',
    'lac__D_e',
    'gly_e',
    'glcn_e',
    'glyald_e',
    'fum_e',
 'ura_e',
    'akg_e',
 'trp__L_e',
 'succ_e',
 'pro__L_e',
 'co2_e',
 '4hthr_e',
 'dha_e',
    'glyc_e',
 'acald_e',
 'for_e',
 'xan_e',
 'pnto__R_e',
 'asp__L_e',
 'xtsn_e']

cols_new=[k for k in cols if "label_" not in k]
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(30, 30))
# fig, ax = plt.subplots(figsize=(30, 10))
no_categ = 5
no_categ = 7
# df_pivoted=(data_coop_uptakes+data_coop_secretions)[cols_new] #to add them
df_pivoted=data_coop_secretions[cols_new]
# df_pivoted=data_coop_uptakes[cols_new]

cols=[k for k in df_pivoted.columns]
heatmap = cm.get_cmap("tab20b", lut=no_categ)
# heatmap = cm.get_cmap("Accent", lut=no_categ)

# heatmap = cm.Greens
alpha_chosen = 0.7
current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2 = ax.imshow(df_pivoted, cmap=current_cmap, vmin=df_pivoted.min().min(), vmax=df_pivoted.max().max(),
#                  alpha=alpha_chosen)

cax2 = ax.imshow(df_pivoted, cmap=current_cmap, vmin=1, vmax=7,
                 alpha=alpha_chosen)
# for the yield only
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=0.0,vmax=1.0,alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax.xaxis.tick_top()
ax.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                          df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax.set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
                          df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])

ax.set_xticklabels([all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols], rotation=90, fontsize=16,
                   style='italic')

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


ax.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax)
# ticks = [k + 0.25 / 2 for k in np.arange(0.0, 1.0, 0.25)]
ticks= [k + 0.25 / 2 for k in np.arange(0.0, 1.0, 0.25/2)]
# ticks= [k+0.1 for k in np.linspace(0,1,8)]
# cbar.set_ticks(ticks)
# labels = [k for k in np.unique(df_pivoted) if not math.isnan(k)]
# cbar.set_ticklabels(labels)
# end puttticks and labels
# cbar.set_ticklabels(np.linspace(1, 7, 7))
plt.tight_layout()
# plt.savefig(output_file_figures+'/081123_data_uptakes_secretions_coop_7_member.svg')

plt.savefig(output_file_figures+'/081123_data_secretions_coop_7_member.svg')

# plt.savefig(output_file_figures+'/081123_data_uptakes_coop_7_member.svg')
df_ess=df_pivoted.count()[df_pivoted.count()>=df_pivoted.count().max()]
essential_mets=[k for k in df_ess.index]

"put pie chard "
biggclasses=pd.read_csv('/remind/projects/bee_project/stats/bigg_classes_from_machado.tsv', sep="\t")
s=biggclasses[biggclasses.bigg.isin([k.replace('_e', '') for k in cols])].groupby('bigg').first().\
    reset_index()[['bigg', 'sub_class', 'name','class','super_class']]['sub_class'].value_counts()
plt.figure()
heatmap = cm.get_cmap("Set3", lut=len(cols))
s.plot(kind='pie',colormap=heatmap,fontsize=12,autopct='%1.0f%%')
plt.tight_layout()
plt.savefig(output_file_figures+'/coop_mets_081123_pie_7member_abiotic.svg')

s=biggclasses[biggclasses.bigg.isin([k.replace('_e', '') for k in essential_mets])].groupby('bigg').first().\
    reset_index()[['bigg', 'sub_class', 'name','class','super_class']]['sub_class'].value_counts()
plt.figure()
heatmap = cm.get_cmap("tab20", lut=len(essential_mets))
s.plot(kind='pie',colormap=heatmap,fontsize=12,autopct='%1.0f%%')
plt.tight_layout()
plt.savefig(output_file_figures+'/081123_pie_7member_abiotic_essential.svg')

alternating=[k for k in cols if k not in essential_mets]
s=biggclasses[biggclasses.bigg.isin([k.replace('_e', '') for k in alternating])].groupby('bigg').first().\
    reset_index()[['bigg', 'sub_class', 'name','class','super_class']]['sub_class'].value_counts()
plt.figure()
heatmap = cm.get_cmap("tab20", lut=len(alternating))
s.plot(kind='pie',colormap=heatmap,fontsize=12,autopct='%1.0f%%')
plt.tight_layout()
plt.savefig(output_file_figures+'/081123_pie_7member_abiotic_alternating.svg')

"""Creating the dataframe for primary and supplementary c sources """
# list possible and met_pool taken from dimes script
met_pool=[k for k in frame_union1.metabolites.unique()]
met_pool+=['EX_uaagmda_e' ,'EX_adn_e']
df_mets=pd.DataFrame()
list_possible=[k for k in frame_possible.annotation_bigg if "EX_" in k]
df_mets['reactions']=list_possible +met_pool
df_mets['type']=["primary"]*len(list_possible)+["supplementary"]*len(met_pool)
biggclasses=pd.read_csv('/remind/projects/bee_project/stats/bigg_classes_from_machado.tsv', sep="\t")
class_name=[]
sub_class_name=[]
met_list=[]
for rxn in df_mets.reactions.values:
    met=rxn.replace('EX_','')
    met_list+=[met]
    met=met.replace('_e','')
    if met in biggclasses.bigg.unique():
        class_name+=[biggclasses[biggclasses.bigg==met]['class'].iloc[0]]
        sub_class_name+=[biggclasses[biggclasses.bigg==met]['sub_class'].iloc[0]]

    else:
        class_name +=["none"]
        sub_class_name +=["none"]

df_mets['metabolite']=met_list
df_mets['metabolite_name']=[all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in met_list]
df_mets['sub_class']=sub_class_name
df_mets['class']=class_name
df_mets.to_csv('081123_extracellular_metabolite_list_for_honeybeegut.csv')

#cols being metabolites in the abiotic environment


"""end 100823"""
dummyarray = np.empty((frame_int.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols)
for k in range(frame_int.alternative.nunique()):
    for met in cols:
        if met in frame_int.iloc[k].pos_int:
            if met in frame_int.Salvi_secretions.iloc[k]:
                if met in frame_int.Gapic_uptake.iloc[k]:
                    data_coop.iloc[k][met] = 0.5
            elif met in frame_int.Gapic_secretions.iloc[k]:
                if met in frame_int.Salvi_uptake.iloc[k]:
                    data_coop.iloc[k][met] = 0.0
        # else:
        #     if (met in frame_int.Gapic_uptake.iloc[k]) and (met in frame_int.Salvi_uptake.iloc[k]):
        #         data_coop.iloc[k][met] = 1.0

    if 'label_interaction' in frame_int.columns:
        data_coop['label_interaction'].iloc[k] = frame_int.iloc[k].label_interaction

helper = data_coop.groupby('label_interaction').apply(lambda gr: gr.nunique())
# helper['label_interaction']=data_coop.groupby('label_interaction').first()

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

size = data_coop.shape[0]
data_coop_ap = data_coop.append(data_coop.where(data_coop.notna()).sum(), ignore_index=True)
data_coop_ap = data_coop_ap.sort_values(by=size, axis=1, ascending=False)
data_coop = data_coop_ap.drop(index=size, axis=1)
cols = [k for k in data_coop.columns]

data_coop_label=data_coop.copy()
data_coop_label['label_interaction']=frame_int_unique.label_interaction
data_coop_label['size']=data_coop.count(axis=1)
data_coop_label=data_coop_label.sort_values('size',ascending=False)
data_coop_label=data_coop_label.drop('size',axis=1)
#frame_int_yield is the yield directional positive interaction dataframe
fr=pd.DataFrame()
gr='label_interaction'
fr['Gapic_yield_min']=frame_int_yield.groupby(gr).Gapicola_yield.min()
fr['Gapic_yield_max']=frame_int_yield.groupby(gr).Gapicola_yield.max().values

fr['Salvi_yield_min']=frame_int_yield.groupby(gr).Salvi_yield.min().values
fr['Salvi_yield_max']=frame_int_yield.groupby(gr).Salvi_yield.max().values


#here complicated
frame_int_yield_unique=frame_int_yield.groupby('label_interaction').first().reset_index()

#do reindexing
data_coop_label.set_index('label_interaction', inplace=True)
frame_int_yield_unique.set_index('label_interaction', inplace=True)
frame_int_yield_unique=frame_int_yield_unique.reindex(data_coop_label.index)


fr=fr.reindex(data_coop_label.index)

"real coop"

frame_int_unique=frame_int.groupby('label_interaction').first().reset_index()
dummyarray = np.empty((frame_int_unique.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols)
for k in range(frame_int_unique.alternative.nunique()):
    for met in cols:
        if met in frame_int_unique.iloc[k].pos_int:
            if met in frame_int_unique.Salvi_secretions.iloc[k]:
                if met in frame_int_unique.Gapic_uptake.iloc[k]:
                    data_coop.iloc[k][met] = 0.5
            elif met in frame_int_unique.Gapic_secretions.iloc[k]:
                if met in frame_int_unique.Salvi_uptake.iloc[k]:
                    data_coop.iloc[k][met] = 0.0
        elif met in frame_int_unique.iloc[k].neg_int:
            if (met in frame_int_unique.Gapic_uptake.iloc[k]) and (met in frame_int.Salvi_uptake.iloc[k]):
                data_coop.iloc[k][met] = 1.0

    if 'label_interaction' in frame_int.columns:
        data_coop['label_interaction'].iloc[k] = frame_int_unique.iloc[k].label_interaction

###
data_coop=data_coop.drop('label_interaction',axis=1)

frame_int_unique['Gapic_yield_min']=frame_int.groupby('label_interaction').Gapicola_yield.min().values
frame_int_unique['Gapic_yield_max']=frame_int.groupby('label_interaction').Gapicola_yield.max().values

frame_int_unique['Salvi_yield_min']=frame_int.groupby('label_interaction').Salvi_yield.min().values
frame_int_unique['Salvi_yield_max']=frame_int.groupby('label_interaction').Salvi_yield.max().values

cols=[k for k in data_coop.columns]

"here a helper dataframe is needed "
all_exchanges = pd.read_csv("/remind/projects/bee_project/stats/all_exchanges_stats.csv")
all_exchanges = all_exchanges.groupby('mets').first().reset_index()
met_name = [all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols]
# now heatmap
df_pivoted = data_coop


df_pivoted=data_distance
df_pivoted=data_coop
# figsize=(5,10)
# fig, ax = plt.subplots(figsize=(20, 10))
fig, ax = plt.subplots(figsize=(20, 50))
# heatmap=cm.Paired
heatmap=cm.Accent
# heatmap=cm.Dark
# get how many colors u want
# number of categories
no_categ = 7
# heatmap = cm.get_cmap("tab10", lut = no_categ)
heatmap = cm.get_cmap("tab20b", lut=no_categ)
# heatmap = cm.Greens
alpha_chosen = 0.9
# cax1=ax.imshow(df_pivoted,cmap=heatmap,vmin=0,vmax=1,alpha=0.5)
# cax1=ax.imshow(df_pivoted,cmap=heatmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# ax.imshow(df_pivoted,cmap=current_cmap,vmin=0,vmax=1,alpha=1.0)
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
plt.savefig(output_file_figures + '/081123_131023_100823_7_species_abiotic_new.svg')

# plt.savefig(output_file_figures + '/7_species_pos_abiotic.svg')

# plt.savefig(output_file_figures + '/7_species_nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs, interested_col))
# plt.savefig(output_file_figures+'/SALVI_yield_nn_7230922_2_member_min_abiot_{}__{}.svg'.format(subs,interested_col))
plt.close('all')

"""for abiotic data 2 species"""


def add_label(df, group=[""],name='abiotic', overshoot=0):
    df['label_{}'.format(name)] = df.groupby(by=group).grouper.group_info[0] + overshoot


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
        elif met in frame_int_all.Gapic_uptake.iloc[k]:
            if not (met in frame_int_all.pos_int.iloc[k]):
                data_coop.iloc[k][met] = 0.0
        if (met in frame_int_all.Gapic_uptake.iloc[k]) and (met in frame_int_all.Salvi_uptake.iloc[k]):
            data_coop.iloc[k][met] = 1.0

    if 'label_abiotic' in cols_all:
        data_coop['label_abiotic'].iloc[k] = frame_int_all.iloc[k].label_abiotic


helper = data_coop.groupby('label_abiotic').apply(lambda gr: gr.nunique())
cols_all = cols
frame_int = frame_int_all.groupby('abiotic').first().reset_index()
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
            elif met in frame_int.Gapic_uptake.iloc[k]:
                if not (met in frame_int.pos_int.iloc[k]):
                    data_coop.iloc[k][met] = 0.3
            if (met in frame_int.Gapic_uptake.iloc[k]) and (met in frame_int.Salvi_uptake.iloc[k]):
                data_coop.iloc[k][met] = 1.0
        else:
            data_coop.iloc[k][met] = 0.0

size = data_coop.shape[0]
data_coop_ap = data_coop.append(data_coop.sum(), ignore_index=True)
data_coop_ap = data_coop_ap.sort_values(by=size, axis=1, ascending=False)
data_coop = data_coop_ap.drop(index=size, axis=1)
cols = [k for k in data_coop.columns]

data_coop_label=data_coop.copy()
data_coop_label['label_abiotic']=frame_int_all.groupby('abiotic').first().label_abiotic.values


data_coop_label=data_coop_label.iloc[index_new].reset_index().drop('index',axis=1)
data_coop=data_coop_label.drop('label_abiotic',axis=1)

"reindex"
data_coop_label.set_index('label_abiotic', inplace=True)
frame_int.set_index('label_abiotic', inplace=True)

frame_int=frame_int.reindex(data_coop_label.index)


Z = hierarchy.linkage(data_distance, method='average',metric='cityblock')
# metric='cityblock')
plt.figure()
plt.title("Dendrograms")
# Dendrogram plotting using linkage matrix
dendrogram = hierarchy.dendrogram(Z,labels=[k+1 for k in data_distance.index])
plt.savefig(output_file_figures+'/trial_dendrogram_7.svg')
index_new=dendrogram['leaves']

data_distance_fornow=data_distance.iloc[index_new]
data_distance_fornow=data_distance_fornow.reindex(columns=cols)


"try also for three-members "
"this is just to put everything to 1 for the distance matrix"
frame_int_unique=frame_int_all.groupby('abiotic').first().reset_index()
dummyarray = np.empty((frame_int_unique.alternative.nunique(), len(cols_all)))
dummyarray[:] = np.nan
data_distance = pd.DataFrame(dummyarray, columns=cols_all)
for k in range(frame_int_unique.alternative.nunique()):
    for met in cols:
        if (met in frame_int_unique.abiotic.iloc[k]):
            data_distance.iloc[k][met] = int(1)
        #
        else:
            data_distance.iloc[k][met] = int(0)


size = data_distance.shape[0]
data_distance_ap = data_distance.append(data_distance.where(data_distance.notna()).sum(), ignore_index=True)
data_distance_ap = data_distance_ap.sort_values(by=size, axis=1, ascending=False)
data_distance = data_distance_ap.drop(index=size, axis=1)
cols = [k for k in data_distance.columns]



#data_pos_int
"try also for three-members "
"this is just to put everything to 1 for the distance matrix"
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
cols_all=[k[0] for k in s_.index]
cols=cols_all
frame_int_unique=frame_int_all.groupby('pos_int').first().reset_index()
dummyarray = np.empty((frame_int_unique.alternative.nunique(), len(cols_all)))
dummyarray[:] = np.nan
data_distance = pd.DataFrame(dummyarray, columns=cols_all)
for k in range(frame_int_unique.alternative.nunique()):
    for met in cols:
        if (met in frame_int_unique.pos_int.iloc[k]):
            data_distance.iloc[k][met] = int(1)
        #
        # else:
        #     data_distance.iloc[k][met] = int(0)


""" abiotic data for 7 species"""
f_ = []
for k in frame_int.abiotic.unique():
    f_ += list(k)

f__ = pd.DataFrame(f_).value_counts()

cols = [k[0] for k in f__.index]
dummyarray = np.empty((frame_int.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols)
for k in range(frame_int.alternative.nunique()):
    for met in cols:
        if (met in frame_int.abiotic.iloc[k]):
            if (met in frame_int.neg_int.iloc[k]):
                data_coop.iloc[k][met] = 1.0

        if (met in frame_int.abiotic.iloc[k]):
            if not (met in frame_int.neg_int.iloc[k]):
                data_coop.iloc[k][met] = 0.5

# one way
data_coop_ap = data_coop.append(data_coop.sum(), ignore_index=True)
data_coop_ap = data_coop_ap.sort_values(by=24, axis=1, ascending=False)
data_coop = data_coop_ap.drop(index=24, axis=1)
cols = [k for k in data_coop.columns]


data_coop=data_coop.iloc[index_new].reset_index().drop('index',axis=1)

data_coop_label=data_coop_label.iloc[index_new].reset_index().drop('index',axis=1)
data_coop=data_coop_label.drop('label_abiotic',axis=1)

# second way

data_coop_ap = data_coop.append((data_coop < 0.6).where(data_coop.notna()).all(), ignore_index=True)
data_coop = data_coop_ap.sort_values(by=24, axis=1, ascending=False)

# third way
data_coop_ap = data_coop.append((data_coop > 0.6).where(data_coop.notna()).sum(), ignore_index=True)
data_coop_ap = data_coop_ap.sort_values(by=24, axis=1, ascending=False)
data_coop = data_coop_ap.drop(index=24, axis=1)
cols = [k for k in data_coop.columns]




"abiotic 7 species how many do uptake"
""" abiotic data for 7 species"""
"here coloring based on how many uptakes it"
f_ = []
for k in frame_int.abiotic.unique():
    f_ += list(k)

f__ = pd.DataFrame(f_).value_counts()

cols = [k[0] for k in f__.index]
dummyarray = np.empty((frame_int.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols)
for k in range(frame_int.alternative.nunique()):
    for met in cols:
        if (met in frame_int.abiotic.iloc[k]):
            if (met in frame_int.neg_int.iloc[k]):
                # print(frame_int.iloc[k].upt_counts[met])
                data_coop.iloc[k][met] = frame_int.iloc[k].upt_counts[met]

        if (met in frame_int.abiotic.iloc[k]):
            if not (met in frame_int.neg_int.iloc[k]):
                data_coop.iloc[k][met] = frame_int.iloc[k].upt_counts[met]

# multiply with number of not nan

data_coop_ap = data_coop.append((data_coop.sum()) * (1000 - data_coop.isna().sum() * 1000), ignore_index=True)

data_coop_ap = data_coop.append((data_coop.sum()) * (data_coop.count()+data_coop.notna().all()*2000), ignore_index=True)

# wrtlast_column
id = data_coop_ap.index.max()
data_coop_ap = data_coop_ap.sort_values(by=id, axis=1, ascending=False)
data_coop = data_coop_ap.drop(index=id, axis=1)
cols = [k for k in data_coop.columns]

"yield plot"
frame_int_yield = frame_int[['Gapicola_yield', 'Salvi_yield']]
df_pivoted = frame_int_yield
fig, ax = plt.subplots()
heatmap = cm.Reds
cax1 = ax.imshow(df_pivoted, cmap=heatmap, vmin=0, vmax=1, alpha=0.5)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
ax.imshow(df_pivoted, cmap=current_cmap, vmin=0, vmax=1, alpha=0.8)
ax.xaxis.tick_top()

ax.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                          df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax.set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
                          df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
                   style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.index.to_list()],rotation=0,fontsize=8,style='italic')
# ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in yield_df.met_name.unique()],rotation=0,fontsize=8,style='italic')
ax.set_yticklabels([k + 1 for k in df_pivoted.index.to_list()], rotation=0, fontsize=10)

#
for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):
        text = (df_pivoted.iloc[y, x])
        if not math.isnan(text):
            text = str(int(100 * np.round(float(text), 1))) + '%'
            plt.text(x, y, text,
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=6)

ax.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax.grid(which="minor", color="black", linestyle='-', linewidth=0.8)
# cbar = fig.colorbar(cax1)#, ticks=[0, 0.5, 1])

plt.tight_layout()

plt.savefig(
    output_file_figures + '/YIELD_MIN_ABIOTIC_COOPERATIONS_table2 memberbee_{}__{}.svg'.format(subs, interested_col))

"For competition data"
cols = [k[0] for k in s_.index]
cols_neg = [k[0] for k in s__.index]
cols = cols_neg
dummyarray = np.empty((frame_int.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols)
for k in range(frame_int.alternative.nunique()):
    for met in cols:
        if met in frame_int.Salvi_secretions.iloc[k]:
            if met in frame_int.Gapic_uptake.iloc[k]:
                data_coop.iloc[k][met] = 0.5
        elif met in frame_int.Gapic_secretions.iloc[k]:
            if met in frame_int.Salvi_uptake.iloc[k]:
                data_coop.iloc[k][met] = 0.0
        if (met in frame_int.Gapic_uptake.iloc[k]) and (met in frame_int.Salvi_uptake.iloc[k]):
            data_coop.iloc[k][met] = 1.0

''''''''
#TODO2SPECIES SUBPLOT
"""SUBPLOT TRIAL for 2 species"""
all_exchanges = pd.read_csv("/remind/projects/bee_project/stats/all_exchanges_stats.csv")
all_exchanges = all_exchanges.groupby('mets').first().reset_index()
met_name = [all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols]
# now heatmap
heatmap_abiotic = "tab20b"
heatmap_abiotic = "Accent"
heatmap_abiotic = "Paired"
heatmap_yield = "rainbow"

df_pivoted = data_coop
# df_pivoted = data_distance_fornow
fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(30, 20)) #20 was10 before

# fig, ax = plt.subplots(figsize=(30, 10))
no_categ = 5
# no_categ = 5
heatmap = cm.get_cmap(heatmap_abiotic, lut=no_categ)
# heatmap = cm.Greens
alpha_chosen = 0.7
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
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax3)
plt.tight_layout()
plt.savefig(output_file_figures + '/2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
fig.tight_layout()
plt.close('all')







"""get the interaction map"""


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
frame_dime = pd.read_hdf(data_path)

plt.rcParams.update({
    'font.family': 'Avenir'})
# chosen min alternative
interested_label = 2 #before 0 and 4
interested_alternative = 0
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
    output_file_figures + '/LABEL_{}_ALTERNATIVE_{}interaction_1_map_graph_function_bee_all_above_all_yield.svg'.format(
        interested_label, interested_alternative))
plt.close()


def get_coordinates_in_circle(n):
    return_list = []
    for i in range(n):
        theta = float(i) / n * 2 * 3.141592654
        x = np.cos(theta)
        y = np.sin(theta)
        return_list.append((x, y))
    return return_list



"this works to create video animation"

from matplotlib.animation import FFMpegWriter

metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=0.5, metadata=metadata)

#put label_abiotic
interested_label = 10
interested = frame_int_all[frame_int_all.label_abiotic == interested_label]

fig, ax = plt.subplots(1, 2, figsize=(12, 12), gridspec_kw={'width_ratios': [1.5, 10]})
with writer.saving(fig, "writer_test_abiotic_161022.mp4", 100):
    for i in range(interested.shape[0]):


        ax[0].clear()
        ax[1].clear()
        interested_alternative = i
        intr = interested.iloc[interested_alternative]
        frame_s = frame_dime[
            (frame_dime.model == 'Snodgrassella_alvi_wkB2GF') & (frame_dime.yield_perc == intr.Salvi_yield) & (
                    frame_dime.alternative == intr.Salvi_alt)]
        frame_g = frame_dime[
            (frame_dime.model == 'Gilliamella_apicola_wkB1GF') & (frame_dime.yield_perc == intr.Gapicola_yield) & (
                    frame_dime.alternative == intr.Gapicola_alt)]
        frame_fig = frame_s.append(frame_g)
        print(frame_fig)
        get_interaction_map_video(frame_fig, positions, intr,i+1,background_black=True)

        writer.grab_frame()

        # plt.close()
        # frames.append(fig)



"for all sizes of positive interaction 2 species model"

frame_int_all = frame_int
add_label(frame_int_all, "pos_int")

frame_int_all=frame_int_all.sort_values('label_pos_int')
frame_int_all['len_pos_int']=[len(k) for k in frame_int_all.pos_int]
frame_int_all=frame_int_all.sort_values(['len_pos_int','label_pos_int'])

interested_label = 0

interested = frame_int_all[frame_int_all.label_pos_int == interested_label]
alts=np.linspace(0,1,interested.shape[0])
interested_alternative = int(alts[0])
intr = interested.iloc[interested_alternative]
frame_s = frame_dime[(frame_dime.model == 'Snodgrassella_alvi_wkB2GF') & (frame_dime.yield_perc == intr.Salvi_yield) & (
            frame_dime.alternative == intr.Salvi_alt)]
frame_g = frame_dime[
    (frame_dime.model == 'Gilliamella_apicola_wkB1GF') & (frame_dime.yield_perc == intr.Gapicola_yield) & (
                frame_dime.alternative == intr.Gapicola_alt)]
frame_fig = frame_s.append(frame_g)
fig = get_interaction_map_positive_int(frame_fig, positions, intr)
plt.savefig(
    output_file_figures + '/POS_INT_LABEL_{}_ALTERNATIVE_{}interaction_1_map_graph_function_bee_all_above_all_yield.svg'.format(
        interested_label, interested_alternative))
plt.close()




from .figure_functions import *
interested_label=[k for k in frame_int_all.label_pos_int.unique()]
from matplotlib.animation import FFMpegWriter

metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=0.5, metadata=metadata)

fig, ax = plt.subplots(1, 2, figsize=(20, 20), gridspec_kw={'width_ratios': [1, 10]})
with writer.saving(fig, "writer_test_pos_int_lol.mp4", 100):
    for lab in interested_label:
        interested = frame_int_all[frame_int_all.label_pos_int == lab]
        for i in range(interested.shape[0]):
            ax[0].clear()
            ax[1].clear()
            interested_alternative = i
            intr = interested.iloc[interested_alternative]
            frame_s = frame_dime[
                (frame_dime.model == 'Snodgrassella_alvi_wkB2GF') & (frame_dime.yield_perc == intr.Salvi_yield) & (
                        frame_dime.alternative == intr.Salvi_alt)]
            frame_g = frame_dime[
                (frame_dime.model == 'Gilliamella_apicola_wkB1GF') & (frame_dime.yield_perc == intr.Gapicola_yield) & (
                        frame_dime.alternative == intr.Gapicola_alt)]
            frame_fig = frame_s.append(frame_g)
            print(frame_fig)
            get_interaction_map_positive_int(frame_fig, positions, intr)

            writer.grab_frame()

        # plt.close()
        # frames.append(fig)


#to only highlight positive interactions
fig, ax = plt.subplots(1, 2, figsize=(20, 20), gridspec_kw={'width_ratios': [1, 10]})
with writer.saving(fig, "writer_test_pos_int_only.mp4", 100):
    for lab in interested_label:
        interested = frame_int_all[frame_int_all.label_pos_int == lab]
        for i in range(interested.shape[0]):
            ax[0].clear()
            ax[1].clear()
            interested_alternative = i
            intr = interested.iloc[interested_alternative]
            frame_s = frame_dime[
                (frame_dime.model == 'Snodgrassella_alvi_wkB2GF') & (frame_dime.yield_perc == intr.Salvi_yield) & (
                        frame_dime.alternative == intr.Salvi_alt)]
            frame_g = frame_dime[
                (frame_dime.model == 'Gilliamella_apicola_wkB1GF') & (frame_dime.yield_perc == intr.Gapicola_yield) & (
                        frame_dime.alternative == intr.Gapicola_alt)]
            frame_fig = frame_s.append(frame_g)

            rows_pos_int=['EX_{}'.format(k) for k in intr.pos_int]
            frame_fig=frame_fig[frame_fig.metabolites.isin(rows_pos_int)]
            print(frame_fig)
            get_interaction_map_positive_int(frame_fig, positions, intr)

            writer.grab_frame()

        # plt.close()
        # frames.append(fig)

"for pairwise analysis"
""" this is for the matrix """
"first create the data frame with nan values 1-1"
d=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_090922_combin2.csv",index_col=0)
d3=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_090922_combin3.csv",index_col=0)
d4=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_090922_combin4.csv",index_col=0)
d5=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_090922_combin5.csv",index_col=0)
d6=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_090922_combin6.csv",index_col=0)
d7=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_090922_combin7.csv",index_col=0)



d=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_080323_combin2.csv",index_col=0)
d3=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_080323_combin3.csv",index_col=0)
d4=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_080323_combin4.csv",index_col=0)
d5=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_080323_combin5.csv",index_col=0)
d6=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_080323_combin6.csv",index_col=0)
d7=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_080323_combin7.csv",index_col=0)




d=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_111023_combin2.csv",index_col=0)
d3=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_111023_combin3.csv",index_col=0)
d4=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_111023_combin4.csv",index_col=0)
d5=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_111023_combin5.csv",index_col=0)
d6=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_111023_combin6.csv",index_col=0)
d7=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_111023_combin7.csv",index_col=0)


import ast #toconvert the string to list

dall=pd.concat([d,d3,d4,d5,d6,d7])
dall['models']=[ ast.literal_eval(k) for k in dall.pair]
dall['pair_len']=[len(k) for k in dall.models]
interested_cols=['min_competition', 'max_cooperation', 'min_abiotic sources']
interested_col=interested_cols[2]
d['models']=[ ast.literal_eval(k) for k in d.pair]


plt.rcParams.update({'font.size':16,
                      'font.family':'Calibri',
                     'legend.fontsize':16,
                   'xtick.labelsize' : 16,
                   'ytick.labelsize' : 16
  })
plt.rcParams['svg.fonttype'] = 'none'

# models_int=['Snodgrassella_alvi_wkB2GF','Gilliamella_apicola_wkB1GF']
# mod_of_int=[k for k in range(dall.shape[0]) if len(list(set(dall.models.iloc[k] ) & set(models_int)))==2 ]
#
# dall=dall.iloc[mod_of_int]
#consider only the ones have snod and gilliamella
plt.figure()
plt.scatter(dall['max_cooperation'],dall['min_abiotic sources'],c=dall.pair_len,cmap=cm.rainbow,alpha=0.5)
plt.xlabel('max cooperation')
plt.ylabel('min abiotic')
# plt.yscale('log')
# plt.xscale('log')
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file_figures+'/061123_241023_pairwise_study.svg')

plt.savefig(output_file_figures+'/270823_pairwise_study.svg')

#better version of figure




dall=dall.sort_values('min_abiotic sources')
dall_selected=dall.groupby('pair_len').first()
dall_selected=dall.groupby('pair_len').apply(lambda gr: gr[gr['min_abiotic sources']==gr['min_abiotic sources'].min()])


fig, ax = plt.subplots(figsize=(8, 6))
no_categ=6
heatmap = cm.get_cmap('rainbow', lut=no_categ)
cax=ax.scatter(dall['max_cooperation'],dall['min_abiotic sources'],c=dall.pair_len,cmap=heatmap,alpha=0.5,s =70)
# cax=ax.scatter(dall['max_cooperation'],dall['min_competition'],c=dall.pair_len,cmap=heatmap,alpha=0.5,s =70)
# ax.scatter(dall_selected['max_cooperation'],dall_selected['min_abiotic sources'],facecolor='none',edgecolor='black',s=180,linewidth=2)
ax.set_xlabel('maximum cross-fed metabolites')
ax.set_ylabel('minimal nutritional requirements')
# ax.set_ylabel('minimum competition')

# plt.yscale('log')
# plt.xscale('log')
# ax.colorbar()
cbar = fig.colorbar(cax, ax=ax)
ax.set_xticks(np.arange(3, 22, step=2))
plt.xlim([2,22])
plt.tight_layout()
plt.savefig(output_file_figures+'/061123_111023_pairwise_study_coop_abiot.svg')

plt.savefig(output_file_figures+'/270823_pairwise_study_mincompet.svg')
plt.close()
# ,s=[k+50 for k in dall_selected.index.tolist()]
# s=[k+30 for k in dall.index.tolist()]



#fetch data for bar chart
path='/remind/models/bee_models/core_members/tfa_real_010922'

data_dir_cobra='/remind/models/bee_models/core_members/carveme_models'

allFiles = glob.glob(path + "/*"+"json")
model_name=[]
genes=[]
rxns=[]
mets=[]
gene_assoc_rxns=[]
rxn_assoc_genes=[]
for k in range(len(allFiles)):
    model = load_json_model(allFiles[k])
    model_name+=[(model.name.split('tfa_')[1])]
    model=cobra.io.read_sbml_model(join(data_dir_cobra, "{}.xml".format((model.name.split('tfa_')[1]))))
    genes+=[len(model.genes)]
    rxns+=[len(model.reactions)]
    mets+=[len(model.metabolites)]
    gene_assoc_rxns += [len([k for k in model.reactions if len(k.gene_reaction_rule) >= 1])]
    rxn_assoc_genes += [len([k for k in model.genes if len(k.reactions)>=1])]

data_models=pd.DataFrame()
data_models['model']=model_name
data_models['genes']=genes
data_models['rxns']=rxns
data_models['metabolites']=mets
data_models['gene_assoc_rxns']=gene_assoc_rxns
data_models['rxn_assoc_genes']=rxn_assoc_genes



alts_no=frame_combined_all.groupby(['model','yield_perc']).alternative.nunique().reset_index(name='altsno')
alts_no.groupby('model').altsno.sum()

model_dict={'Bifidobacterium_asteroides_PRL2011GF': "B.asteroides",
        'Gilliamella_apicola_wkB1GF':"G.apicola",
        'Lactobacillus_apis_Hma11GF':"L.apis",
        'Lactobacillus_kullabergensis_Biut2GF':"L.kullabergensis",
        'Lactobacillus_mellifer_Bin4': "L.mellifer",
        'Lactobacillus_mellis_Hon2':"L.mellis",
        'Snodgrassella_alvi_wkB2GF':"S.alvi"}

fig, ax = plt.subplots(figsize=(8, 6))
data_models_sorted=data_models.sort_values(['genes'],ascending=False)
ticks=[model_dict[k] for k in data_models_sorted.model.unique()]
ax=data_models_sorted.plot.bar(x='model',y=['genes'])
ax.set_xticklabels(ticks,rotation=90,style='italic')
plt.tight_layout()
plt.savefig(output_file_figures+'/barchart_genes.svg')

fig, ax = plt.subplots(figsize=(8, 6))
data_models_treated=data_models_sorted[['model','gene_assoc_rxns','genes']]
data_models_treated=data_models_treated.set_index('model')
data_models_treated=data_models_treated.sort_values(['genes'],ascending=False)
ticks=[model_dict[k] for k in data_models_treated.index.to_list()]
ax=data_models_treated.plot.bar()
ax.set_xticklabels(ticks,rotation=90,style='italic')
plt.tight_layout()
plt.savefig(output_file_figures+'/barchart_genes_and_assoc_rxns.svg')






fig, ax = plt.subplots(figsize=(8, 6))
no_categ=6
heatmap = cm.get_cmap('rainbow', lut=no_categ)
ax.set_xlabel('maximum cross-fed metabolites')
ax.set_ylabel('minimal nutritional requirements')
# plt.yscale('log')
# plt.xscale('log')
# ax.colorbar()
cbar = fig.colorbar(cax, ax=ax)
ax.set_xticks(np.arange(2, 18, step=2))
plt.tight_layout()
plt.savefig(output_file_figures+'/pairwise_study.svg')
plt.close()





plt.figure()
plt.scatter(dall['max_cooperation'],dall['min_competition'],c=dall.pair_len,cmap=cm.rainbow,alpha=0.5)
plt.xlabel('max cooperation')
plt.ylabel('min competition')
# plt.yscale('log')
# plt.xscale('log')
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file_figures+'/pairwise_study_2.svg')

dall=dall[dall.pair_len==2]
plt.figure()
plt.scatter(dall['max_cooperation'],dall['min_competition'],cmap=cm.rainbow,alpha=0.5)
plt.xlabel('max cooperation')
plt.ylabel('min competition')
# plt.yscale('log')
# plt.xscale('log')
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file_figures+'/pairwise_study_3.svg')

dall_2_3=dall[dall.pair_len.isin([2,3])]
dall_3=dall[dall.pair_len.isin([3])]
dall_4=dall[dall.pair_len.isin([4])]
dall_5=dall[dall.pair_len.isin([5])]

dall_2=dall[dall.pair_len.isin([2])]
dall_3['models_tup'] = [tuple(k) for k in dall_3.models]
dall_4['models_tup'] = [tuple(k) for k in dall_4.models]
dall_5['models_tup'] = [tuple(k) for k in dall_5.models]

first=[]
pair1=[]
pair2=[]
pair3=[]
pair4=[]
model=[]

int_col='min_abiotic sources'

# int_col='max_cooperation'
def stratify_data(dall_4,dall_3,limit=4.0):
    first = []
    pair1 = []
    pair2 = []
    pair3 = []
    pair4 = []
    pair5=[]
    model = []
    for mod in dall_4.models:
        print([k for k in dall_3.models if len(list(set(k) & set(mod)))>=limit ])

        int=[tuple(k) for k in dall_3.models if len(list(set(k) & set(mod)))>=limit ]

        dall_3['models_tup']=[tuple(k) for k in dall_3.models]
        dall_3[dall_3.models_tup.isin(int)]
        model_int=[tuple(mod)]
        data=dall_3[dall_3.models_tup.isin(int)][int_col]
        pair1.append(data.iloc[0])
        pair2.append(data.iloc[1])
        pair3.append(data.iloc[2])
        pair4.append(data.iloc[3])
        pair5.append(data.iloc[4])
        model.append(model_int)
        first.append(dall_4[dall_4.models_tup.isin(model_int)][int_col].iloc[0])


    d=pd.DataFrame()
    d['model']=model
    d['orig']=first
    d['pair1']=pair1
    d['pair2']=pair2
    d['pair3']=pair3
    d['pair4']=pair4
    d['pair5']=pair5

    return d



#creating dendrogram

from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram

# d = hierarchy.distance.pdist(data_distance,mask = ~np.isfinite(data_distance))   # vector of (100 choose 2) pairwise distances


# L = sch.linkage(d, method='complete')
# Z = hierarchy.linkage(frame_int.abiotic, method='complete')
Z = hierarchy.linkage(data_distance, method='average',metric='cityblock')
# metric='cityblock')
plt.figure()
plt.title("Dendrograms")
# Dendrogram plotting using linkage matrix
dendrogram = hierarchy.dendrogram(Z,labels=[k+1 for k in data_distance.index])
plt.savefig(output_file_figures+'/trial_dendrogram_3.png')


data_dist

# #trial2 cart tree
# import matplotlib.pyplot as plt
# from sklearn.datasets import load_iris
# from sklearn.tree import DecisionTreeClassifier
# from sklearn import tree
#
# X, y = load_iris(return_X_y=True)
#
# # Make an instance of the Model
# clf = DecisionTreeClassifier()
#
# # Train the model on the data
# clf.fit(X, y)
#
# fn=['sepal length (cm)','sepal width (cm)','petal length (cm)','petal width (cm)']
# cn=['setosa', 'versicolor', 'virginica']
#
# # Setting dpi = 300 to make image clearer than default
# fig, axes = plt.subplots(nrows = 1,ncols = 1,figsize = (4,4), dpi=300)
#
# tree.plot_tree(clf,
#            feature_names = fn,
#            class_names=cn,
#            filled = True);
#
# fig.savefig(output_file_figures+'/imagename.png')


import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree

# X, y = load_iris(return_X_y=True)

# Make an instance of the Model
clf = DecisionTreeClassifier()
classes= np.linspace(0,11,12)
# Train the model on the data
clf.fit(f.values, classes)

# fn=['sepal length (cm)','sepal width (cm)','petal length (cm)','petal width (cm)']
fn=[k for k in data_distance.columns]
cn=[str(k) for k in classes]

# Setting dpi = 300 to make image clearer than default
fig, axes = plt.subplots(nrows = 1,ncols = 1,figsize = (4,4), dpi=300)

tree.plot_tree(clf,
           feature_names = fn,
           class_names=cn,
           filled = True);

fig.savefig(output_file_figures+'/lolipop_imagename.png')


cols=[k for k in data_distance.columns]

f=data_distance[cols]. apply(lambda x: x. astype('category'))


from figure_functions import *
frame_int_all=get_3_species_data(d)


frame_int = frame_int_all.groupby('abiotic').first().reset_index()

frame_int['Gapic_yield_min']=frame_int_all.groupby('abiotic').Gapicola_yield.min().values
frame_int['Gapic_yield_max']=frame_int_all.groupby('abiotic').Gapicola_yield.max().values

frame_int['Salvi_yield_min']=frame_int_all.groupby('abiotic').Salvi_yield.min().values
frame_int['Salvi_yield_max']=frame_int_all.groupby('abiotic').Salvi_yield.max().values

frame_int['Bifido_yield_min']=frame_int_all.groupby('abiotic').Bifido_yield.min().values
frame_int['Bifido_yield_max']=frame_int_all.groupby('abiotic').Salvi_yield.max().values

frame_int_all['Gapic_dimes'] = frame_int_all.Gapic_uptake + frame_int_all.Gapic_secretions
frame_int_all['Salvi_dimes'] = frame_int_all.Salvi_uptake + frame_int_all.Salvi_secretions
frame_int_all['Gapic_dimes'] = frame_int_all.Gapic_uptake + frame_int_all.Gapic_secretions
frame_int_all['Bifido_dimes'] = frame_int_all.Bifido_uptake + frame_int_all.Bifido_secretions

interested = frame_int_all[frame_int_all.label_abiotic == interested_label]
intr = interested.iloc[interested_alternative]
frame_s = frame_dime[(frame_dime.model == 'Snodgrassella_alvi_wkB2GF') & (frame_dime.yield_perc == intr.Salvi_yield) & (
            frame_dime.alternative == intr.Salvi_alt)]
frame_g = frame_dime[
    (frame_dime.model == 'Gilliamella_apicola_wkB1GF') & (frame_dime.yield_perc == intr.Gapicola_yield) & (
                frame_dime.alternative == intr.Gapicola_alt)]

frame_b = frame_dime[
    (frame_dime.model == 'Bifidobacterium_asteroides_PRL2011GF') & (frame_dime.yield_perc == intr.Bifido_yield) & (
                frame_dime.alternative == intr.Bifido_alt)]

frame_fig = frame_s.append(frame_g)
frame_fig=frame_fig.append(frame_b)
fig = get_interaction_map_3_species(frame_fig, positions, intr)
plt.savefig(
    output_file_figures + '/3_species_LABEL_{}_ALTERNATIVE_{}interaction_1_map_graph_function_bee_all_above_all_yield.svg'.format(
        interested_label, interested_alternative))
plt.close()

mets_=[]
for k in frame_int_all.Salvi_dimes.unique():
    mets_ += list(k)
for k in frame_int_all.Gapic_dimes.unique():
    mets_ += list(k)
for k in frame_int_all.Bifido_dimes.unique():
    mets_ += list(k)

mets_ = list(set(mets_))
no_mets = len(mets_)

pos = get_coordinates_in_circle(no_mets)
positions = dict(zip(mets_, pos))


#make the videp
frame_fig = frame_s.append(frame_g)
frame_fig=frame_fig.append(frame_b)
fig = get_interaction_map_3_species(frame_fig, positions, intr)
plt.savefig(
    output_file_figures + '/3_species_LABEL_{}_ALTERNATIVE_{}interaction_1_map_graph_function_bee_all_above_all_yield.svg'.format(
        interested_label, interested_alternative))
plt.close()



"3 species video"
from matplotlib.animation import FFMpegWriter

metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=0.5, metadata=metadata)

# fig, ax = plt.subplots(1, 2, figsize=(12, 12), gridspec_kw={'width_ratios': [1, 10]})
fig, ax = plt.subplots(1, 2, figsize=(20, 20), gridspec_kw={'width_ratios': [1, 10]})

with writer.saving(fig, "writer_test_video_3_soecies.mp4", 100):
    for i in range(interested.shape[0]):
        ax[0].clear()
        ax[1].clear()
        interested_alternative = i
        intr = interested.iloc[interested_alternative]
        frame_s = frame_dime[
            (frame_dime.model == 'Snodgrassella_alvi_wkB2GF') & (frame_dime.yield_perc == intr.Salvi_yield) & (
                    frame_dime.alternative == intr.Salvi_alt)]
        frame_g = frame_dime[
            (frame_dime.model == 'Gilliamella_apicola_wkB1GF') & (frame_dime.yield_perc == intr.Gapicola_yield) & (
                    frame_dime.alternative == intr.Gapicola_alt)]
        frame_b = frame_dime[
            (frame_dime.model == 'Bifidobacterium_asteroides_PRL2011GF') & (
                        frame_dime.yield_perc == intr.Bifido_yield) & (
                    frame_dime.alternative == intr.Bifido_alt)]
        frame_fig = frame_s.append(frame_g)
        frame_fig = frame_fig.append(frame_b)
        print(frame_fig)
        get_interaction_map_3_species_video(frame_fig, positions, intr)

        writer.grab_frame()



dime=("26dap__M_e", "acmana_e", "arg__L_e", "cmp_e", "glc__D_e", "glu__L_e", "his__L_e", "ins_e", "nmn_e", "pydx_e", "thm_e", "asp__L_e", "dha_e", "lac__D_e")

interested_cut = interested[interested.Gapic_dimes == dime]
interested=interested_cut





"reindex"
data_coop_label.set_index('label_abiotic', inplace=True)
frame_int.set_index('label_abiotic', inplace=True)

frame_int=frame_int.reindex(data_coop_label.index)



"Subplot 7 species"


"""SUBPLOT TRIAL for 2 species"""
all_exchanges = pd.read_csv("/remind/projects/bee_project/stats/all_exchanges_stats.csv")
all_exchanges = all_exchanges.groupby('mets').first().reset_index()
met_name = [all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols]
# now heatmap
heatmap_abiotic = "tab20b"
# heatmap_abiotic = "Accent"
# heatmap_abiotic = "Paired"
heatmap_yield = "rainbow"

# df_pivoted = data_coop

df_pivoted = data_distance
# df_pivoted = data_distance_fornow
fig, (ax, ax2, ax3,ax4,ax5,ax6,ax7,ax8) = plt.subplots(1, 8, figsize=(80, 20)) #20 was10 before

# fig, ax = plt.subplots(figsize=(30, 10))
no_categ = 5
# no_categ = 5
heatmap = cm.get_cmap(heatmap_abiotic, lut=no_categ)
# heatmap = cm.Greens
alpha_chosen = 0.7
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
frame_int=frame_int_unique
# frame_int=fr
df_pivoted = frame_int[['Gapicola_yield_min', 'Gapicola_yield_max']]

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
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax3)
# plt.tight_layout()
# plt.savefig(output_file_figures + '/2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
# fig.tight_layout()
# plt.close('all')


# forth figure
df_pivoted = frame_int[['Bifido_yield_min', 'Bifido_yield_max']]
no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.summer
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax4.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax4.xaxis.tick_top()
ax4.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax4.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax4.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
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

ax4.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax4.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax4.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax4)
# plt.tight_layout()
# plt.savefig(output_file_figures + '/2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
# fig.tight_layout()
# plt.close('all')


# fifth figure
df_pivoted = frame_int[['Lapis_yield_min', 'Lapis_yield_max']]
no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.summer
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax5.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax5.xaxis.tick_top()
ax5.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax5.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax5.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
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

ax5.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax5.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax5.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax5)
# plt.tight_layout()
# plt.savefig(output_file_figures + '/2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
# fig.tight_layout()
# plt.close('all')



# sixth figure
df_pivoted = frame_int[['Lkulla_yield_min', 'Lkulla_yield_max']]
no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.summer
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax6.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax6.xaxis.tick_top()
ax6.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax6.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax6.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
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

ax6.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax6.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax6.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax6)
# plt.tight_layout()
# plt.savefig(output_file_figures + '/2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
# fig.tight_layout()
# plt.close('all')


# seventh figure
df_pivoted = frame_int[['Lmellifer_yield_min', 'Lmellifer_yield_max']]
no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.summer
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax7.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax7.xaxis.tick_top()
ax7.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax7.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax7.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
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

ax7.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax7.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax7.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax7)
# plt.tight_layout()
# plt.savefig(output_file_figures + '/2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
# fig.tight_layout()
# plt.close('all')

# eigth figure
df_pivoted = frame_int[['Lmellis_yield_min', 'Lmellis_yield_max']]
no_categ = 11
# heatmap = cm.get_cmap("tab20b", lut = no_categ)
# heatmap = cm.summer
heatmap = cm.get_cmap(heatmap_yield, lut=no_categ)

current_cmap = matplotlib.cm.get_cmap(heatmap)
current_cmap.set_bad(color='white')
# cax2=ax.imshow(df_pivoted,cmap=current_cmap,vmin=df_pivoted.min().min(),vmax=df_pivoted.max().max(),alpha=alpha_chosen)
# for the yield only
cax2 = ax8.imshow(df_pivoted, cmap=current_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)

# cax_cont=ax.imshow(contr,cmap=heatmap,vmin=contr.min(),vmax=contr.max(),alpha=1.0)

ax8.xaxis.tick_top()
ax8.set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                           df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
ax8.set_yticks(np.linspace(0,df_pivoted.index.nunique()-1,df_pivoted.index.nunique()))#,[str(k) for k in df_pivoted.columns.to_list()])
# ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
# trial
ax8.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
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

ax8.set_xticks(np.arange(-.5, df_pivoted.columns.nunique() - 1, 1), minor=True)
ax8.set_yticks(np.arange(-.5, df_pivoted.index.nunique() - 1, 1), minor=True)
ax8.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
# ticks_ = np.linspace(0, 6,7, endpoint=True)
cbar = fig.colorbar(cax2, ax=ax8)
plt.tight_layout()
plt.savefig(output_file_figures + '/7_member_pos_int_2_all_pos_int_interactions_same_col_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
fig.tight_layout()
plt.close('all')



