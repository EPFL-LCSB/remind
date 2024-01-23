import pandas as pd
import numpy as np
from figure_functions import *


# data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'
# frame_dime = pd.read_hdf(data_path)

df=pd.read_hdf('/remind/projects/bee_project/ilp_solutions_070823_nplusone_2_member_bee/obj_num_2_modelofint_6_alternatives_interact\
ion_positive/sols_active_alternative_for_diff_obj_2_alternation_interaction_positive_altno_1283_.h5')

df=pd.read_hdf("/remind/projects/bee_project/131023_ilp_solutions_280823_nplusone_2_member_bee/\
obj_num_6_modelofint_0_alternatives_abiotic_sources_yield/sols_active_alternative_for_diff_obj_6_alternation_abiotic_sources_yield_altno_3_.h5")


df_un=pd.read_hdf("/remind/projects/bee_project/131023_ilp_solutions_280823_nplusone_2_member_bee/obj_num_6_modelofint_0_alternatives_abiotic_sources/sols_active_alternative_for_diff_obj_6_alternation_abiotic_sources_altno_1_.h5")


models_dict={
             "model1":"Gapicola",
             "model2":"Salvi",
            "model3":"Bifido"}

frame_int=get_species_data(df,models_dict=models_dict,len_models=3)


#get d from the other script

def add_label(df, group=[""], overshoot=0):
    df['label_{}'.format(group)] = df.groupby(by=group).grouper.group_info[0] + overshoot

add_label(frame_int, "abiotic")
frame_int_all = frame_int
add_label(frame_int_all, "abiotic")

# add label


label = False
f_ = []
for k in frame_int.abiotic.unique():
    f_ += list(k)

f__ = pd.DataFrame(f_).value_counts()

cols = [k[0] for k in f__.index]
if label:
    cols_all = cols + ['label_abiotic']
else:
    cols_all = cols

""" first to get the number of uptakes"""
# cols = [k[0] for k in f__.index]
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

helper=data_coop

cols_all = cols
dummyarray = np.empty((frame_int_all.alternative.nunique(), len(cols_all)))
dummyarray[:] = np.nan
data_coop = pd.DataFrame(dummyarray, columns=cols_all)
for k in range(frame_int_all.alternative.nunique()):
    # label_abiot = frame_int.iloc[k].label_abiotic
    for met in cols:
        # now helper is added
        if helper.iloc[k][met] == 1.0:
            if (met in frame_int.Salvi_uptake.iloc[k]):
                if not (met in frame_int.pos_int.iloc[k]):
                    data_coop.iloc[k][met] = 0.4
            if met in frame_int.Gapicola_uptake.iloc[k]:
                if not (met in frame_int.pos_int.iloc[k]):
                    data_coop.iloc[k][met] = 0.3
            if met in frame_int.Bifido_uptake.iloc[k]:
                if not (met in frame_int.pos_int.iloc[k]):
                    data_coop.iloc[k][met] = 0.2
        if helper.iloc[k][met] == 2.0:
            if (met in frame_int.Gapicola_uptake.iloc[k]) and (met in frame_int.Salvi_uptake.iloc[k]):
                data_coop.iloc[k][met] = 0.5

            if (met in frame_int.Gapicola_uptake.iloc[k]) and (met in frame_int.Bifido_uptake.iloc[k]):
                data_coop.iloc[k][met] = 0.6

            if (met in frame_int.Salvi_uptake.iloc[k]) and (met in frame_int.Bifido_uptake.iloc[k]):
                data_coop.iloc[k][met] = 0.7

        if helper.iloc[k][met] == 3.0:
                data_coop.iloc[k][met] =0.8
        if math.isnan(helper.iloc[k][met]) :
                data_coop.iloc[k][met] = np.nan

data_coop['label_abiotic']=frame_int['label_abiotic']


"""here coloring is done if there is unique value colored accordingly if multiple all colored one"""
dummyarray2 = np.empty((frame_int_all.label_abiotic.nunique(), len(cols_all+['label_abiotic'])))
dummyarray2[:] = np.nan
data_coop_abiot_3 = pd.DataFrame(dummyarray2, columns=cols_all+['label_abiotic'])
cols_count=[]
for label in data_coop.label_abiotic.unique():
    d=data_coop[data_coop.label_abiotic==label]
    for col in d:
        col_count=len(d[col].unique())
        cols_count+=([col_count])
        label_abiot = d['label_abiotic'].unique()[0]
        data_coop_abiot_3.iloc[label_abiot]['label_abiotic'] = label_abiot
        if col_count<=1:
            if col != 'label_abiotic':
                data_coop_abiot_3.iloc[label_abiot][col]=d[col].unique()[0]
                print(col,d[col].unique()[0])
        else:
            data_coop_abiot_3.iloc[label_abiot][col] = 0.1
            # data_coop_abiot_3.iloc[k]['met']=d[col].unique()[0]
        # print(len(d[col].unique()))
    # print(col,d[col].unique())

np.unique(data_coop_abiot_3[cols_all].values)

# for col in data_coop:
#     print(col,data_coop[col].unique())

"""end"""
frame_int = frame_int_all.groupby('label_abiotic').first().reset_index()

frame_int['Gapic_yield_min']=frame_int_all.groupby('abiotic').Gapicola_yield.min().values
frame_int['Gapic_yield_max']=frame_int_all.groupby('abiotic').Gapicola_yield.max().values

frame_int['Salvi_yield_min']=frame_int_all.groupby('abiotic').Salvi_yield.min().values
frame_int['Salvi_yield_max']=frame_int_all.groupby('abiotic').Salvi_yield.max().values

frame_int['Bifido_yield_min']=frame_int_all.groupby('abiotic').Bifido_yield.min().values
frame_int['Bifido_yield_max']=frame_int_all.groupby('abiotic').Bifido_yield.max().values #was erroneous



"data distance for heatmap"
"try also for three-members "
"this is just to put everything to 1 for the distance matrix"
frame_int_unique=frame_int_all.groupby('abiotic').first().reset_index()
dummyarray = np.empty((frame_int_unique.alternative.nunique(), len(cols)))
dummyarray[:] = np.nan
data_distance = pd.DataFrame(dummyarray, columns=cols_all)
for k in range(frame_int_unique.alternative.nunique()):
    for met in cols:
        if (met in frame_int_unique.abiotic.iloc[k]):
            data_distance.iloc[k][met] = int(1)
        #
        # else:
        #     data_distance.iloc[k][met] = int(0)
size = data_distance.shape[0]
data_distance_ap = data_distance.append(data_distance.where(data_distance.notna()).sum(), ignore_index=True)
data_distance_ap = data_distance_ap.sort_values(by=size, axis=1, ascending=False)
data_distance = data_distance_ap.drop(index=size, axis=1)
cols = [k for k in data_distance.columns]

data_distance_label=data_distance.copy()
data_distance_label['label_abiotic']=frame_int_all.groupby('abiotic').first()['label_abiotic'].values

"""new"""
"reindex"

data_coop_plot=data_coop_abiot_3
data_coop_plot.set_index('label_abiotic', inplace=True)
data_coop_plot=data_coop_plot[cols]
data_coop_plot=data_coop_plot.reindex(data_distance_label.index)

frame_int.set_index('label_abiotic', inplace=True)
frame_int=frame_int.reindex(data_distance_label.index)


"""new 090823"""
new_cols=['xylan4_e',
 'ribflv_e',
 'gln__L_e',
'ala_B_e',
'pydx_e',
'acmana_e',
 'nmn_e',
 'his__L_e',
 '26dap__M_e',
 'cmp_e',
'thm_e',
'arg__L_e',
 'fol_e',

 'ade_e',
 'ins_e',
 'hxan_e',
 'glu__L_e',
 'akg_e',
 'fru_e',
 'glc__D_e',
 'sucr_e',
 'dha_e']

data_coop_plot=data_coop_plot[new_cols]
cols=new_cols

"data distance for dendrogram"
"try also for three-members "
"this is just to put everything to 1 for the distance matrix"
frame_int_unique=frame_int_all.groupby('abiotic').first().reset_index()
dummyarray = np.empty((frame_int_unique.alternative.nunique(), len(cols_all)))
dummyarray[:] = np.nan
data_distance_dend = pd.DataFrame(dummyarray, columns=cols_all)
for k in range(frame_int_unique.alternative.nunique()):
    for met in cols:
        if (met in frame_int_unique.abiotic.iloc[k]):
            data_distance_dend.iloc[k][met] = int(1)
        #
        else:
            data_distance_dend.iloc[k][met] = int(0)



from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram
Z = hierarchy.linkage(data_distance_dend, method='average',metric='cityblock')
plt.figure()
plt.title("Dendrograms")
dendrogram = hierarchy.dendrogram(Z,labels=[k+1 for k in data_distance_dend.index])
plt.savefig(output_file_figures+'/trial_dendrogram_3.svg')

index_new=dendrogram['leaves']
#put this to abiotic
data_distance_label=data_distance_label.iloc[index_new].reset_index().drop('index',axis=1)

data_distance=data_distance_label.drop('label_abiotic',axis=1)


"reindex"
data_distance_label.set_index('label_abiotic', inplace=True)
frame_int.set_index('label_abiotic', inplace=True)
frame_int=frame_int.reindex(data_distance_label.index)
# data_distance_dend=data_distance_dend.iloc[index_new].reset_index().drop('index',axis=1)

"subplot three species"


"""SUBPLOT TRIAL for 2 species"""
all_exchanges = pd.read_csv("/remind/projects/bee_project/stats/all_exchanges_stats.csv")
all_exchanges = all_exchanges.groupby('mets').first().reset_index()
met_name = [all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols]
# now heatmap
heatmap_abiotic = "tab20b"
heatmap_yield = "rainbow"


df_pivoted = data_distance


#plot on ax0 dendrogram

Z = hierarchy.linkage(data_distance_dend, method='average',metric='cityblock')
plt.figure()
plt.title("Dendrograms")
dendrogram = hierarchy.dendrogram(Z,labels=[k+1 for k in data_distance_dend.index])
plt.savefig(output_file_figures+'/trial_dendrogram_3.svg')



fig, (ax, ax2, ax3,ax4) = plt.subplots(1, 4, figsize=(40, 10))
# fig, ax = plt.subplots(figsize=(30, 10))
no_categ = 5
no_categ = 8
df_pivoted=data_coop_plot
heatmap = cm.get_cmap("Paired", lut=no_categ)
# heatmap = cm.Greens
alpha_chosen = 0.6
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

ax2.set_xticklabels([str(k.split('_')[0]) for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12,
                    style='italic')

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
cbar = fig.colorbar(cax2, ax=ax3)


"forth figure"

# third figure
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
cbar = fig.colorbar(cax2, ax=ax4)

plt.tight_layout()
plt.savefig(output_file_figures + '/061123_281023_090823_3_new_member_subplot_nn_7230922_2_member_min_abiot_{}__{}.svg')
plt.close('all')








interested_label=0
interested_alternative=0
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
    output_file_figures + '/061123_3_species_LABEL_{}_ALTERNATIVE_{}interaction_1_map_graph_function_bee_all_above_all_yield.svg'.format(
        interested_label, interested_alternative))
plt.close()
