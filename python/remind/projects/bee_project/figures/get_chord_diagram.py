import pandas as pd
import itertools
import matplotlib.pyplot as plt
import random
from remind.utils.postprocessing import get_unique_mets,get_substrates,get_products
import os
from itertools import combinations
from matplotlib import cm
import math
# data_path = "/remind/projects/bee_project/stats/combined_dimes_unique_180522_bee.h5"
data_path='/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'
data_path='/remind/projects/bee_project/stats/combined_dimes_unique_141122_bee.h5'


data_path='/remind/projects/bee_project/stats/combined_dimes_unique_111023_bee.h5'


plt.rcParams['svg.fonttype'] = 'none'
frame_combined_all=pd.read_hdf(data_path)


frame_combined=frame_combined_all.groupby('model').apply(lambda gr:get_unique_mets(gr)).reset_index(name='dimes')

gr='model'
frame_dimes_all = frame_combined_all.groupby(gr).apply(lambda gr: get_unique_mets(gr)).reset_index(name='dimes')
frame_dimes_subs = frame_combined_all.groupby(gr).apply(lambda gr: get_substrates(gr)).reset_index(name='dimes_subs')
frame_dimes_prods = frame_combined_all.groupby(gr).apply(lambda gr: get_products(gr)).reset_index(name='dimes_prods')

frame_dimes = frame_dimes_all
frame_dimes['dimes_subs'] = frame_dimes_subs['dimes_subs']
frame_dimes['dimes_prods'] = frame_dimes_prods['dimes_prods']

#get all pairs
all_pairs=list(combinations(frame_combined.model.unique(), 2))
# sample_dict={2:2000,
#              3:2000,
#              4:2000,
#              5:2000,
#              6:2000,
#              7:2000,}
#do it for each pair

counts = []
counts_subs = []
counts_prods = []
models=[]
counts_scaled=[]
counts_subs_scaled=[]
counts_prods_scaled=[]
count_unique_subs=[]
for index_pair in range(len(all_pairs)):
    models_of_int=list(all_pairs[index_pair])
    print("performing comparison for the community {}".format(models_of_int))
    """get the frame of interest"""
    frame_of_int=frame_dimes[frame_dimes.model.isin(models_of_int)]
    #
    # frame_dimes=frame_dimes_all
    # frame_dimes['dimes_subs']=frame_dimes_subs['dimes_subs']
    # frame_dimes['dimes_prods']=frame_dimes_prods['dimes_prods']

    no_models=frame_of_int.model.nunique()
    # possible_combin=
    list_dimes=[]
    index_dimes=[]
    d=pd.DataFrame()
    list_subs = []
    list_prods = []
    # counts=[]
    # counts_subs=[]
    # counts_prods=[]
    threshold = 2
    # n_sample=sample_dict[no_models]
    for model in frame_of_int.model.unique():
        f=frame_of_int.groupby('model').get_group(model)
        list_dimes+=f.dimes.iloc[0]
        list_subs += f.dimes_subs.iloc[0]
        list_prods += f.dimes_prods.iloc[0]
        #do s to p and p to s comparison


    list_mets_stats = pd.Series(list_dimes).value_counts()
    list_subs_stats = pd.Series(list_subs).value_counts()
    list_prods_stats = pd.Series(list_prods).value_counts()

    counts.append(list_mets_stats[list_mets_stats >= threshold].shape[0])
    counts_subs.append(list_subs_stats[list_subs_stats >= threshold].shape[0])
    counts_prods.append(list_prods_stats[list_prods_stats >= threshold].shape[0])
    counts_scaled.append(list_mets_stats[list_mets_stats >= threshold].shape[0]/len(list(set(list_dimes))))
    counts_subs_scaled.append(list_subs_stats[list_subs_stats >= threshold].shape[0]/len(list(set(list_subs))))
    counts_prods_scaled.append(list_prods_stats[list_prods_stats >= threshold].shape[0]/len(list(set(list_prods))))
    count_unique_subs.append(len(list(set(list_subs))))
    models.append(models_of_int)

        # idx=f.index.to_list()
        # index_dimes.append(idx)
        # d[model+'_index']=random.choices(idx,k=n_sample)

d=pd.DataFrame()
d['models']=models
d['count_mets_scaled']=counts_scaled
d['count_mets']=counts
d['count_subs_scaled']=counts_subs_scaled
d['count_prods_scaled']=counts_prods_scaled
d['count_prods']=counts_prods
d['model1']=[k[0] for k in d.models]
d['model2']=[k[1] for k in d.models]
d['count_subs']=counts_subs
d['count_unique_subs']=count_unique_subs
d['to_the_chord']=[round(k*100) for k in counts_subs_scaled]

scaling_factor=100

# d['to_the_chord']=[round(k*scaling_factor) for k in counts_subs]
# scaling_factor=5
#d['to_the_chord']=[round(k*100) for k in counts_subs_scaled]

chord_frame=d[['model1','model2','to_the_chord']]
# from holoviews.plotting.mpl.element
from holoviews.plotting.mpl import element
import holoviews as hv
from holoviews import opts, dim
from holoviews.util import opts

cities = list(set(chord_frame["model1"].unique().tolist() + chord_frame["model2"].unique().tolist()))
cities_dataset = hv.Dataset(pd.DataFrame(cities, columns=["model_names"]))
#define colors for the nodes
node_colors=['indigo','blueviolet','darkgreen','limegreen','sandybrown','mediumpurple','darkgoldenrod']
hv.extension("matplotlib")
fig,ax=plt.subplots()
colormap='coolwarm'
alpha_chosen=0.8
#['blue','red','orange','grey','black','yellow','pink']
chord = hv.Chord((chord_frame,cities_dataset))#.opts(colorbar=True, xaxis=None, yaxis=None)
bins=cities
label=list(np.linspace(0,1.0,7))
chord.opts(opts.Chord(labels='model_names',
                      # cmap='Tab10', #Category20'
                    edge_cmap=colormap,
                   # edge_color=dim('to_the_chord').norm(),
                      edge_color=dim('to_the_chord'),
                      edge_linewidth=3,
                      node_size=70,
                       clim=(chord_frame.to_the_chord.min(),chord_frame.to_the_chord.max()),
                      # colorbar=True,
                      # cbar_extend='both',
                      # colorbar_opts={cmap:'coolwarm',clim:(25,51)},
                        #node_linewidth=30,
                      edge_alpha=alpha_chosen,
                      normalize=True,
                      padding=500,
                      show_legend=True,
                      symmetric=True,
                    node_color=node_colors
                      # node_color=dim('model_names').str()
)
)

hv.output(fig='svg', size=800)
# plt.tight_layout()
# fig.colorbar(cax)
#Create another chord plot and add the colorbar



scaling_factor=1
min_value=chord_frame.to_the_chord.min()
max_value=chord_frame.to_the_chord.max()
fig1,ax=plt.subplots()
fig1 = hv.render(chord)
img = ax.imshow(np.array([[chord_frame.to_the_chord.min()/scaling_factor,chord_frame.to_the_chord.max()/scaling_factor]]), cmap=colormap,alpha=alpha_chosen,vmin = min_value, vmax = max_value,)
img.set_visible(False)
cbar = fig1.colorbar(img,shrink=0.55)
cbar.ax.tick_params(labelsize=50,width=5,length=10)
cbar.set_ticks(np.arange(min_value, max_value, 5))
# cbar.set_clim(20, 61)
plt.tight_layout()
output_file_figures='/remind/projects/bee_project/figures'
output_file_figures='/remind/projects/bee_project/constrained_medium/analysis/figures_trial'
output_file_figures='/remind/projects/bee_project/constrained_medium/analysis/figures_pray'

# fig1.savefig(output_file_figures+'/newest_130622_new_cmap_chord_trial_plot_not_scaled_substrates.svg')
fig1.savefig(output_file_figures+'/231023_scaled_chord_TRY_PAD_new_cmap_pivoted_common_mets_competition.svg')



#to save only the chord
# hv.save(chord, output_file_figures+'/new_130622_chord_trial_plot_all_scaled_allmets.svg', fmt='svg')



""" this is for the matrix """
"first create the data frame with nan values 1-1"
d=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_090922_combin2.csv",index_col=0)
d=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/pairwise_studies/built_bee_pairwise_objectives_080323_combin2.csv",index_col=0)

d3=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/built_bee_pairwise_objectives_090922_combin3.csv",index_col=0)
d4=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/built_bee_pairwise_objectives_090922_combin4.csv",index_col=0)
d5=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/built_bee_pairwise_objectives_090922_combin5.csv",index_col=0)
d6=pd.read_csv("/remind/projects/bee_project/constrained_medium/analysis/data/built_bee_pairwise_objectives_090922_combin6.csv",index_col=0)


d3['models']=[ ast.literal_eval(k) for k in d3.pair]
import ast #toconvert the string to list
interested_cols=['min_competition', 'max_cooperation', 'min_abiotic sources']
interested_col=interested_cols[2]
d['models']=[ ast.literal_eval(k) for k in d.pair]
d['model1']=[k[0] for k in d.models]
d['model2']=[k[1] for k in d.models]


d_n=pd.DataFrame()
d_n['model1']=frame_combined_all.model.unique()
d_n['model2']=frame_combined_all.model.unique()
d_n[interested_col]=np.nan
heatmap_d=d[['model1','model2',interested_col]]
heatmap_d=heatmap_d.append(d_n)


#to have combined
"from here its to have combined matrix new"
heatmap_1=d[['model1','model2',interested_cols[0]]]
heatmap_1['int_col']=heatmap_1[interested_cols[0]]
heatmap_1['kind']=interested_cols[0]

heatmap_2=d[['model1','model2',interested_cols[1]]]
heatmap_2['int_col']=heatmap_2[interested_cols[1]]

heatmap_2_copy=pd.DataFrame()
heatmap_2_copy['model1']=heatmap_2.model2
heatmap_2_copy['model2']=heatmap_2.model1
heatmap_2_copy['int_col']=heatmap_2.int_col
heatmap_2_copy['kind']=interested_cols[1]


heatmap_d1=heatmap_1[['model1','model2','int_col','kind']]
heatmap_d2=heatmap_2_copy[['model1','model2','int_col','kind']]

heatmap_d=heatmap_d1.append(heatmap_d2)

d_n=pd.DataFrame()
d_n['model1']=frame_combined_all.model.unique()
d_n['model2']=frame_combined_all.model.unique()
d_n["int_col"]=np.nan
d_n["kind"]="none"
heatmap_d=heatmap_d.append(d_n)
heatmap_d=heatmap_d.reset_index()
interested_col='int_col'
df_pivoted = heatmap_d.pivot(columns='model1', index='model2', values=interested_col)


df_pivoted_1 = heatmap_d[heatmap_d.kind.isin(['none',interested_cols[0]])].pivot(columns='model1', index='model2', values=interested_col)
df_pivoted_2 = heatmap_d[heatmap_d.kind.isin(['none',interested_cols[1]])].pivot(columns='model1', index='model2', values=interested_col)

#cm.PuBuGn
fig,ax=plt.subplots()
cax1=ax.imshow(df_pivoted_1,cmap=cm.Reds,vmin= heatmap_d[heatmap_d.kind.isin(['none',interested_cols[0]])][interested_col].min(),vmax= heatmap_d[heatmap_d.kind.isin(['none',interested_cols[0]])][interested_col].max(),alpha=alpha_chosen)
cax2=ax.imshow(df_pivoted_2,cmap=cm.Greens,vmin= heatmap_d[heatmap_d.kind.isin(['none',interested_cols[1]])][interested_col].min(),vmax= heatmap_d[heatmap_d.kind.isin(['none',interested_cols[1]])][interested_col].max(),alpha=alpha_chosen)
ax.set_yticks(np.linspace(0,6,7))#, [str(k) for k in df_pivoted.index.to_list()])
ax.set_xticks(np.linspace(0,6,7))#,[str(k) for k in df_pivoted.columns.to_list()])
ax.set_xticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.columns.to_list()],rotation=90,fontsize=8,style='italic')
ax.set_yticklabels([str(k.split('_')[0]+' '+k.split('_')[1]) for k in df_pivoted.index.to_list()],fontsize=8,style='italic')

#
for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):
        text=(df_pivoted.iloc[y, x])
        if not math.isnan(text):
            text=int(text)
        plt.text(x , y ,  text,
                 horizontalalignment='center',
                 verticalalignment='center',
                fontsize=12 )

ax.set_xticks(np.arange(-.5, 6, 1), minor=True)
ax.set_yticks(np.arange(-.5, 6, 1), minor=True)
ax.grid(which="minor", color="black", linestyle='-', linewidth=0.5)
cbar = fig.colorbar(cax1)#, ticks=[0, 0.5, 1])
cbar2 = fig.colorbar(cax2)#, ticks=[0, 0.5, 1])

plt.tight_layout()

plt.savefig(output_file_figures+'/NAnew_TRY_140323_cmapinterested_col_oted_common_mets_competition_{}.svg'.format(interested_col))




















"---------"


from math import ceil
from holoviews.plotting.util import process_cmap

colormaps = hv.plotting.list_cmaps()
spacing = np.linspace(0, 1, 64)[np.newaxis]
opt_kwargs = dict(aspect=6, xaxis=None, yaxis=None, sublabel_format='')

def filter_cmaps(category):
    return hv.plotting.util.list_cmaps(records=True,category=category,reverse=False)

def cmap_examples(category,cols=4):
    cms = filter_cmaps(category)
    n = len(cms)*1.0
    c=ceil(n/cols) if n>cols else cols
    bars = [hv.Image(spacing, ydensity=1, label="{0} ({1})".format(r.name,r.provider))\
            .opts(cmap=process_cmap(r.name,provider=r.provider), **opt_kwargs)
           for r in cms]
    return hv.Layout(bars).opts(vspace=0.1, hspace=0.1, transpose=(n>cols)).cols(c)



#create the matrix

#kullabergensis
names=['L.kullabergensis','L.apis','L.mellis','L.mellifer','B.asteroides','G.apicola','S.alvi']
matrix = [
    [0, 0.5, 0.37, 0.506, 0.42, 0.35, 0.30],
    [0.5, 0, 0.37, 0.46, 0.34, 0.28, 0.25],
    [0.37, 0.37, 0, 0.41, 0.36, 0.27, 0.35],
    [0.506, 0.46, 0.41, 0, 0.45, 0.33, 0.31],
    [0.42, 0.34, 0.34, 0.45, 0, 0.37, 0.28],
    [0.35, 0.28, 0.27, 0.33, 0.37, 0, 0.27],
    [0.30, 0.25, 0.35, 0.31, 0.28, 0.27, 0],
]


#try with holoviews
