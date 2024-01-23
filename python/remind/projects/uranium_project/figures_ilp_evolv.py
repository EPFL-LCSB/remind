import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

rep = ' // '

sheet_to_df_map = pd.read_excel('ilp_solutions_Geobacter_Rhodoferax/_physiology/physiology.xlsx', sheet_name=None,
                                index_col=0, usecols=[0,1,2,3,4,5])
inorganics = ['fe2_e', 'fe3_e', 's_e', 'h2_e', 'co2_e', 'ss_e',
              'h_e', 'h2o_e']
nitro = ['nh4_e', 'n2_e']

name_dict = {
    'ac' : 'Acetate',
    'cit': 'Citrate',
    'fru': 'Fructose',
    'fum': 'Fumarate',
    'lac-L': 'Lactate',
    'mal-L': 'Malate',
    'ppa' : 'Propionate',
    'pyr' : 'Pyruvate',
    'n2' : 'N$_{2}$',
    'nh4': 'Ammonia',
    'co2': 'CO$_{2}$',
    'fe2': 'Fe (II)',
    'fe3' : 'Fe (III)',
    'h2' : 'H$_{2}$',
    's' : 'Sulfur',
    'ss': 'Disulfur',
    }

#positive
l=[]
height_subplots = []
for key, frame_int in sheet_to_df_map.items():
    new_frame = frame_int.drop_duplicates()
    height_subplots += [len(new_frame)]
    
    for col in frame_int.columns:
        for i,k in frame_int[col].iteritems():
            try: 
                mets = [x.strip() for x in k.split(rep)] 
                l+=[m for m in mets if m not in l]
            except AttributeError: # it is empty or it is the order
                frame_int[col][i] = ''
    #negative
    # l_=[]
    # for k in frame_int.secretions:
    #     mets = [x.strip() for x in k.split(',')]
    #     l_+=[m for m in mets if m not in l_]
    
s_=pd.Series(l).value_counts()
# s__=pd.Series(l_).value_counts()
# cols = [x for x in interactions] # to put the interactions first
cols = {}
found = [] #  a list to avoid duplicates
c_source = [x for x in set([k for k in s_.index]) if x not in inorganics and x not in nitro] # C sources 
cols['c source'] =  sorted(c_source)
found += c_source
n_sources = [x for x in set([k for k in s_.index]) if x not in inorganics and x not in found] # N sources
cols['n source'] = sorted(n_sources)
found += n_sources
the_rest = [x for x in set([k for k in s_.index]) if x not in found]
cols['the rest'] = sorted(the_rest)

width_subplots = [len(c_source), len(n_sources), len(the_rest)]

# frame_int = sheet_to_df_map['min_tot_upt']
fig,axes = plt.subplots(ncols=3, nrows=len(sheet_to_df_map),gridspec_kw={'height_ratios': height_subplots,
                                                                'width_ratios' : width_subplots},
                        figsize = (30,100))
for index, (key, frame_int) in enumerate(sheet_to_df_map.items()):

    frame_int.columns = ['Comp', 'Geo', 'Rhod', 'G__R', 'R__G']
    #
    # plt.rcParams.update({'font.size':16,
    #                       'font.family':'Calibri',
    #                      'legend.fontsize':16,
    #                    'xtick.labelsize' : 14,
    #                    'ytick.labelsize' : 14
    #   })
    frame_int = frame_int.drop_duplicates()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.subplots_adjust(hspace=0.2, wspace=0.05)
    
    for ind , (k, col) in enumerate(cols.items()):
        dummyarray = np.empty((len(frame_int.Comp),len(col)))
        dummyarray[:] = np.nan
        data_coop=pd.DataFrame(dummyarray,columns=col)
        for k in range(len(frame_int.Comp)):
            for met in col:
                if met in [x.strip() for x in frame_int.Comp.iloc[k].split(rep)]:
                    data_coop.iloc[k][met] = 1
                if met in [x.strip() for x in frame_int.Geo.iloc[k].split(rep)]:
                    data_coop.iloc[k][met] = 0
                if met in [x.strip() for x in frame_int.Rhod.iloc[k].split(rep)]:
                    data_coop.iloc[k][met] = 0.4
                if met in [x.strip() for x in frame_int.G__R.iloc[k].split(rep)]:
                    data_coop.iloc[k][met] = 0.7
                if met in [x.strip() for x in frame_int.R__G.iloc[k].split(rep)]:
                    data_coop.iloc[k][met] = 0.2
        
        
        ###
        
        "plot_df_pivoted"
        df_pivoted=data_coop
        # figsize=(5,10)
        ax = axes[index, ind]
        heatmap=matplotlib.cm.Dark2
        cax1=ax.imshow(df_pivoted,cmap=heatmap,vmin=0,vmax=1,alpha=0.5)
        
        current_cmap = matplotlib.cm.get_cmap(heatmap)
        current_cmap.set_bad(color='white')
        ax.imshow(df_pivoted,cmap=current_cmap,vmin=0,vmax=1,alpha=1.0)
        
        
        ax.set_yticks([])
        ax.xaxis.tick_top()
        
        labels = [name_dict[str(k.split('_')[0])] for k in df_pivoted.columns.to_list()]
        ax.set_xticks(np.linspace(0,df_pivoted.columns.nunique()-1,df_pivoted.columns.nunique()))#, [str(k) for k in df_pivoted.index.to_list()])
        if index ==0:
            ax.set_xticklabels(labels,rotation=90,fontsize=40 #,style='italic'
                               )
        else:
            ax.set_xticklabels([],rotation=90,fontsize=30 #,style='italic'
                               )    
            
        # ax.set_yticklabels([k+1 for k in df_pivoted.index.to_list()],rotation=0,fontsize=16)
        
        #in case for annotation
        for y in range(df_pivoted.shape[0]):
            for x in range(df_pivoted.shape[1]):
                text=(df_pivoted.iloc[y, x])
                if not math.isnan(text):
                    text=str(np.round(float(text),1))
                    # plt.text(x , y ,  text,
                    #      horizontalalignment='center',
                    #      verticalalignment='center',
                    #     fontsize=3 )
        
        
        # ax.xaxis.set_minor_locator(AutoMinorLocator())
        # ax.yaxis.set_minor_locator(AutoMinorLocator())
        
        #
        ax.set_xticks(np.arange(-.5, df_pivoted.columns.nunique()-1, 1), minor=True)
        ax.set_yticks(np.arange(-.5, df_pivoted.index.nunique()-1, 1), minor=True)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=1.0)
        # cbar = fig.colorbar(cax1)#, ticks=[0, 0.5, 1])
        
        # plt.tight_layout()
        
        #define output
        output_file_figures="define"
plt.savefig('ilp_solutions_Geobacter_Rhodoferax/_physiology/physiology.svg')

