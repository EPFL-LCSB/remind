import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

frames = dict()
frames['Geobacter'] = pd.read_csv('ilp_solutions_Geobacter_Rhodoferax/_data/Geobacter.csv')
frames['Rhodoferax'] = pd.read_csv('ilp_solutions_Geobacter_Rhodoferax/_data/Rhodoferax.csv')

interactions = ['NH4+','SO42-','Acetate','Fe3+','PO43-']

name_dict = {
    'Acetate' : 'Acetate',
    'SO42-': 'Sulfate',
    'PO43-': 'Phosphate',
    'NH4+': 'Ammonia',
    'Fe2+': 'Fe (II)',
    'Fe3+' : 'Fe (III)',
    'H2' : 'H$_{2}$',
    'S' : 'Sulfur',
    'S2-': 'Disulfur',
    'CO2': 'CO$_{2}$',
    'H2O' : 'H$_{2}$O',
    'H2S' : 'H$_{2}$S',
    'Formate': 'Formate',
    'H+'   : 'H$^{+}$',
    'glycolate' : 'Glycolate',
    'O2'  : 'O$_{2}$'
    }

lengths = []
for k,frame_int in frames.items():
    frame_int.columns = ['uptakes', 'secretions', 'yields']
    #positive
    l=[]
    for k in frame_int.uptakes:
        mets = [x.strip() for x in k.split(',')] 
        l+=[m for m in mets if m not in l]
    #negative
    l_=[]
    for k in frame_int.secretions:
        mets = [x.strip() for x in k.split(',')]
        l_+=[m for m in mets if m not in l_]
    
    s_=pd.Series(l).value_counts()
    s__=pd.Series(l_).value_counts()
    columns = [x for x in interactions] # to put the interactions first
    columns += [x for x in set([k for k in s_.index] + [k for k in s__.index]) \
              if x not in  interactions
             ]
        
    lengths.append(len(columns))
    
width_subplots = [lengths[0], 1, lengths[1], 1]


fig,axes = plt.subplots(ncols=4, nrows=1, gridspec_kw={#'height_ratios': height_subplots,
                                                       'width_ratios' : width_subplots},
                        figsize = (15,40))


for ind, (k,frame_int) in enumerate(frames.items()):
    frame_int.columns = ['uptakes', 'secretions', 'yields']
    #
    # plt.rcParams.update({'font.size':16,
    #                       'font.family':'Calibri',
    #                      'legend.fontsize':16,
    #                    'xtick.labelsize' : 14,
    #                    'ytick.labelsize' : 14
    #   })
    
    plt.rcParams['svg.fonttype'] = 'none'
    
    #positive
    l=[]
    for k in frame_int.uptakes:
        mets = [x.strip() for x in k.split(',')] 
        l+=[m for m in mets if m not in l]
    #negative
    l_=[]
    for k in frame_int.secretions:
        mets = [x.strip() for x in k.split(',')]
        l_+=[m for m in mets if m not in l_]
    
    s_=pd.Series(l).value_counts()
    s__=pd.Series(l_).value_counts()
    cols = [x for x in interactions] # to put the interactions first
    cols += [x for x in set([k for k in s_.index] + [k for k in s__.index]) \
              if x not in  interactions
             ]    
    
    dummyarray = np.empty((len(frame_int.uptakes),len(cols)))
    dummyarray[:] = np.nan
    data_coop=pd.DataFrame(dummyarray,columns=cols)
    for k in range(len(frame_int.uptakes)):
        for met in cols:
            if met in [x.strip() for x in frame_int.uptakes.iloc[k].split(',')]:
                if met in interactions:
                    data_coop.iloc[k][met] = 0.5
                else:
                    data_coop.iloc[k][met] = 0
            elif met in [x.strip() for x in frame_int.secretions.iloc[k].split(',')]:
                data_coop.iloc[k][met] = 1.0
    
    
    ###
    
    "plot_df_pivoted"
    df_pivoted=data_coop
    ax = axes[2*ind]
    # figsize=(5,10)
    # fig,ax=plt.subplots(figsize=(5,10))
    heatmap=matplotlib.cm.cividis
    cax1=ax.imshow(df_pivoted,cmap=heatmap,vmin=0,vmax=1,alpha=0.5)
    
    current_cmap = matplotlib.cm.get_cmap(heatmap)
    current_cmap.set_bad(color='white')
    ax.imshow(df_pivoted,cmap=current_cmap,vmin=0,vmax=1,alpha=1.0)
    ax.xaxis.tick_top()
    
    ax.set_xticks(np.linspace(0,df_pivoted.columns.nunique()-1,df_pivoted.columns.nunique()))#, [str(k) for k in df_pivoted.index.to_list()])
    ax.set_yticks([])
    
    labels = [name_dict[str(k.split('_')[0])] for k in df_pivoted.columns.to_list()]
    ax.set_xticklabels(labels,rotation=90,fontsize=30,
                       # style='italic'
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
    
    plt.tight_layout()
    
    #define output
    output_file_figures="define"
    
# plotting the yields
for ind, (k,frame_int) in enumerate(frames.items()):
    frame_int.columns = ['uptakes', 'secretions', 'yields']
    #
    # plt.rcParams.update({'font.size':16,
    #                       'font.family':'Calibri',
    #                      'legend.fontsize':16,
    #                    'xtick.labelsize' : 14,
    #                    'ytick.labelsize' : 14
    #   })
    
    plt.rcParams['svg.fonttype'] = 'none'
        
    
    dummyarray = np.empty((len(frame_int.yields),1))
    dummyarray[:] = np.nan
    data_coop=pd.DataFrame(dummyarray,columns=['Yield'])
    for k in range(len(frame_int.yields)):
        y = float(frame_int.yields.iloc[k].replace('%', ''))/100
        data_coop.iloc[k]['Yield'] = y
    
    
    ###
    
    "plot_df_pivoted"
    df_pivoted=data_coop
    ax = axes[2*ind+1]
    # figsize=(5,10)
    # fig,ax=plt.subplots(figsize=(5,10))
    heatmap=matplotlib.cm.Reds
    cax1=ax.imshow(df_pivoted,cmap=heatmap,vmin=0,vmax=1,alpha=0.5)
    
    current_cmap = matplotlib.cm.get_cmap(heatmap)
    current_cmap.set_bad(color='white')
    ax.imshow(df_pivoted,cmap=current_cmap,vmin=0,vmax=1,alpha=1.0)
    ax.xaxis.tick_top()
    
    ax.set_xticks(np.linspace(0,df_pivoted.columns.nunique()-1,df_pivoted.columns.nunique()))#, [str(k) for k in df_pivoted.index.to_list()])
    ax.set_yticks([])
    
    ax.set_xticklabels(['Yield (of the maximum)'],rotation=90,fontsize=30,
                       # style='italic'
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
    
    plt.tight_layout()
    
    #define output
    output_file_figures="define"
    
    
    
plt.savefig('ilp_solutions_Geobacter_Rhodoferax/_data/data_integ_rhodof.svg')

