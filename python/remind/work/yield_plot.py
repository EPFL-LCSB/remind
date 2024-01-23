import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import glob
import numpy as np

#"read your data"
data=pd.read_csv('min_req/min_req_toy.csv',index_col=0)


# data=pd.read_csv('/remind/projects/bee_project/figures/data_yield/min_req_bee_2_species_170622.csv',index_col=0)
#put your output
output_file_figures= 'min_req/'

plt.rcParams.update({'font.size':16,
                      'font.family':'Calibri',
                     'legend.fontsize':16,
                   'xtick.labelsize' : 14,
                   'ytick.labelsize' : 14
  })

plt.rcParams['svg.fonttype'] = 'none'
#put first column name that goes on the x axis
column1='Geobacter Yield' #name of your first column as in the data frame
#put second column name
value_col='Minimal Abiotic Nutrients' #name of your value column as in the data frame
column2='Rhodoferax Yield' #name of your first column as in the data frame
xlabel_text="G. sulfurreducens" #user defined put what name u want it there to appear
ylabel_text="R. ferrireducens"#user defined put what name u want it there to appear

unique_yields=data[column1].nunique() #your unique yield cuts
fig,ax=plt.subplots()
# plt.figure()
# df = pd.read_csv(StringIO(data_str), names=['x', 'y', 'z'], delim_whitespace=True)
df_pivoted = data.pivot(columns=column1, index=column2, values=value_col)
#here define your colormap and choose in this one it is coolwarm
cax=ax.imshow(df_pivoted,cmap=cm.PuRd,vmin=data[value_col].min(),vmax=data[value_col].max(),alpha=0.8)


ax.set_yticks(np.linspace(0,unique_yields-1,unique_yields))#, [str(k) for k in df_pivoted.index.to_list()])
ax.set_xticks(np.linspace(0,unique_yields-1,unique_yields))#,[str(k) for k in df_pivoted.columns.to_list()])
ax.set_xticklabels([str(int(round(k,1)*100))+' %' for k in df_pivoted.columns.to_list()],rotation=90)
ax.set_yticklabels([str(int(round(k,1)*100))+' %' for k in df_pivoted.index.to_list()])
# plt.xlabel(column1.split('_Yield')[0],style='italic')
# plt.ylabel(column2.split('_Yield')[0],style='italic')
plt.xlabel(xlabel_text,style='italic')
plt.ylabel(ylabel_text,style='italic')
ax.set_xticks(np.arange(-.5, unique_yields-1, 1), minor=True)
ax.set_yticks(np.arange(-.5, 4, 1), minor=True)

#     plt.annotate(label, # this is the text
#                  (df_st['A'].iloc[i],df_st['B'].iloc[i]), # this is the point to label
#                  textcoords="offset points", # how to position the text
#                  xytext=(0,10), # distance from text to points (x,y)
#                  ha='center') # horizontal alignment can be left, right or center

for y in range(df_pivoted.shape[0]):
    for x in range(df_pivoted.shape[1]):

        plt.text(x , y ,  round(df_pivoted.iloc[y, x]),
                 horizontalalignment='center',
                 verticalalignment='center',
                 fontsize=12)

# plt.yticks(np.linspace(0,4,5),data["Geobacter_Yield"]) # sets yticks
# plt.xticks(np.linspace(0,4,5),data["Rhodoferax_Yield"]) # sets yticks
ax.grid(which="minor", color="black", linestyle='-', linewidth=0.5)
cbar = fig.colorbar(cax, ticks=[x for x in range(round(data[value_col].min()),round(data[value_col].max())+1,1)])
# cbar.ax.set_yticklabels([str(0), '0.5', str(1)])
plt.gca().invert_yaxis()
plt.tight_layout()
#todo put desired name and output_file
plt.savefig(output_file_figures+'/lollipop_community.svg')
# plt.savefig(output_file_figures+'/toy_system_heatmap.png')

#
# ax = sns.heatmap(data=df_pivoted, annot=True, fmt='d', cmap='RdYlGn', cbar=True, cbar_kws={'label': 'z'}, square=True)
# ax.tick_params(labelrotation=0)
# plt.show()


