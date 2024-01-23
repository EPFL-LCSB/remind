import pandas as pd
import numpy as np
from itertools import combinations


data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'

data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_111023_bee.h5'
frame_combined_all = pd.read_hdf(data_path)
from itertools import combinations


df_bg=pd.read_hdf('/remind/projects/bee_project/060623_ilp_solutions_nplusone_2_member_bee/obj_num_2_modelofint_0_alternatives_interaction_positive_directional_general_yield/sols_active_alternative_for_diff_obj_2_alternation_interaction_positive_directional_general_yield_altno_606_.h5')

df_bs=pd.read_hdf("/remind/projects/bee_project/060623_ilp_solutions_nplusone_2_member_bee/obj_num_2_modelofint_5_alternatives_interaction_positive_directional_general_yield/sols_active_alternative_for_diff_obj_2_alternation_interaction_positive_directional_general_yield_altno_386_.h5")

df_sg=pd.read_hdf("/remind/projects/bee_project/060623_ilp_solutions_nplusone_2_member_bee/obj_num_2_modelofint_10_alternatives_interaction_positive_directional_general_yield/sols_active_alternative_for_diff_obj_2_alternation_interaction_positive_directional_general_yield_altno_530_.h5")

df_bsg=pd.read_hdf("/remind/projects/bee_project/ilp_solutions_070823_nplusone_2_member_bee/obj_num_2_modelofint_6_alternatives_interaction_positive/sols_active_alternative_for_diff_obj_2_alternation_interaction_positive_altno_1283_.h5")


models_dict_bg={"model1": "Bifido",
               "model2": "Gapicola"}
models_dict_bs={"model1": "Bifido",
               "model2": "Salvi"}

models_dict_sg={"model1": "Salvi",
               "model2": "Gapicola"}

models_dict_bsg={"model1": "Bifido",
                "model2": "Salvi",
               "model3": "Gapicola"}


frame_bg=get_species_data(df_bg,models_dict=models_dict_bg,len_models=2)
frame_bs=get_species_data(df_bs,models_dict=models_dict_bs,len_models=2)
frame_sg=get_species_data(df_sg,models_dict=models_dict_sg,len_models=2)
frame_bsg=get_species_data(df_bsg,models_dict=models_dict_bsg,len_models=2)



def give_union(df):
    f_ = []
    for k in df.pos_int.unique():
        f_ += list(k)

    f__ = pd.DataFrame(f_).value_counts()

    cols = [k[0] for k in f__.index]

    return cols

bs_un=[all_exchanges[all_exchanges.mets==k].name_met.iloc[0] for k in give_union(frame_bs)]
bg_un=[all_exchanges[all_exchanges.mets==k].name_met.iloc[0] for k in give_union(frame_bg)]
sg_un=[all_exchanges[all_exchanges.mets==k].name_met.iloc[0] for k in give_union(frame_sg)]
bsg_un=[all_exchanges[all_exchanges.mets==k].name_met.iloc[0] for k in give_union(frame_bsg)]

all_inter=[k for k in list(set(bs_un) & set(bg_un) & set(sg_un))]

bs_bg_inter=[k for k in list(set(bs_un) & set(bg_un))]

bg_sg_inter=[k for k in list(set(sg_un) & set(bg_un))]
bs_sg_inter=[k for k in list(set(sg_un) & set(bs_un))]

all_inter_2= list(set(bs_bg_inter + bg_sg_inter+bs_sg_inter))

bs_diff=[k for k in bs_un if k not in all_inter_2]
bg_diff=[k for k in bg_un if k not in all_inter_2]
gs_diff=[k for k in sg_un if k not in all_inter_2]


import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

# Line style: can be 'dashed' or 'dotted' for example
plt.figure()
v=venn3(subsets = (1, 5, 1, 2,1,2,8), set_labels = ('B+S', 'B+G', 'S+G'))
c=venn3_circles(subsets = (1, 5, 1, 2,1,2,8), linestyle='solid', linewidth=1, color="grey")
plt.savefig(output_file_figures+'/3_venn_comparison_subs_prods_overall.svg')

plt.figure()
v=venn3(subsets = (1,5,1,4,2,1,10), set_labels = ('B+S', 'B+G', 'S+G'))
c=venn3_circles(subsets = (1,5,1,4,2,1,10), linestyle='solid', linewidth=1, color="grey")
plt.savefig(output_file_figures+'/231023_3_venn_comparison_subs_prods_overall.svg')


bs_un=[k for k in give_union(frame_bs)]
bg_un=[k for k in give_union(frame_bg)]
sg_un=[k for k in give_union(frame_sg)]
bsg_un=[k for k in give_union(frame_bsg)]

all_inter=[k for k in list(set(bs_un) & set(bg_un) & set(sg_un))]

bs_bg_inter=[k for k in list(set(bs_un) & set(bg_un))]

bg_sg_inter=[k for k in list(set(sg_un) & set(bg_un))]
bs_sg_inter=[k for k in list(set(sg_un) & set(bs_un))]

all_inter_2= list(set(bs_bg_inter + bg_sg_inter+bs_sg_inter))

bs_diff=[k for k in bs_un if k not in all_inter_2]
bg_diff=[k for k in bg_un if k not in all_inter_2]
gs_diff=[k for k in sg_un if k not in all_inter_2]


all_mets=give_union(frame_bsg)


#second way
frame_bs=frame_combined_all[frame_combined_all.model.isin(['Bifidobacterium_asteroides_PRL2011GF','Snodgrassella_alvi_wkB2GF'])]

frame_bg=frame_combined_all[frame_combined_all.model.isin(['Bifidobacterium_asteroides_PRL2011GF','Gilliamella_apicola_wkB1GF'])]

frame_sg=frame_combined_all[frame_combined_all.model.isin(['Snodgrassella_alvi_wkB2GF','Gilliamella_apicola_wkB1GF'])]


frame=get_subs_prods_per_group(frame_combined_all,gr=['model']).reset_index()


all_pairs=list(combinations(frame_combined_all.model.unique(), 2))
all_pairs=[('Bifidobacterium_asteroides_PRL2011GF', 'Gilliamella_apicola_wkB1GF'),
            ('Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'),
            ('Gilliamella_apicola_wkB1GF', 'Snodgrassella_alvi_wkB2GF'),]
all_exchanges = pd.read_csv("/remind/projects/bee_project/stats/all_exchanges_stats.csv")
all_exchanges = all_exchanges.groupby('mets').first().reset_index()
# met_name = [all_exchanges[all_exchanges.mets == k].name_met.values[0] for k in cols]

coop1_2_list=[]
coop2_1_list=[]
competition_list=[]
coop_union_list=[]
for index_pair in range(len(all_pairs)):
    model1=all_pairs[index_pair][0]
    model2=all_pairs[index_pair][1]
    frame[frame.model.isin([model1])]
    subs1=    frame[frame.model.isin([model1])].subs_list.values[0]
    subs2=    frame[frame.model.isin([model2])].subs_list.values[0]

    prods1=    frame[frame.model.isin([model1])].prods_list.values[0]
    prods2=    frame[frame.model.isin([model2])].prods_list.values[0]

    coop1_2=[all_exchanges[all_exchanges.mets == k.replace('EX_','')].name_met.values[0]  for k in prods1 if k in subs2]
    coop2_1=[all_exchanges[all_exchanges.mets == k.replace('EX_','')].name_met.values[0]  for k in prods2 if k in subs1]
    competition=[all_exchanges[all_exchanges.mets == k.replace('EX_','')].name_met.values[0] for k in subs2 if k in subs1]

    coop1_2=[all_exchanges[all_exchanges.mets == k.replace('EX_','')].mets.values[0].replace('_e','')  for k in prods1 if k in subs2]
    coop2_1=[all_exchanges[all_exchanges.mets == k.replace('EX_','')].mets.values[0].replace('_e','')  for k in prods2 if k in subs1]
    competition=[all_exchanges[all_exchanges.mets == k.replace('EX_','')].mets.values[0].replace('_e','')  for k in subs2 if k in subs1]
    coop_union=list(set(coop1_2+coop2_1))

    # df=pd.DataFrame()
    # df['coop1_2']=coop1_2
    # df['model1'] = model1
    # df['model2'] = model2
    # df['coop2_1']=coop2_1
    # df['competition']=competition

    coop1_2_list.append(tuple(coop1_2))
    coop2_1_list.append(tuple(coop2_1))
    competition_list.append(tuple(competition))
    coop_union_list.append(tuple(coop_union))


model1=[k[0] for k in all_pairs]
model2=[k[1] for k in all_pairs]
data_annot=pd.DataFrame()
data_annot['model1']=model1
data_annot['model2']=model2


data_annot['coop1_2']=coop1_2_list
data_annot['coop2_1']=coop2_1_list
data_annot['competition']=competition_list
data_annot['coop_union'] = coop_union_list


bg=data_annot.coop_union.iloc[0]
bs=data_annot.coop_union.iloc[1]
sg=data_annot.coop_union.iloc[2]

bg_bs=[k for k in bg if k in bs]

bg_sg=[k for k in bg if k in sg]

bs_sg=[k for k in bs if k in sg]

all_int=[k for k in bg_bs if k in sg]

all_inter_2= list(set(bg_bs + bg_sg+bs_sg))

bs_diff=[k for k in bs if k not in all_inter_2]
bg_diff=[k for k in bg if k not in all_inter_2]
gs_diff=[k for k in sg if k not in all_inter_2]
