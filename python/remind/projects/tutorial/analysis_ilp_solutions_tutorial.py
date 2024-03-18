"""
@author: ReMIND Team
"""
from remind.utils.postprocessing import *
import os
import glob
import pandas as pd

comm_name = '2_member_bee'

data_path="/remind/projects/tutorial/bee_tutorial_ilp_solutions_2_member_bee"


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

if "modelofint" in df.columns:
    gb = df.groupby(['modelofint'])
    keys = [k for k in gb.groups.keys()]
    keys_int = [k for k in keys]
    grouping = ['alternation', 'modelofint','obj_num']
    gb = df.groupby(grouping)
    keys = [k for k in gb.groups.keys()]
    keys_int = [k for k in keys]

#here all data is merged and from keys_int  you can check for which one to analyse

d_int_coop=gb.get_group(keys_int[0])
d_int_coop_direct=gb.get_group(keys_int[1])
d_int_coop_direct_yield=gb.get_group(keys_int[2])


d_int_abiot=gb.get_group(keys_int[3])
d_int_abiot_yield=gb.get_group(keys_int[4])


d_int_comp=gb.get_group(keys_int[5])
d_int_comp_direct=gb.get_group(keys_int[6])
d_int_comp_direct_yield=gb.get_group(keys_int[7])


models_dict={
             "model1":"Gapicola",
             "model2":"Salvi"}

#to get the analysis to reorder columns
frame_int_coop=get_species_data(d_int_coop,models_dict=models_dict,len_models=2)
frame_int_coop_direct=get_species_data(d_int_coop_direct,models_dict=models_dict,len_models=2)
frame_int_coop_direct_yield=get_species_data(d_int_coop_direct_yield,models_dict=models_dict,len_models=2)


frame_int_abiot=get_species_data(d_int_abiot,models_dict=models_dict,len_models=2)
frame_int_abiot_yield=get_species_data(d_int_abiot_yield,models_dict=models_dict,len_models=2)


frame_int_comp=get_species_data(d_int_comp,models_dict=models_dict,len_models=2)
frame_int_comp_direct=get_species_data(d_int_comp_direct,models_dict=models_dict,len_models=2)
frame_int_comp_direct_yield=get_species_data(d_int_comp_direct_yield,models_dict=models_dict,len_models=2)


#columns description
#from here on you can check the dataframes created and check for patterns with pos_int
#pos_int==> positive_interaction pattern #cross-feeding
#neg_int==> negative interaction pattern #competition
#yield ==> yield of the indicated model supporting this interaction
#secretions==> secretion of the indiciated model
#uptake ==> uptakes of the indicated model
#dimes ==> dimes of the models
#upt_counts==number of species that uptake this metabolite
#secretion_counts==> number of species that secrete this metabollite