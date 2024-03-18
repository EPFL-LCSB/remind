#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a community model
"""
import pandas as pd
from remind.core.interaction import InteractionModel
from remind.core.objects import CommunitySpecies
from remind.utils.postprocessing import *
import os
import glob
import time
import cobra
from remind.io.json import save_json_model, load_json_model
from pytfa.io.json import load_json_model as tfa_load_json_model


"here merge your DiMEs into a combined dataframe"

path="/remind/projects/tutorial/tutorial_bee_DiMEs_True_limited_carbon"
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
frame_combined_all=pd.concat(list_removed, ignore_index=True) #here removed the repeating solutions across different
# yield cuts

output_file_all='/remind/projects/tutorial/dimes_stats'
if not os.path.exists(output_file_all):
    os.makedirs(output_file_all)
frame_combined_all.to_hdf(output_file_all+'/combined_dimes_unique_tutorial_bee.h5',key='s')

"""prepare for model building"""
models_of_int=['Gilliamella_apicola_wkB1GF','Snodgrassella_alvi_wkB2GF']
frame_ilp=frame_combined_all[frame_combined_all.model.isin(models_of_int)]
frame_ilp['yield'] = frame_ilp['yield_perc']
frame_ilp['element'] = ['ONE'] * frame_ilp.shape[0]

model_path='/remind/models/bee_models/core_members/tfa_real_101023/{}.json'

# defining the setup
setup = {
    'sense': 'any',  # the sense of the objective for the interaction
    'diff_interaction': True,  # to differentiate between cooperation and competition
}
species_data = {
    'Gapicola':
        {'model_name_json': 'tmodel_Gilliamella_apicola_wkB1GF',
         'model_name':'Gilliamella_apicola_wkB1GF',
         'growth': 0.2},
    'Salvi':
        {'model_name_json': 'tmodel_Snodgrassella_alvi_wkB2GF',
          'model_name': 'Snodgrassella_alvi_wkB2GF',
         'growth': 0.2}

}
# Take the models and DIMEs
model_dict = dict()
dime_dict = dict()
for species, data in species_data.items():

    model_dict[species] = tfa_load_json_model(
        model_path.format(data['model_name_json']))
    dime_dict[species] = frame_ilp[frame_ilp.model.isin([data['model_name']])]

# Covert the models to Species model
species_dict = dict()
for species, mod in model_dict.items():
    # the name of the model
    mod.name = species
    mod.id = species
    # take the coupling dict
    coupling_dict = get_dime_coupling(dime_dict[species])
    # convert and add the DIMEs
    new = CommunitySpecies(mod, inplace=False)
    new.add_dimes(coupling_dict)
    species_dict[species] = new

# initiate the model
model = InteractionModel()

t = time.time()
'generate untill max alternative'

# add the species to the community
model.add_species(species_dict.values(), **setup)
elapsed = time.time() - t
print('time for adding species is ', elapsed)
#### save the model
filepath = ("/remind/projects/tutorial/ilp_model_tutorial/model_2_member_180324.json")
save_json_model(model, filepath)

