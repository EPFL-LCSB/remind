#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a community model
"""
import pandas as pd
from pytfa.io.json import load_json_model
from pytfa.optim.utils import symbol_sum

from remind.core.interaction import InteractionModel
from remind.core.objects import CommunitySpecies
from remind.utils.postprocessing import get_dime_coupling
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, DIMEVariable
from remind.optim.alternatives import find_alternative_solutions
from remind.utils.merge_files import merge_files_for_sliced_yield
from cobra.io import read_sbml_model
import cobra
import os
from os.path import join
import cobra.test
from cobra.io import read_sbml_model
import glob
import time
import cobra
from remind.io.json import save_json_model, load_json_model
from pytfa.io.json import load_json_model_tmodel_struct, save_json_model_tmodel_struct

from pytfa.io.json import load_json_model as tfa_load_json_model
from pytfa.io.json import save_json_model as tfa_savejson_model
# data_path = "/remind/projects/bee_project/stats/combined_dimes_unique_180522_bee.h5"

# data_path='/remind/projects/bee_project/stats/combined_dimes_unique_050822_bee.h5'
# data_path='/remind/projects/bee_project/stats/combined_dimes_unique_080822_bee.h5'

# data_path='/remind/projects/bee_project/stats/combined_dimes_unique_120822_bee_new.h5'
# data_path='/remind/projects/bee_project/stats/combined_dimes_unique_120822_bee_new.h5'
# data_path='/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'
data_path='/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'


data_path='/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'

# data_path='/remind/projects/bee_project/stats/combined_dimes_unique_081023_bee.h5'


data_path='/remind/projects/bee_project/stats/combined_dimes_unique_111023_bee.h5'
frame_combined_all=pd.read_hdf(data_path)


# data_path='/remind/projects/bee_project/stats/combined_dimes_not_unique_111023_bee.h5'
# frame_combined_nu=pd.read_hdf(data_path)


models_of_int=['Gilliamella_apicola_wkB1GF','Snodgrassella_alvi_wkB2GF']
frame_ilp=frame_combined_all[frame_combined_all.model.isin(models_of_int)]
# alternatives_of_int=[1,2,3,4,5]
# frame_ilp=frame_ilp[frame_ilp.alternative.isin(alternatives_of_int)]
frame_ilp['yield'] = frame_ilp['yield_perc']
frame_ilp['element'] = ['ONE'] * frame_ilp.shape[0]


# model_path='/remind/models/bee_models/core_members/tfa_structures_corrected/{}.json'
model_path='/remind/models/bee_models/core_members/tfa_real_010922/{}.json'
model_path='/remind/models/bee_models/core_members/tfa_latest_111023/{}.json'

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
    # model_dict[species] = load_json_model_tmodel_struct(
    #     model_path.format(data['model_name_json']))

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
filepath = ("/remind/projects/bee_project/ilp_model/model_2_member_111023.json")

# filepath = ("/remind/projects/bee_project/ilp_model/model_2_member_070922.json")
save_json_model(model, filepath)


# model_ilp = load_json_model(filepath)
