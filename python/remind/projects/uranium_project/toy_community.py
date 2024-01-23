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
from remind.optim.variables import UptakeActivation, SecretionActivation

from remind.io.json import save_json_model


data_path = "ime_alternatives/imes_toy/DIME_{}_alternative_mets.csv"
model_path = "models/{}.json" # ThermoModel objects

# defining the setup
setup = {
         'sense'           : 'any', # the sense of the objective for the interaction
         'diff_interaction': True, # to differentiate between cooperation and competition
         }

# defining the members of the community
species_data = {
                'Geobacter':
                {'model_name'  : 'Geobacter',
                  # 'alternatives': 1000,
                  'growth'      : 0.2},
                'Rhodoferax':
                {'model_name'  : 'Rhodoferax',
                  # 'alternatives': 1000,
                  'growth'      : 0.2},
                # 'Shewanella':
                # {'model_name'  : 'Shewanella',
                #   # 'alternatives': 326,
                #   'growth'      : 0.2},
                }

medium = pd.read_csv('ime_alternatives/imes_toy/medium.csv',index_col=0).to_dict(orient='index')
met_type = {'_{}'.format(k):v['type'] for k,v in medium.items()}
setup['met_type'] = met_type
    
# Take the models and DIMEs
model_dict = dict()
dime_dict = dict()
for species, data in species_data.items():
    model_dict[species] = load_json_model(
        model_path.format(data['model_name']))
    dime_dict[species] = pd.read_csv(
        data_path.format(data['model_name'], #data['growth']
                         ))

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

# add the species to the community
model.add_species(species_dict.values(), **setup)

# setting the objective function
num_rxns = symbol_sum(model.get_variables_of_type(UptakeActivation)) + \
    symbol_sum(model.get_variables_of_type(SecretionActivation))
           
####
model.objective = num_rxns 
model.objective_direction = 'min'


save_json_model(model, 'community_models/Geobacter_Rhodoferax_community.json')