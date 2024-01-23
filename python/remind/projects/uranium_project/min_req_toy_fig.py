#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from remind.io.json import load_json_model
import pandas as pd

from pytfa.optim.utils import symbol_sum
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse
    

model = load_json_model('../projects/community_models/Geobacter_Rhodoferax_community.json')


pos_int = model.get_variables_of_type(PositiveInteraction)
neg_int = model.get_variables_of_type(NegativeInteraction)
# exchanges = model.get_variables_of_type(ExchangeActivation)
neg_intcs = symbol_sum([v for v in neg_int]) 
pos_intcs = symbol_sum([v for v in pos_int])
num_upts = symbol_sum(model.get_variables_of_type(UptakeActivation))
num_scrs = symbol_sum(model.get_variables_of_type(SecretionActivation))
    
model.objective = num_upts - pos_intcs - neg_intcs
model.objective_direction = 'min'


yield_vars = model.get_variables_of_type(YieldUse)
yield_labels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
geo_labels  = ['Geobacter_{}'.format(x) for x in yield_labels]
rhod_labels = ['Rhodoferax_{}'.format(x) for x in yield_labels]

sol_dict = {}
ind = 0
for geo_id in geo_labels:
    geo_yield = [var for var in yield_vars if geo_id in var.name][0]
    geo_yield.variable.lb = 1
    for rhod_id in rhod_labels:
        rhod_yield = [var for var in yield_vars if rhod_id in var.name][0]
        rhod_yield.variable.lb = 1
        model.slim_optimize()
        min_req = model.objective.value
        sol_dict[ind] = {'Geobacter Yield' : geo_id.replace('Geobacter_',''),
                         'Rhodoferax Yield': rhod_id.replace('Rhodoferax_',''),
                         'Minimal Abiotic Nutrients' : min_req}
        ind += 1
        rhod_yield.variable.lb = 0
    geo_yield.variable.lb = 0
        
sol_data = pd.DataFrame.from_dict(sol_dict, orient = 'index')
sol_data.to_csv('min_req_toy.csv')
