#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 10:05:01 2022

@author: ReMIND Team
"""
import os
from remind.io.json import load_json_model
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse
from remind.optim.alternatives import find_alternative_solutions
    
from pytfa.optim.utils import symbol_sum
from remind.core.medium import constrain_abiotics


from sys import argv


comm_name = 'Geobacter_Rhodoferax'
model = load_json_model('../../projects/community_models/{}_community.json'.format(comm_name))

### setting different objective functions for the community
# the number of competitions and cooperations in the community
pos_int = model.get_variables_of_type(PositiveInteraction)
neg_int = model.get_variables_of_type(NegativeInteraction)
competitions = symbol_sum([v for v in neg_int]) 
cooperations = symbol_sum([v for v in pos_int])

# total number of reactions
num_upt = symbol_sum(model.get_variables_of_type(UptakeActivation))
num_sec = symbol_sum(model.get_variables_of_type(SecretionActivation))
        
# take out solutions that support higher or lower yields
yield_vars = model.get_variables_of_type(YieldUse)
yield_cuts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
organisms = ['Geobacter', 'Rhodoferax']
sum_dict = dict()
for org in organisms:
    sum_dict[org] = list()
    for y in yield_cuts:
        var = [var for var in yield_vars if var.species == org \
               and str(y) in var.name][0]
        sum_dict[org] += [y*var]
yield_expr_geo = symbol_sum(sum_dict['Geobacter'])  
yield_expr_rhod = symbol_sum(sum_dict['Rhodoferax'])

alternation_opts=['reaction','interaction', 'interaction_directional']
# alternation = 'reaction'
# alternation = 'interaction'

_, objective_no,alt=(1, 9, 2)
obj_num=int(objective_no)
alternation=alternation_opts[int(alt)]

### do the simulation based on the objective that you choose

## Physiology
# # First, try to find the minimal number of active uptakes for the community (no interactiion)
if obj_num == 1:
    model.objective = num_upt
    model.objective_direction = 'min'
    mod = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)


# # Second, try to minimize the total abiotic uptake by the community (uptakes - interactions)
if obj_num == 2:    
    model.objective = num_upt - competitions - cooperations
    model.objective_direction = 'min'
    mod = 'ilp_solutions/{}_min_abiotic_upt_{}.csv'.format(comm_name,alternation)


# # Third, try to minimize the total abiotic uptake by the community (uptakes - cooperations)
if obj_num == 3:
    model.objective = num_upt - cooperations
    model.objective_direction = 'min' 
    mod = 'ilp_solutions/{}_min_upt_compt_{}.csv'.format(comm_name,alternation)

# # Forth, try to maximize the cooperations and yield of both organisms
if obj_num == 4:
    model.objective = yield_expr_geo + yield_expr_rhod + cooperations
    model.objective_direction = 'max' 
    mod = 'ilp_solutions/{}_max_yield_coop_{}.csv'.format(comm_name,alternation)

# # Fifth, try to minimize the competitions and maximize yield of both organisms
if obj_num == 5:
    model.objective = yield_expr_geo + yield_expr_rhod - competitions
    model.objective_direction = 'max' 
    mod = 'ilp_solutions/{}_max_yield_min_compt_{}.csv'.format(comm_name,alternation)

# # Sixth, try to minimize the uptakes and maximize yield of both organisms
if obj_num == 6:
    model.objective = 10*yield_expr_geo + 10*yield_expr_rhod - num_upt
    model.objective_direction = 'max' 
    mod = 'ilp_solutions/{}_max_yield_min_upt_{}.csv'.format(comm_name,alternation)


## Engineering
# # First, try to maximize the yield of Geobacter and minimize Rhodoferax
if obj_num == 7:
    model.objective = yield_expr_geo - yield_expr_rhod 
    model.objective_direction = 'max' 
    mod = 'ilp_solutions/{}_max_geo_min_rhod_{}.csv'.format(comm_name,alternation)
    
# # Second, try to maximize the yield of Geobacter and minimize Rhodoferax and minimize uptakes
if obj_num == 8:
    model.objective = 10*yield_expr_geo - 10*yield_expr_rhod - num_upt
    model.objective_direction = 'max' 
    mod = 'ilp_solutions/{}_max_geo_min_rhod_min_upt_{}.csv'.format(comm_name,alternation)
    
# # Third, try to maximize the benefit of Geobacter from Rhodoferax
if obj_num == 9:
    model.objective = cooperations - yield_expr_rhod + 5*yield_expr_geo
    model.objective_direction = 'max' 
    mod = 'ilp_solutions/{}_max_benefit_{}.csv'.format(comm_name,alternation)


## Data Integration
# # the observed interactions are competition for Fe(III), ammonium and acetate
if obj_num == 10:
    comp_ids = ['_ac_e', '_fe3_e', '_nh4_e'] # these interactions should be one and the rest should be zero
    expr = symbol_sum([1-2*neg_int.get_by_id(x) for x in comp_ids]) # multiply by -2 as it appears once in competitions, so the final coefficient will be -1
    model.objective = expr + competitions + cooperations
    model.objective_direction = 'min' 
    mod = 'ilp_solutions/{}_data_integration_{}.csv'.format(comm_name,alternation)
    
# # Eleventh, try to minimize the abiotic uptakes and maximize yield of both organisms
if obj_num == 11:
    constrain_abiotics(model)
    abiot_expr = model.objective.expression
    model.objective = 10*yield_expr_geo + 10*yield_expr_rhod - abiot_expr
    model.objective_direction = 'max' 
    mod = 'ilp_solutions/{}_max_yield_min_abiot_upt_{}.csv'.format(comm_name,alternation)


# output file to save the results
data_path = mod
# if not os.path.exists(data_path):
#     os.makedirs(data_path)

### solving and finding alternatives
# checking for feasibility and generating alternative solutions
growth_ids = {
    'Geobacter' : 'agg_GS13m',
    'Rhodoferax': 'BIO_Rfer3',
    # 'Shewanella': 'Biomass',
    }
growth_limits = {
    'Geobacter' : 0.2,
    'Rhodoferax': 0.2,
    # 'Shewanella': 0.2,
    }

sols,_ = find_alternative_solutions(model, growth_ids, growth_limits, 
                                alternation = alternation, max_alts=50, criteria=0)

sols.to_csv(data_path)

 


