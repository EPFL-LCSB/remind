#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: ReMIND Team
"""

from remind.io.json import load_json_model
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse,AbioticResource
from remind.optim.alternatives import find_alternative_solutions

from pytfa.optim.utils import symbol_sum
from remind.core.medium import constrain_abiotics
#
import os
from sys import argv

#normally for optimal
# criteria_lim=0
#until the point they are connected



# criteria_lim=0


# model = load_json_model('../projects/community_models/{}_community.json'.format(comm_name))

alternation_opts=['reaction','interaction','interaction_positive','interaction_negative','abiotic_sources',\
                  'reaction_yield','interaction_yield','interaction_positive_yield','interaction_negative_yield','abiotic_sources_yield',"abiotic_sources_uptake",
                  'interaction_positive_directional_general']

_, objective_no,alt,crit=argv
obj_num=int(objective_no)
alternation=alternation_opts[int(alt)]
criteria_lim=int(crit)


comm_name = '7_member_bee'
max_alternative=5000

# obj_num=6
# alternation='abiotic_sources'
# filepath = ("/remind/projects/bee_project/ilp_model/model_7_member_190722_trial.json")
# filepath=("/remind/projects/bee_project/ilp_model/model_7_member_050822_correct.json")
# filepath=("/remind/projects/bee_project/ilp_model/model_7_member_120822.json")
# filepath=("/remind/projects/bee_project/ilp_model/model_7_member_070922.json")

#before
# filepath=("/remind/projects/bee_project/ilp_model/model_7_member_051022.json")
#after
filepath=("/remind/projects/bee_project/ilp_model/model_7_member_120623.json")

filepath=("/remind/projects/bee_project/ilp_model/model_7_member_111023.json")

model= load_json_model(filepath)
### setting different objective functions for the community
# the number of competitions and cooperations in the community
pos_int = model.get_variables_of_type(PositiveInteraction)
neg_int = model.get_variables_of_type(NegativeInteraction)
competitions = symbol_sum([v for v in neg_int])
cooperations = symbol_sum([v for v in pos_int])

# total number of reactions
num_upt = symbol_sum(model.get_variables_of_type(UptakeActivation))
num_sec = symbol_sum(model.get_variables_of_type(SecretionActivation))


### do the simulation based on the objective that you choose

# # First, try to find the minimal number of active uptakes for the community (no interactiion)
if obj_num==1:
    model.objective = competitions
    model.objective_direction = 'min'
    # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)

if obj_num==2:
    model.objective = cooperations
    model.objective_direction = 'max'
    # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)


# if obj_num==3:
#     model.objective = competitions
#     model.objective_direction = 'max'
#     # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)
#
# if obj_num==4:
#     model.objective = cooperations
#     model.objective_direction = 'min'
#     # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)

#new objective
# if obj_num==5:
#     model.objective = 5*cooperations-competitions
#     model.objective_direction = 'max'


if obj_num==6:
    constrain_abiotics(model)

# output file to save the results
# data_path = '/remind/projects/bee_project/ilp_solutions_051022_correct_{}/obj_num_{}_alternatives_{}'.format(comm_name, obj_num,

data_path = '/remind/projects/bee_project/121023_ilp_solutions_290623_{}/obj_num_{}_alternatives_{}'.format(comm_name, obj_num, alternation)
if not os.path.exists(data_path):
    os.makedirs(data_path)


### solving and finding alternatives
# checking for feasibility and generating alternative solutions
#shouldnt be necessary
growth_ids = {
    'Gapicola': 'Growth',
    'Salvi': 'Growth',
    'Bifido': 'Growth',
    'Lapis': 'Growth',
    'Lkulla': 'Growth',
    'Lmellifer': 'Growth',
    'Lmellis': 'Growth',
}


growth_limits = {
    'Gapicola': 0.2,
    'Salvi': 0.2,
    'Bifido': 0.2,
    'Lapis': 0.2,
    'Lkulla': 0.2,
    'Lmellifer': 0.2,
    'Lmellis': 0.2,
}



sols_all,sols_active = find_alternative_solutions(model, growth_ids, growth_limits,
                                  alternation=alternation, max_alts=max_alternative,check_feasibility = False,stop_cond=True,criteria=criteria_lim)




altsno=sols_active.alternative.nunique()

sols_all.to_hdf(data_path+'/sols_all_alternative_for_diff_obj_{}_alternation_{}_altno_{}_.h5'.format(obj_num,alternation,altsno),key='s')


sols_active.to_hdf(data_path+'/sols_active_alternative_for_diff_obj_{}_alternation_{}_altno_{}_.h5'.format(obj_num,alternation,altsno),key='s')


