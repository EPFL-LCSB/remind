#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: ReMIND Team
"""

from remind.io.json import load_json_model
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse
from remind.optim.alternatives import find_alternative_solutions

from pytfa.optim.utils import symbol_sum
#
import os
from sys import argv
from remind.core.medium import constrain_abiotics

# stop_crit=2
#default it is 0 if only wanted to stay in optimal solutions

comm_name = '2_member_bee'
# model = load_json_model('../projects/community_models/{}_community.json'.format(comm_name))

# alternation_opts=['reaction','interaction']
alternation_opts=['reaction','interaction','interaction_positive','interaction_negative','abiotic_sources',\
                  'reaction_yield','interaction_yield','interaction_positive_yield','interaction_negative_yield','abiotic_sources_yield',
                  'interaction_positive_directional_2','interaction_positive_directional_general','interaction_positive_directional_general_yield',
                  'interaction_negative_directional_general','interaction_negative_directional_general_yield']

max_alternative=5000

_, objective_no,alt=argv
obj_num=int(objective_no)
alternation=alternation_opts[int(alt)]

# filepath = ("/remind/projects/bee_project/ilp_model/model_2_member_170622_correct.json")

# filepath = ("/remind/projects/bee_project/model_2_member_060822_correct.json")
# filepath = ("/remind/projects/bee_project/ilp_model/model_2_member_060822_correct.json")
# filepath = ("/remind/projects/bee_project/ilp_model/model_2_member_120822.json")
filepath = ("/remind/projects/bee_project/ilp_model/model_2_member_070922.json")
filepath = ("/remind/projects/bee_project/ilp_model/model_2_member_081023.json")
filepath = ("/remind/projects/bee_project/ilp_model/model_2_member_111023.json")

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


# take out solutions that support higher or lower yields
# yield_vars = model.get_variables_of_type(YieldUse)
# yield_cuts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# organisms = ['Geobacter', 'Rhodoferax']
# sum_dict = dict()
# for org in organisms:
#     sum_dict[org] = list()
#     for y in yield_cuts:
#         var = [var for var in yield_vars if var.species == org \
#                and str(y) in var.name][0]
#         sum_dict[org] += [y * var]
# yield_expr_geo = symbol_sum(sum_dict['Geobacter'])
# yield_expr_rhod = symbol_sum(sum_dict['Rhodoferax'])

# alternation = 'reaction'
# alternation = 'interaction'

### do the simulation based on the objective that you choose

# # First, try to find the minimal number of active uptakes for the community (no interactiion)


if obj_num == 1:
    model.objective = competitions
    model.objective_direction = 'min'
    # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)

if obj_num == 2:
    model.objective = cooperations
    model.objective_direction = 'max'
    # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)

if obj_num==6: #('minimize abiotic nutrients for two species')
    # model.objective = num_upt - competitions - cooperations
    # model.objective_direction = 'min'
    constrain_abiotics(model)


sol = model.optimize()
solution = sol.objective_value
    # model.objective_direction = 'max'
    # solution_max=sol.objective_value
stop_crit = int(solution - 1)
print("MAX COOP NUMBER IS {}".format(solution))
# if obj_num == 3:
#     model.objective = competitions
#     model.objective_direction = 'max'
#     # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)
#
# if obj_num == 4:
#     model.objective = cooperations
#     model.objective_direction = 'min'
#     # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)
#

# if obj_num==1:
#     model.objective = num_upt
#     model.objective_direction = 'min'
#     # data_path = 'ilp_solutions/{}_min_tot_upt_{}.csv'.format(comm_name,alternation)
#
#
# # # Second, try to minimize the total abiotic uptake by the community (uptakes - interactions)
# if obj_num==2:
#     model.objective = num_upt - competitions - cooperations
#     model.objective_direction = 'min'
#     # data_path = 'ilp_solutions/{}_min_abiotic_upt_{}.csv'.format(comm_name,alternation)
#
#
# # # Third, try to minimize the total abiotic uptake by the community (uptakes - cooperations)
# if obj_num==3:
#     model.objective = num_upt - cooperations
#     model.objective_direction = 'min'
#     # data_path = 'ilp_solutions/{}_min_upt_compt_{}.csv'.format(comm_name,alternation)

# # # Forth, try to maximize the cooperations and yield of both organisms
# if obj_num==4:
#     model.objective = yield_expr_geo + yield_expr_rhod + cooperations
#     model.objective_direction = 'max'
#     # data_path = 'ilp_solutions/{}_max_yield_coop_{}.csv'.format(comm_name,alternation)
#
# # # Fifth, try to minimize the competitions and maximize yield of both organisms
# # model.objective = yield_expr_geo + yield_expr_rhod - competitions
# # model.objective_direction = 'max'
# # data_path = 'ilp_solutions/{}_max_yield_min_compt_{}.csv'.format(comm_name,alternation)
#
# # Sixth, try to minimize the uptakes and maximize yield of both organisms
# model.objective = 10 * yield_expr_geo + 10 * yield_expr_rhod - num_upt
# model.objective_direction = 'max'
# data_path = 'ilp_solutions/{}_max_yield_min_upt_{}.csv'.format(comm_name, alternation)


# output file to save the results
data_path = '/remind/projects/bee_project/111023_2member_ilp_solutions_031022_{}/obj_num_{}_alternatives_{}'.format(comm_name, obj_num,
                                                                                                 alternation)
if not os.path.exists(data_path):
    os.makedirs(data_path)

#
# data_path = '/remind/projects/bee_project/ilp_solutions_060822_{}'.format(comm_name)
# if not os.path.exists(data_path):
#     os.makedirs(data_path)


### solving and finding alternatives
# checking for feasibility and generating alternative solutions
#shouldnt be necessary
growth_ids = {
    'Gapicola': 'Growth',
    'Salvi': 'Growth',

}


growth_limits = {
    'Gapicola': 0.2,
    'Salvi': 0.2,

}

#
# sols = find_alternative_solutions(model, growth_ids, growth_limits,
#                                   alternation=alternation, max_alts=100,check_feasibility=False)
#
#
# sols.to_hdf(data_path+'/alternative_for_diff_obj_{}_alternation_{}_.h5'.format(obj_num,alternation),key='s')
#



sols_all,sols_active = find_alternative_solutions(model, growth_ids, growth_limits,
                                  alternation=alternation, max_alts=max_alternative,check_feasibility = False,stop_cond=True,criteria=stop_crit)


altsno=sols_active.alternative.nunique()

sols_all.to_hdf(data_path+'/sols_all_alternative_for_diff_obj_{}_alternation_{}_altno_{}_.h5'.format(obj_num,alternation,altsno),key='s')


sols_active.to_hdf(data_path+'/sols_active_alternative_for_diff_obj_{}_alternation_{}_altno_{}_.h5'.format(obj_num,alternation,altsno),key='s')

#

