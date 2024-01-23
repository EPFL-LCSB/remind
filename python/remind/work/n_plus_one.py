#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from remind.io.json import load_json_model
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse
from remind.optim.alternatives import find_alternative_solutions

from remind.utils.utils import find_environment
from remind.analysis.perturbation import add_rem_member
    
from pytfa.optim.utils import symbol_sum


'''
This script performs a specific analysis where a new member is added to the 
community and its impact on the interactions is investigated. The main assumption
is that the medium provided to the n community is minimal in the sense that all
metabolites available in the environment is being used. Thus, when we move to
the n+1 community, the environment must remain the same.
However, we do not apply this assumption as a hard constraint and instead we formulate 
it as part of the objective function.

The script has three parts:
    1) find the solution for the n community (ref. solution)
    2) find the abiotic metabolites (specify the environment)
    3) add new varibales and constraints and formulate the objective for the n+1 community
'''


### finding the ref. solution
n_comm_name = 'Geobacter_Rhodoferax'
n_model = load_json_model('../projects/community_models/{}_community.json'.format(n_comm_name))

pos_int = n_model.get_variables_of_type(PositiveInteraction)
neg_int = n_model.get_variables_of_type(NegativeInteraction)
competitions = symbol_sum([v for v in neg_int]) 
cooperations = symbol_sum([v for v in pos_int])
# total number of reactions
num_upt = symbol_sum(n_model.get_variables_of_type(UptakeActivation))
num_sec = symbol_sum(n_model.get_variables_of_type(SecretionActivation))     
# take out solutions that support higher or lower yields
yield_vars = n_model.get_variables_of_type(YieldUse)
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
# Assumed objective function: try to minimize the total abiotic uptake by the community (uptakes - interactions)
n_model.objective = num_upt - competitions - cooperations
n_model.objective_direction = 'min'
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
sols = find_alternative_solutions(n_model, growth_ids, growth_limits, 
                                 max_alts=1)

### find the boitic resources
environment_ids = find_environment(n_model, sols)

###
n_plus_1_comm_name = 'Geobacter_Rhodoferax_Shewanella'
n_plus_1_model = load_json_model('../projects/community_models/{}_community.json'.format(n_plus_1_comm_name))
add_rem_member(n_plus_1_model, environment_ids)