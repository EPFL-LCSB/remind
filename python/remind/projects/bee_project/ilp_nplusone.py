#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: ReMIND Team
"""
#build models
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
from sys import argv

from pytfa.io.json import load_json_model as tfa_load_json_model
from pytfa.io.json import save_json_model as tfa_savejson_model
"Should be generalized for any size"
import pandas as pd
from remind.io.json import load_json_model
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse
from remind.optim.alternatives import find_alternative_solutions

from pytfa.optim.utils import symbol_sum
#
import os
from sys import argv
from remind.core.medium import constrain_abiotics


# _, number=argv

number=0
def build_ilp_model_pair(models_of_int):
    no_combin=len(models_of_int)
    #TODO YOU SHOULD NOT PUT SPECIES ID WITH UNDERSCORE
    #TODO IN THE SPECIES DATA IF YES IT WILL GIVE ERROR!!!!!
    models_dict={'Bifidobacterium_asteroides_PRL2011GF':"Bifido",
                 'Gilliamella_apicola_wkB1GF':"Gapicola",
                 'Lactobacillus_apis_Hma11GF':"Lapis",
                 'Lactobacillus_kullabergensis_Biut2GF':"Lkulla",
                 'Lactobacillus_mellifer_Bin4':"Lmellifer",
                 'Lactobacillus_mellis_Hon2':"Lmellis",
                 'Snodgrassella_alvi_wkB2GF':"Salvi"}


    " models of int as list"
    # data_path='/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'
    data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_141122_bee.h5'
    data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_111023_bee.h5'

    frame_combined_all=pd.read_hdf(data_path)
    "this will be the input to the function"
    frame_ilp=frame_combined_all[frame_combined_all.model.isin(models_of_int)]
    frame_ilp['yield'] = frame_ilp['yield_perc']
    frame_ilp['element'] = ['ONE'] * frame_ilp.shape[0]
    model_path='/remind/models/bee_models/core_members/tfa_real_010922/{}.json'
    model_path = '/remind/models/bee_models/core_members/tfa_latest_111023/{}.json'

    # defining the setup
    setup = {
        'sense': 'any',  # the sense of the objective for the interaction
        'diff_interaction': True,  # to differentiate between cooperation and competition
    }

    species_data={}
    for k in range(no_combin):
        model_id=models_dict[models_of_int[k]]
        species_data[model_id]={'model_name_json': 'tmodel_{}'.format(models_of_int[k]),
             'model_name':'{}'.format(models_of_int[k]),
             'growth': 0.2}

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
    return model



comm_name = '2_member_bee'
# model = load_json_model('../projects/community_models/{}_community.json'.format(comm_name))

# alternation_opts=['reaction','interaction']
alternation_opts=['reaction','interaction','interaction_positive','interaction_negative','abiotic_sources',\
                  'reaction_yield','interaction_yield','interaction_positive_yield','interaction_negative_yield','abiotic_sources_yield','interaction_positive_directional_general', 'interaction_positive_directional_general_yield']

model_list=[
            ['Gilliamella_apicola_wkB1GF', 'Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'],
            # ['Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'],
            # ['Gilliamella_apicola_wkB1GF', 'Bifidobacterium_asteroides_PRL2011GF'],
            # ['Gilliamella_apicola_wkB1GF', 'Snodgrassella_alvi_wkB2GF'],
            ]

max_alternative=5000
num=int(number)
models_of_int=model_list[num]


#for min abiotic sources and witht he corresponding integer cut
obj_num=2
alternation=alternation_opts[10]
# obj_num=6
# # alternation=alternation_opts[9]
# alternation=alternation_opts[4]

model = build_ilp_model_pair(models_of_int)
print("Model built for pair {}".format(models_of_int))

### setting different objective functions for the community
# the number of competitions and cooperations in the community
pos_int = model.get_variables_of_type(PositiveInteraction)
neg_int = model.get_variables_of_type(NegativeInteraction)
competitions = symbol_sum([v for v in neg_int])
cooperations = symbol_sum([v for v in pos_int])

# total number of reactions
num_upt = symbol_sum(model.get_variables_of_type(UptakeActivation))
num_sec = symbol_sum(model.get_variables_of_type(SecretionActivation))


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


#solve to store
sol=model.optimize()
solution=sol.objective_value
stop_cond=int(7)
#
# data_path = '/remind/projects/bee_project/ilp_solutions_111022_nplusone_{}/obj_num_{}_modelofint_{}_alternatives_{}'.format(comm_name, obj_num,
#                                                                                                  num,alternation)
#


data_path = '/remind/projects/bee_project/131023_ilp_solutions_280823_nplusone_{}/obj_num_{}_modelofint_{}_alternatives_{}'.format(comm_name, obj_num,
                                                                                                 num,alternation)


if not os.path.exists(data_path):
    os.makedirs(data_path)

growth_ids = {
    'Gapicola': 'Growth',
    'Salvi': 'Growth',
    'Bifido': 'Growth',

}


growth_limits = {
    'Gapicola': 0.2,
    'Salvi': 0.2,
    'Bifido': 0.2,

}



sols_all,sols_active = find_alternative_solutions(model, growth_ids, growth_limits,
                                  alternation=alternation, max_alts=max_alternative,check_feasibility = False,stop_cond=True,criteria=stop_cond)


altsno=sols_active.alternative.nunique()

sols_all.to_hdf(data_path+'/sols_all_alternative_for_diff_obj_{}_alternation_{}_altno_{}_.h5'.format(obj_num,alternation,altsno),key='s')


sols_active.to_hdf(data_path+'/sols_active_alternative_for_diff_obj_{}_alternation_{}_altno_{}_.h5'.format(obj_num,alternation,altsno),key='s')

#

