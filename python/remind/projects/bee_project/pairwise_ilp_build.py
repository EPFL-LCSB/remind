"put the function"
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
from sys import argv

from pytfa.io.json import load_json_model as tfa_load_json_model
from pytfa.io.json import save_json_model as tfa_savejson_model
"Should be generalized for any size"

_, combination=argv
combin=int(combination)



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
    # data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_141122_bee.h5'
    data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_081023_bee.h5'

    frame_combined_all=pd.read_hdf(data_path)
    "this will be the input to the function"
    frame_ilp=frame_combined_all[frame_combined_all.model.isin(models_of_int)]
    frame_ilp['yield'] = frame_ilp['yield_perc']
    frame_ilp['element'] = ['ONE'] * frame_ilp.shape[0]
    model_path='/remind/models/bee_models/core_members/tfa_real_010922/{}.json'

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

    # species_data = {
    #     models_dict[models_of_int[0]]:
    #         {'model_name_json': 'tmodel_{}'.format(models_of_int[0]),
    #          'model_name':'{}'.format(models_of_int[0]),
    #          'growth': 0.2},
    #     models_dict[models_of_int[1]]:
    #         {'model_name_json': 'tmodel_{}'.format(models_of_int[1]),
    #           'model_name': '{}'.format(models_of_int[1]),
    #          'growth': 0.2}
    # }
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
    return model

import pandas as pd

from pytfa.optim.utils import symbol_sum
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse
from remind.core.medium import *
"read the frame"
# data_path='/remind/projects/bee_project/stats/combined_dimes_unique_070922_bee.h5'
data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_141122_bee.h5'

# data_path = '/remind/projects/bee_project/stats/combined_dimes_unique_141122_bee.h5'

frame_combined_all = pd.read_hdf(data_path)
from itertools import combinations
# combin=2
all_pairs = list(combinations(frame_combined_all.model.unique(), combin))
all_pairs=[list(k) for k in all_pairs]

"read the model"


sol_dict = {}
ind = 0
for index_pair in range(len(all_pairs)):
    # with model as model_copy:
    models_of_int=list(all_pairs[index_pair])
    #perform with model.copy
    model_copy=build_ilp_model_pair(models_of_int)
    print("Copy model built for pair {}".format(models_of_int))
    pos_int = model_copy.get_variables_of_type(PositiveInteraction)
    neg_int = model_copy.get_variables_of_type(NegativeInteraction)
    competitions = symbol_sum([v for v in neg_int])
    cooperations = symbol_sum([v for v in pos_int])
    num_upt = symbol_sum(model_copy.get_variables_of_type(UptakeActivation))
    num_sec = symbol_sum(model_copy.get_variables_of_type(SecretionActivation))

    print("Currently performing for the pair {}".format(models_of_int))

    "1) minimize competition"
    model_copy.objective = competitions
    model_copy.objective_direction = 'min'
    sol_opt = model_copy.optimize()
    min_compet = sol_opt.objective_value
    print("Min competition for pair {} is {}".format(models_of_int,min_compet))



    "2) maximize cooperation"
    model_copy.objective = cooperations
    model_copy.objective_direction = 'max'
    sol_opt = model_copy.optimize()
    max_coop = sol_opt.objective_value
    print("Max cooperation for pair {} is {}".format(models_of_int, max_coop))


    # "3) minimize abiotic sources"
    # model_copy.objective = num_upt - competitions - cooperations
    # model_copy.objective_direction = 'min'
    # sol_opt = model_copy.optimize()
    # min_req = sol_opt.objective_value
    # print("Min abiotic sources for pair {} is {}".format(models_of_int, min_req))

    "3) minimize abiotic sources 2"
    # model_copy.objective = competitions
    # model_copy.objective_direction = 'min'
    constrain_abiotics(model_copy)
    sol_opt = model_copy.optimize()
    min_req = sol_opt.objective_value
    print("Min abiotic sources for pair {} is {}".format(models_of_int,min_req))

    "4) maximize competition"
    model_copy.objective = competitions
    model_copy.objective_direction = 'max'
    sol_opt = model_copy.optimize()
    max_compet = sol_opt.objective_value
    print("Max competition for pair {} is {}".format(models_of_int, max_compet))

    "5) Minimize cooperation"
    model_copy.objective = cooperations
    model_copy.objective_direction = 'min'
    sol_opt = model_copy.optimize()
    min_coop = sol_opt.objective_value
    print("Min cooperation for pair {} is {}".format(models_of_int, min_coop))

    "Save the data as a dictionary"
    sol_dict[ind] = {'pair': models_of_int,
                     'min_competition': min_compet,
                     'max_cooperation': max_coop,
                     'min_abiotic sources': min_req,
                     'max_competition': max_compet,
                     'min_cooperation': min_coop}
    ind += 1

# "convert the data into a dataframe"
sol_data = pd.DataFrame.from_dict(sol_dict, orient='index')
# sol_data.to_csv('/remind/projects/bee_project/constrained_medium/built_bee_pairwise_objectives_090922_combin{}.csv'.format(combin))

sol_data.to_csv('/remind/projects/bee_project/constrained_medium/built_bee_pairwise_objectives_091023_combin{}.csv'.format(combin))
