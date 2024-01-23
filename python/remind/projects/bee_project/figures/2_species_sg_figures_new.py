import pandas as pd
import numpy as np
from figure_functions import *



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
import glob
from sys import argv
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

comm_name = '7_member_bee'
comm_name = '2_member_bee'
 #for 2speciespos int data #031022

model_of_int_list=[['Gilliamella_apicola_wkB1GF','Lactobacillus_mellifer_Bin4'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_mellifer_Bin4','Snodgrassella_alvi_wkB2GF'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_mellis_Hon2'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_mellis_Hon2','Snodgrassella_alvi_wkB2GF'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_apis_Hma11GF'],
            ['Gilliamella_apicola_wkB1GF','Lactobacillus_apis_Hma11GF','Snodgrassella_alvi_wkB2GF'],
            ['Gilliamella_apicola_wkB1GF', 'Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'],
            ['Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'],
            ['Gilliamella_apicola_wkB1GF', 'Bifidobacterium_asteroides_PRL2011GF'],
            ['Gilliamella_apicola_wkB1GF', 'Snodgrassella_alvi_wkB2GF'],
            ]

model_of_int_list=[
            ['Gilliamella_apicola_wkB1GF', 'Bifidobacterium_asteroides_PRL2011GF', 'Snodgrassella_alvi_wkB2GF'],


            ]
def add_label(df, group=[""],name='abiotic', overshoot=0):
    df['label_{}'.format(name)] = df.groupby(by=group).grouper.group_info[0] + overshoot

output_file_figures="/remind/projects/bee_project/constrained_medium/analysis/figures_updated"
if not os.path.exists(output_file_figures):
    os.makedirs(output_file_figures)


data_path = '/remind/projects/bee_project/ilp_solutions_031022_{}'.format(comm_name)



data_path = '/remind/projects/bee_project/241022_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/101122_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/101122_2_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_2_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_3_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_new_data_integ_ilp_solutions_2_member_bee'.format(comm_name)
data_path = '/remind/projects/bee_project/111122_newest_2_data_integ_ilp_solutions_2_member_bee'.format(comm_name)

data_path='/remind/projects/bee_project/060623_ilp_solutions_nplusone_2_member_bee/obj_num_2_modelofint_10_alternatives_interaction_positive_directional_general_yield/'

data_path="/remind/projects/bee_project/ilp_solutions_280823_nplusone_2_member_bee"



data_path="/remind/projects/bee_project/2member_081023_ilp_solutions_031022_2_member_bee"

data_path="/remind/projects/bee_project/abiotic_2member_081023_ilp_solutions_031022_2_member_bee"
allDirs = os.listdir(data_path)

obj_numbers = {1: 'min competition',
               2: 'max cooperation',
               3: 'max competition',
               4: 'min cooperation',
               5: "weighted",
               6: "min abiotic sources"}

obj_num_sep=True
list_ = []
for k in allDirs:
    # k.split('obj_num_')
    if obj_num_sep:
        obj_num = int(k.split('obj_num_')[1][0])
    if "alternatives" in k:
        alternation = k.split('alternatives_')[1]

    elif "alternation" in k:
        alternation = k.split('alternation_')[1]

    path_new = data_path + '/' + k
    allFiles = glob.glob(path_new + "/*" + "h5")
    for file in allFiles:
        if 'sols_active' in file:
            d = pd.read_hdf(file)
            if obj_num_sep:
                d['obj_num'] = [obj_num] * d.shape[0]
            d['alternation'] = [alternation] * d.shape[0]
            if "modelofint" in file:
                d['modelofint_no'] = [int(file.split('modelofint_')[1][0])] * d.shape[0]

                d['modelofint'] = tuple([model_of_int_list[int(file.split('modelofint_')[1][0])]]) * d.shape[0]
            list_.append(d)

# contains merged all groups
df = pd.concat(list_, ignore_index=True)

df['obj_num'] = [obj_numbers[k] for k in df.obj_num]
df.groupby(['obj_num', 'alternation']).objective.unique()
df.groupby(['obj_num', 'alternation']).alternative.nunique()
gb = df.groupby(['obj_num', 'alternation'])
keys = [k for k in gb.groups.keys()]
keys_int = [k for k in keys]
# keys_int=[k for k in keys if 'interaction' in k]

# newgrouping for n+1 studies
if "modelofint" in df.columns:
    gb = df.groupby(['modelofint'])
    keys = [k for k in gb.groups.keys()]
    keys_int = [k for k in keys]

    # groupofints = [keys_int[3], keys_int[4]]
    #
    # groupofints = [keys_int[5], keys_int[6]]
    # d = pd.concat([gb.get_group(gr) for gr in groupofints])
    grouping = ['alternation', 'modelofint','obj_num']
    gb = df.groupby(grouping)
    keys = [k for k in gb.groups.keys()]
    keys_int = [k for k in keys]

d_int=gb.get_group(keys_int[0])
d_int_direct=gb.get_group(keys_int[1])
d_int_yield=gb.get_group(keys_int[2])

models_dict={
             "model1":"Gapicola",
             "model2":"Salvi"}
frame_int=get_species_data(d_int,models_dict=models_dict,len_models=2)
frame_int=get_species_data(d_int_direct,models_dict=models_dict,len_models=2)
frame_int_yield=get_species_data(d_int_yield,models_dict=models_dict,len_models=2)