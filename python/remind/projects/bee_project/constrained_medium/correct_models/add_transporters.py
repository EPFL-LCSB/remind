'''

author:DiME TEAM
'''


from pytfa.thermo.tmodel_struct import ThermoModelStructure
import os
from os.path import join
import cobra.test
from cobra.io import read_sbml_model
import time
from pytfa.optim.utils import symbol_sum
from pytfa.optim.constraints import ModelConstraint
import pandas as pd
import numpy as np
from sys import argv
from pytfa.thermo.utils import check_transport_reaction
from cobra import Reaction
from pytfa.io.json import save_json_model, load_json_model,load_json_model_tmodel_struct,save_json_model_tmodel_struct
from pytfa.redgem.debugging import check_BBB_production
import glob
# todo if not carbon source and if have the same formula dont add them binary variables

_, model_num=argv
model_no=int(model_num)
# #

constrain_fluxes=True
# model_name="Bifidobacterium_asteroides_PRL2011GF"
# model_name="Gilliamella_apicola_wkB1GF"
# model_name="Snodgrassella_alvi_wkB2GF"

# model_name='Lactobacillus_kullabergensis_Biut2GF' #3
# model_name='Lactobacillus_mellifer_Bin4' #4
# model_name='Snodgrassella_alvi_wkB2GF' #6


# model_no=0

biomass_rxn_id='Growth'
growth_limit=0.2
#todo put imm and iee scripts in other folders

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = CPLEX

aerobic=True
#output_file='./imm_alternatives/parallel_salvi_all_030921_manual_curation'

output_file_mat='/remind/models/bee_models/core_members/tfa_structures_corrected_mat/'
if not os.path.exists(output_file_mat):
    os.makedirs(output_file_mat)

output_file_json='/remind/models/bee_models/core_members/tfa_structures_corrected_json/'
if not os.path.exists(output_file_json):
    os.makedirs(output_file_json)

path='/remind/models/bee_models/core_members/tfa_structures'
allFiles = glob.glob(path + "/*"+"json")
biomass_rxn='Growth'

model = load_json_model_tmodel_struct(allFiles[model_no])
allFolders = os.listdir(path)
model_name=model.name.split('tfa_')[-1]
# native_exchanges=len(model.exchanges)

import re
def listcarbonsources_list(model,list_mets) :
    listPotCarbonSources=[]
    # you dont want Cl or cOBALT but c followed by a capital letter or number
    pattern = r"[C]\d|[C][A-Z]"
    for met in list_mets:
        m=model.metabolites.get_by_id(met)
        result=re.findall(pattern,m.formula)
        if any(result):
            listPotCarbonSources.append(m.id)
                   # model.reactions.get_by_id(r.id).lower_bound=-20

    return listPotCarbonSources

bbb_list=[k.id for k in model.reactions.Growth.reactants]
c_bbb=listcarbonsources_list(model,bbb_list)
inorg_bbb=[k for k in bbb_list if k not in c_bbb]
inorg=[k.replace('_c','_e') for k in inorg_bbb]





#for mn+2
#{'MNt2': 'h_e + mn2_e --> h_c + mn2_c'}
if model_name=='Lactobacillus_kullabergensis_Biut2GF':
    new_mn2=Reaction(id='MNt2pp')
    new_mn2.add_metabolites({model.metabolites.h_e:-1,model.metabolites.mn2_p:-1,\
                         model.metabolites.h_c:1,model.metabolites.mn2_c:1})
    model.add_reaction(new_mn2)

if model_name=='Lactobacillus_mellifer_Bin4':
    #todo
    new_zn2=Reaction(id='ZNabcpp')
    new_zn2.add_metabolites({model.metabolites.atp_c:-1,model.metabolites.zn2_p:-1,\
                         model.metabolites.h2o_c:-1,model.metabolites.adp_c:1,model.metabolites.zn2_c:1,\
                         model.metabolites.h_c:1,model.metabolites.pi_c:1})

    # new_zn2.add_metabolites({model.metabolites.zn2_p: -1, model.metabolites.zn2_c: 1})
    model.add_reaction(new_zn2)

if model_name == 'Snodgrassella_alvi_wkB2GF':
    # todo
    new_ca2=Reaction(id='CA2abcpp')
    new_ca2.add_metabolites({model.metabolites.atp_c:-1,model.metabolites.ca2_p:-1,\
                         model.metabolites.h2o_c:-1,model.metabolites.adp_c:1,model.metabolites.ca2_c:1,\
                         model.metabolites.h_c:1,model.metabolites.pi_c:1})

    model.add_reaction(new_ca2)

# "script to add reactions"
# # # 'add ca2 transporter by abc system potassium can also be an option'
# # new=Reaction(id='CA2abc')
# # new.add_metabolites({mytfa.metabolites.atp_c:-1,mytfa.metabolites.ca2_e:-1,\
# #                      mytfa.metabolites.h2o_c:-1,mytfa.metabolites.adp_c:1,mytfa.metabolites.ca2_c:1,\
# #                      mytfa.metabolites.h_c:1,mytfa.metabolites.pi_c:1})
# # model.add_reaction(new)

if constrain_fluxes:
    for k in model.exchanges:
        if k.id!='EX_o2_e':
            k.bounds=(-25,25)



if model_name=="Snodgrassella_alvi_wkB2GF":
    model.reactions.EX_o2_e.bounds=(-25,-1)

if 'EX_glc__D_e' in model.reactions:
    model.reactions.get_by_id('EX_glc__D_e').upper_bound=0

if 'EX_fru_e' in model.reactions:
    model.reactions.get_by_id('EX_fru_e').upper_bound=0
#
#
filepath_mat=output_file_mat +'/tmodel_{}.mat'.format(model_name)
filepath_json=output_file_json +'/tmodel_{}.json'.format(model_name)

save_json_model_tmodel_struct(model,filepath_json)

from pytfa.io.base import write_matlab_model_struct
write_matlab_model_struct(model, filepath_mat, varname='tmodel')
