'''

author:DiME TEAM
'''


from pytfa.thermo.tmodel_struct import ThermoModelStructure
from pytfa.thermo.tmodel import ThermoModel
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
from pytfa.io import import_matlab_model, load_thermoDB, \
    read_lexicon, annotate_from_lexicon, \
    read_compartment_data, apply_compartment_data

import pytfa

# _, model_number,aerobic_=argv
# model_no=int(model_number)
# aerobic=bool(float(aerobic_))


model_no=0
aerobic=True

biomass_rxn_id='Growth'
growth_limit=0.2
constrain_fluxes=True
#todo put imm and iee scripts in other folders

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = CPLEX

biomass_rxn='Growth'
#output_file='./imm_alternatives/parallel_salvi_all_030921_manual_curation'

output_file='/remind/models/bee_models/core_members/tfa_real_010922/'
output_file='/remind/models/bee_models/core_members/tfa_latest_110123/'

if not os.path.exists(output_file):
    os.makedirs(output_file)

data_dir= '/remind/models/bee_models/core_members/carveme_models'
allFolders = os.listdir(data_dir)
model=cobra.io.read_sbml_model(join(data_dir,"{}".format(allFolders[model_no])))
model_name=allFolders[model_no].split(".xml")[0]
t=time.time()






"check atpm"
print("ATPM lower bound is {}".format(model.reactions.ATPM.lower_bound))

model.reactions.ATPM.lower_bound=3.0
# for k in model.exchanges:
#     k.bounds=(-1000,1000)

if constrain_fluxes:
    for k in model.exchanges:
        if k.id!='EX_o2_e':
            k.bounds=(-25,25)

else:
    for k in model.exchanges:
        k.bounds=(-1000,1000)

if 'EX_glc__D_e' in model.reactions:
    model.reactions.get_by_id('EX_glc__D_e').upper_bound=0

if 'EX_fru_e' in model.reactions:
    model.reactions.get_by_id('EX_fru_e').upper_bound=0


#new additional info
"it is known that it cannot catabolize carbohydrates"
if model_name=="Snodgrassella_alvi_wkB2GF":
    model.reactions.get_by_id('EX_fru_e').lower_bound = 0
    model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
    #to force o2 uptake
    model.reactions.EX_o2_e.bounds = (-25, -1)



"check length of exchanges and if there are sink or demands"
print('Number of exchange reactions is {}'.format(len(model.exchanges)))
rxn_exch=[k.id for k in model.reactions if len(k.metabolites)==1]
diff=[k for k in rxn_exch if k not in model.exchanges]
print('There are {} Demand or Sink Reactions which are {}'.format(len(diff),diff))

"aerobic or anaerobic"
if not aerobic:
    model.reactions.EX_o2_e.lower_bound=0.0

"add required transporters"
if model_name=='Lactobacillus_kullabergensis_Biut2GF':
    new_mn2=Reaction(id='MNt2pp')
    new_mn2.add_metabolites({model.metabolites.h_e:-1,model.metabolites.mn2_p:-1,\
                         model.metabolites.h_c:1,model.metabolites.mn2_c:1})
    model.add_reaction(new_mn2)


if model_name=='Lactobacillus_mellifer_Bin4':
    new_zn2=Reaction(id='ZNabcpp')
    new_zn2.add_metabolites({model.metabolites.atp_c:-1,model.metabolites.zn2_p:-1,\
                         model.metabolites.h2o_c:-1,model.metabolites.adp_c:1,model.metabolites.zn2_c:1,\
                         model.metabolites.h_c:1,model.metabolites.pi_c:1})

    # new_zn2.add_metabolites({model.metabolites.zn2_p: -1, model.metabolites.zn2_c: 1})
    model.add_reaction(new_zn2)



if model_name == 'Snodgrassella_alvi_wkB2GF':
    new_ca2=Reaction(id='CA2abcpp')
    new_ca2.add_metabolites({model.metabolites.atp_c:-1,model.metabolites.ca2_p:-1,\
                         model.metabolites.h2o_c:-1,model.metabolites.adp_c:1,model.metabolites.ca2_c:1,\
                         model.metabolites.h_c:1,model.metabolites.pi_c:1})

    model.add_reaction(new_ca2)



print("model is aerobic {} max uptake of o2 is set in the model to {} ".format(aerobic,model.reactions.EX_o2_e.lower_bound) )
'here thermomodel is created without thermo constraints but use constraints'

#force for Gilliamella that pyruvate secretion is on
#close cellobiose and try

"here build the real thermo model"

lexicon = pd.read_csv('/remind/models/thermo/salvi_gapicola_union.csv',encoding='latin-1',index_col=0)
# lexicon=lexicon.drop("atp")
for k in lexicon.index:
    lexicon.loc[k+'_c']=lexicon.loc[k]
    lexicon.loc[k+'_e']=lexicon.loc[k]
    lexicon.loc[k+'_p']=lexicon.loc[k]

#try to drop

compartment_data = read_compartment_data('/remind/models/thermo/compartment_data_salvi_gapicola.json')
print("Loading thermo data...")

thermo_data = load_thermoDB('/remind/models/thermo/salvi_gapicola.thermodb')
# Initialize the cobra_model
mytfa = pytfa.ThermoModel(thermo_data, model)



annotate_from_lexicon(mytfa, lexicon)
apply_compartment_data(mytfa, compartment_data)

# mytfa =ThermoModelStructure ([], model)
mytfa.name = 'tfa_{}'.format(model_name)
mytfa.solver = solver

# Solver settings
def apply_solver_settings(model, solver = solver):
    model.solver = solver
    # model.solver.configuration.verbosity = 1
    model.solver.configuration.tolerances.feasibility = 1e-9
    model.solver.configuration.tolerances.optimality = 1e-9
    model.solver.configuration.tolerances.integrality = 1e-9

    if solver == 'optlang_gurobi':
        model.solver.problem.Params.NumericFocus = 3
    model.solver.configuration.presolve = True

apply_solver_settings(mytfa)
mytfa.prepare()
mytfa.convert()
fba_solution = model.optimize()
print(fba_solution)
fba_value = fba_solution.objective_value

# Solver settings

## FBA

# fva = flux_variability_analysis(mytfa)

## TFA conversion

# med=['EX_h2o_e',
#  'EX_h_e',
#  'EX_cl_e',
#  'EX_pi_e',
#  'EX_nh4_e',
#  'EX_fe3_e',
#  'EX_k_e',
#  'EX_ca2_e',
#  'EX_cit_e',
#  'EX_mg2_e',
#  'EX_mn2_e',
#  'EX_cobalt2_e',
#  'EX_zn2_e',
#  'EX_cu2_e',
#  'EX_o2_e',
#  'EX_so4_e',
#  'EX_thm_e',
#  'EX_arg__L_e',
#  'EX_fol_e',
#  'EX_dha_e',
#      ]

# for rxn in mytfa.exchanges:
#     if rxn.id not in med:
#         rxn.lower_bound=0.0

# no_list=['EX_meoh_e','EX_mso3_e']#,'EX_ser__L_e']
#
# # no_list=[]
# for rxn in no_list:
#     mytfa.reactions.get_by_id(rxn).lower_bound=0.0

lower_bound_growth=0.4
# lower_bound_growth=fba_value*0.80
## Info on the cobra_model
mytfa.print_info()
from optlang.exceptions import SolverError

try:
    tfa_solution = mytfa.optimize()
    tfa_value = tfa_solution.objective_value
    if tfa_value<=0:
        print('Thermo model has zero growth')
        from pytfa.optim.relaxation import relax_dgo, relax_dgo_err
        #before 0.4
        mytfa.reactions.get_by_id(biomass_rxn).lower_bound = lower_bound_growth
        relaxed_model, slack_model, relax_table = relax_dgo(mytfa, in_place=False)

        original_model, mytfa = mytfa, relaxed_model

        print('Relaxation: ')
        print(relax_table)

        tfa_solution = mytfa.optimize()
        tfa_value = tfa_solution.objective_value
    else:
        print("Thermomodel is feasible")

except SolverError:
    print('Thermo model was not feasible')
    from pytfa.optim.relaxation import relax_dgo,relax_dgo_err

    mytfa.reactions.get_by_id(biomass_rxn).lower_bound = lower_bound_growth
    relaxed_model, slack_model, relax_table = relax_dgo(mytfa,in_place=False)

    original_model, mytfa = mytfa, relaxed_model

    print('Relaxation: ')
    print(relax_table)

    tfa_solution = mytfa.optimize()
    tfa_value = tfa_solution.objective_value


from pytfa.analysis.variability import _variability_analysis_element
mytfa.reactions.Growth.lower_bound=0.0
min_growth=_variability_analysis_element(mytfa,mytfa.reactions.Growth,"min")
max_growth=_variability_analysis_element(mytfa,mytfa.reactions.Growth,"max")
print("ATTENTION!!!!! Min growth is {} max growth is {}".format(min_growth,max_growth))

filepath = output_file + '/tmodel_{}.json'.format(model_name)
# save_json_model(mytfa, filepath)

"save the thermo model"

limit=growth_limit
mytfa.reactions.get_by_id(biomass_rxn_id).upper_bound=limit
mytfa.reactions.get_by_id(biomass_rxn_id).lower_bound=limit
check_BBB=check_BBB_production(mytfa, biomass_rxn_id, verbose = False)
if (check_BBB<=0).any():
    print("NOT ALL BBBs can be produced in this condition")
else:
    print("all BBBs can be produced")
    print(check_BBB)
#     filepath = output_file + '/tmodel_{}.json'.format(model_name)
#     save_json_model(mytfa, filepath)





# dd=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_956.0_pars_c_uptakes_only_C_moles_tol_pars_11.0_tol_lb_10.0.h5")
# dd=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_199.0_pars_c_uptakes_only_C_moles_tol_pars_1.1_tol_lb_1.0.h5")
#
# dd=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/alternatives/alternative_mets_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_199.0_pars_c_uptakes_only_C_moles_tol_pars_1.1_tol_lb_1.0.h5")
#
# dd=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_199.0_pars_c_uptakes_only_C_moles_tol_pars_1.1_tol_lb_1.0.h5")
# dd=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_80.0_pars_c_uptakes_only_C_moles_tol_pars_1.2222222222222223_tol_lb_1.1111111111111112.h5")
#
