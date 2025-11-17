
from pytfa.thermo.tmodel_struct import ThermoModelStructure
import os
from os.path import join
from cobra.io import read_sbml_model
import time
from remind.core.exchanges import imeanalysis_all_size,imeanalysis_once_size
from pytfa.optim.utils import symbol_sum
from pytfa.optim.constraints import ModelConstraint
from remind.core.parsimonious import *
import pandas as pd
import numpy as np
from sys import argv
from pytfa.thermo.utils import check_transport_reaction
from cobra import Reaction
from pytfa.io.json import save_json_model, load_json_model,load_json_model_tmodel_struct,save_json_model_tmodel_struct
import glob
from remind.core.medium import constrain_sources,constrain_uptake
from optlang.exceptions import SolverError
import cobra
# #



# _, model_num,max_flux,=argv
# model_no=int(model_num)
# max_flux_lim=float(max_flux)
# # tolerance_lb=float(tol)
# tolerance=tolerance_lb+1e-8

model_no= 6
max_flux_lim=1.0
tolerance_lb=1.11
# tolerance=tolerance_lb+1e-8


max_alternative=5
growth_limit=0.2
biomass='Growth'


tolerance_lb=1.0
tolerance=tolerance_lb+1e-8


max_secretion=25
max_it=5
max_consump=2
# model_no=1

constrain_fluxes=True
constrain_not_c_source_fluxes=True
#todo put imm and iee scripts in other folders

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = CPLEX


t=time.time()


#if started from an empty list possible
#TODO check_here
parsimonious='all_uptakes_for_list'
#parsimonious='c_uptakes'


# path='/remind/models/bee_models/core_members/tfa_real'
# path='/remind/models/bee_models/core_members/tfa_real_010922'
path='/remind/models/bee_models/core_members/tfa_real_101023'
path='/remind/models/bee_models/core_members/tfa_latest_111023'

# path='/remind/models/bee_models/core_members/tfa_structures_corrected'

allFiles = glob.glob(path + "/*"+"json")
biomass_rxn='Growth'
model = load_json_model(allFiles[model_no])

# model = load_json_model_tmodel_struct(allFiles[model_no])
model_name=model.name.split('tfa_')[-1]
native_exchanges=len(model.exchanges)
data_file='/remind/projects/bee_project/stats/'
# frame_screened=pd.read_csv(data_file+'/frame_screened_w_exp_data_stats_feasibility_checked_removed_duplicates_2.csv',index_col=0)
# frame_possible=pd.read_csv(data_file+'/possible_carbon_sources_bee_empty_emy.csv',index_col=None)
# frame_possible=pd.read_csv(data_file+'/possible_carbon_sources_initial.csv',index_col=None)
# frame_possible=pd.read_csv(data_file+'/possible_carbon_sources_bee.csv',index_col=None)
frame_possible=pd.read_csv(data_file+"possible_c_sources_list_300822_initial.csv",index_col=None)

# frame_possible.to_csv(data_file+'/possible_carbon_sources_bee_empty_{}.csv'.format(model_name))
# possible_c_sources=pd.read_csv(data_file+'/possible_carbon_sources_bee.csv',index_col=None)
# possible_c_sources=pd.read_csv(data_file+'/possible_carbon_sources_bee_3.csv',index_col=None)
possible_c_sources=pd.read_csv(data_file+"possible_c_sources_list_300822_initial.csv",index_col=None)

# possible_c_sources=pd.read_csv(data_file+'/possible_carbon_sources_initial.csv',index_col=None)

# list_all=[k for k in frame_screened.mets.unique()]
list_possible=[k for k in frame_possible.annotation_bigg.unique()]

# essential_subs=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_subs_090822.h5')

# essential_subs=pd.read_csv(data_file+'/not_main_c_source.csv')
# list_es=[k for k in essential_subs.annotation_bigg if k!='no']

essential_subs=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_subs_170822.h5')

list_es=[]
for k in essential_subs.essential_subs:
    list_es+=[i for i in k]

list_essential=list(set(list_es))
not_main_c_source=[k for k in list_essential if k not in possible_c_sources.annotation_bigg.unique()]


# not_main_c_source+=['EX_ser__D_e','EX_ala__D_e', 'EX_glu__L_e',"EX_arg__L_e","EX_ala__L_e"]
# not_main_c_source+=['EX_mbdg_e', 'EX_madg_e', 'EX_gam_e']
# not_main_c_source+=['EX_mbdg_e', 'EX_madg_e',"EX_23cgmp_e"]
# not_main_c_source+=['EX_mbdg_e', 'EX_madg_e', 'EX_gam_e','EX_23cgmp_e']

if model_name=="Snodgrassella_alvi_wkB2GF":
    model.reactions.EX_o2_e.bounds=(-25,-1)
    # max_flux_lim=0.01


# if model_name=='Lactobacillus_mellifer_Bin4':
#     # model.reactions.EX_o2_e.bounds=(-25,-1)
#     max_flux_lim=1.1

#here put elliott carbon sources


it=1
'''here put list remove'''
while it<max_it:
    # model = load_json_model_tmodel_struct(allFiles[model_no])
    model = load_json_model(allFiles[model_no])
    allFolders = os.listdir(path)
    model_name=model.name.split('tfa_')[-1]
    # frame_possible=pd.read_csv('/remind/projects/bee_project/stats/possible_carbon_sources_bee_empty.csv',index_col=0)
   #to do for each model separately
    # frame_possible = pd.read_csv(data_file + '/possible_carbon_sources_initial.csv', index_col=None)

    # frame_possible=pd.read_csv('/remind/projects/bee_project/stats/possible_carbon_sources_bee_empty_3.csv')
    frame_possible=pd.read_csv(data_file+"possible_c_sources_list_300822_initial_initial.csv",index_col=None)

    list_possible = [k for k in frame_possible.annotation_bigg.unique()]



    if constrain_fluxes:
        for k in model.exchanges:
            if k.id!='EX_o2_e':
                k.bounds=(-25,max_secretion)


    # if 'EX_glc__D_e' in model.reactions:
    #     model.reactions.get_by_id('EX_glc__D_e').upper_bound=0
    #
    # if 'EX_fru_e' in model.reactions:
    #     model.reactions.get_by_id('EX_fru_e').upper_bound=0

    # #context dependent
    # if 'EX_cit_e' in model.reactions:
    #     model.reactions.get_by_id('EX_cit_e').upper_bound=0
    #
    # if 'EX_icit_e' in model.reactions:
    #     model.reactions.get_by_id('EX_icit_e').upper_bound=0
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

    apply_solver_settings(model)



    limit=growth_limit
    model.reactions.get_by_id(biomass).upper_bound=limit
    model.reactions.get_by_id(biomass).lower_bound=limit*0.99


    list_carbon_sources = listcarbonsources(model)
    list_catabolite=[k for k in list_possible if k in list_carbon_sources ]
    if list_catabolite:
        constrain_sources(model,list_catabolite, max_cnsmp=max_consump, min_cnsmp=0)

    # 'just for information'
    list_carbon_sources = listcarbonsources(model)
    rxns_to_consider=[model.reactions.get_by_id(k) for k in list_carbon_sources if k not in list_possible]
    # rxns_to_consider=[model.reactions.get_by_id(k) for k in list_carbon_sources]


    if constrain_not_c_source_fluxes:
        for k in rxns_to_consider:
            if k.id!='EX_o2_e':
                k.bounds=(-25,max_secretion)


    if model_name == "Snodgrassella_alvi_wkB2GF":
        model.reactions.get_by_id('EX_fru_e').lower_bound = 0
        model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
        model.reactions.get_by_id('EX_sbt__D_e').lower_bound = 0
        # model.reactions.EX_cit_e.upper_bound = -1.0

    # list_exchanges=[k.id for k in mytfa2.exchanges]
    # list_c_reverse= [mytfa2.reactions.get_by_id(k).reverse_variable.name for k in list_carbon_sources]

    'generate untill max alternative'
    #mytfa2.exchanges

    # my_list,df_all=imeanalysis_all_size(model,model.exchanges,getAlternatives=True,max_alternative=max_alternative,biomass=biomass,stop_cond=True,criteria=2)
    # my_list,df_all=imeanalysis_all_size(model,rxns_to_consider,getAlternatives=True,max_alternative=max_alternative,biomass=biomass,stop_cond=True,criteria=0)
    # obj_min,obj=imeanalysis_once_size(model,rxns_to_consider)
    #
    #
    # #expr
    # from pytfa.optim.variables import ForwardBackwardUseVariable, ModelVariable, ReactionVariable
    #
    # BFUSE_VARS = model._var_kinds[ForwardBackwardUseVariable.__name__]
    # BFUSE_var_names = [x.name for x in BFUSE_VARS]
    #
    # expr_min_essential = sum(BFUSE_VARS)



    # # add model constraint
    # cons0 = model.add_constraint(ModelConstraint,
    #                               model,
    #                               expr_min_essential,
    #                               id_='min_essential_constraint',
    #                                 #fix uptake
    #                               lb=obj*0.99,
    #                               ub=obj)
    # model.repair()

    from remind.core.medium import constrain_uptake
    rxns_c=[model.reactions.get_by_id(k) for k in list_carbon_sources]

    try:
        block_essential_list=[k for k in not_main_c_source if k in list_carbon_sources]
        #convert  block essentials to True
        #here consider all
        # my_list_c_source,df_all_c_source=constrain_uptake(model,rxns_c,getAlternatives=True,max_alternative=50,stop_cond=True,block_essentials=False,essential_list=block_essential_list,biomass='Growth',max_flux_lim=max_flux_lim)
        my_list_c_source,df_all_c_source=constrain_uptake(model,rxns_to_consider,getAlternatives=True,max_alternative=50,stop_cond=True,block_essentials=True,essential_list=block_essential_list,biomass='Growth',max_flux_lim=max_flux_lim)
        # my_list_c_source,df_all_c_source=imeanalysis_all_size(model,rxns_to_consider,getAlternatives=True,max_alternative=max_alternative,biomass=biomass,stop_cond=True,criteria=0)
        frame_c_source = pd.concat(my_list_c_source, ignore_index=True)
        # frame['uptake_min']=np.ones(frame.shape[0])*sol_pars.objective_value
        # frame['yield_perc']=np.ones(frame.shape[0])*1.0/tolerance_lb

        alternative = frame_c_source.alternative.max()
        elapsed = time.time() - t
        print('time for ime analysis', elapsed)

        from remind.utils.postprocessing import *

        add_size(frame_c_source, ['alternative'], 'alt_size')
        add_list=list(frame_c_source[frame_c_source.MFAS == 1].metabolites.unique())
        print("For model {} essential is {}".format(model_name,add_list))
        list_possible += add_list
        it+=1

        if len(add_list)==0:
            it=max_it

        list_pos=pd.DataFrame(list_possible,columns=['annotation_bigg'])
        # list_pos.to_csv('/remind/projects/bee_project/stats/possible_carbon_sources_initial.csv')
        # list_pos.to_csv('/remind/projects/bee_project/stats/possible_c_sources_list_300822_initial.csv')

        # list_pos.to_csv('/remind/projects/bee_project/stats/possible_carbon_sources_bee_empty_3.csv')




    except SolverError:

        max_consump+=1
        print('Max consump constraint for the model {} is increased now {}'.format(model.name,max_consump))



#


print("Max consumption considering CATABOLITE REPRESSION is {}".format(max_consump))
#
# """start"""

### this is to generate essential substrates#
max_uptake=25
max_secretion=25
model = load_json_model(allFiles[model_no])
# add this reaction for snod
#FROM HERE ON DO THE DIFFERENCE
model_name=model.name.split('tfa_')[-1]
data_dir_cobra= '/remind/models/bee_models/core_members/carveme_models'

model_cobra=cobra.io.read_sbml_model(join(data_dir_cobra,"{}.xml".format(model_name)))
diff_rxns=[k for k in model.reactions if k not in model_cobra.reactions]
sink_diff=[k for k in diff_rxns if "sink" in k.id]
model.remove_reactions(sink_diff)
list_carbon_sources = listcarbonsources(model)
inorganics=[k for k in model.exchanges if k.id not in list_carbon_sources]
inorganics=[k for k in inorganics if k.id!="EX_o2_e"]

# sink_carbon=[k.replace('EX_','sink_') for k in listcarbonsources(model)]
# sink_carbon=[k.replace('_e','_c') for k in sink_carbon]
# c_source_all_rxns=listcarbonsources_all_rxns(model)
# sink_close=[k for k in sink_diff if k.id in c_source_all_rxns]
# other_sink=[k for k in sink_diff if k not in sink_close]

apply_solver_settings(model)

if model_name=="Snodgrassella_alvi_wkB2GF":
    model.reactions.EX_o2_e.bounds=(-max_uptake,-1)
if constrain_fluxes:
    for k in model.exchanges:
        if k.id != 'EX_o2_e':
            if k.id in list_carbon_sources:
                k.bounds = (-1*max_flux_lim, max_secretion)
        if k.id in list_possible:
            k.bounds = (-max_uptake, max_secretion)


# model.reactions.EX_thm_e.bounds=(-max_uptake,max_secretion)
for rxn_inor in inorganics:
    rxn_inor.bounds= (-max_uptake,50)

if 'EX_glc__D_e' in model.reactions:
    model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0

if 'EX_fru_e' in model.reactions:
    model.reactions.get_by_id('EX_fru_e').upper_bound = 0


if 'EX_cit_e' in model.reactions:
    model.reactions.get_by_id('EX_cit_e').upper_bound = 0

if 'EX_icit_e' in model.reactions:
    model.reactions.get_by_id('EX_icit_e').upper_bound = 0


if model_name=="Snodgrassella_alvi_wkB2GF":
    model.reactions.get_by_id('EX_fru_e').lower_bound = 0
    model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
    model.reactions.get_by_id('EX_sbt__D_e').lower_bound = 0


# if model_name=="Snodgrassella_alvi_wkB2GF":
#     model.reactions.EX_o2_e.bounds=(-25,-1)
#     max_flux_lim=0.01

limit=growth_limit
model.reactions.get_by_id(biomass).upper_bound=limit
model.reactions.get_by_id(biomass).lower_bound=limit*0.99



list_carbon_sources = listcarbonsources(model)
list_catabolite=[k for k in list_possible if k in list_carbon_sources ]
rxns_to_consider = [model.reactions.get_by_id(k) for k in list_carbon_sources if k not in list_possible]
# max_consump=4
constrain_sources(model,list_catabolite, max_cnsmp=max_consump, min_cnsmp=0)

# model.reacti[ons.EX_meoh_e.bounds=(0,max_secretion)
# model.reactions.EX_mso3_e.bounds=(0,max_secretion)
model.repair()



# my_list,df_all=imeanalysis_all_size(model,model.exchanges,getAlternatives=True,max_alternative=max_alternative,biomass=biomass,stop_cond=True,criteria=2)
# my_list,df_all=imeanalysis_all_size(model,rxns_to_consider,getAlternatives=True,max_alternative=max_alternative,biomass=biomass,stop_cond=True,criteria=0)
obj_min,obj=imeanalysis_once_size(model,rxns_to_consider)


#expr
from pytfa.optim.variables import ForwardBackwardUseVariable, ModelVariable, ReactionVariable

BFUSE_VARS = model._var_kinds[ForwardBackwardUseVariable.__name__]
BFUSE_var_names = [x.name for x in BFUSE_VARS]

expr_min_essential = sum(BFUSE_VARS)

# add model constraint
cons0 = model.add_constraint(ModelConstraint,
                              model,
                              expr_min_essential,
                              id_='min_essential_constraint',
                                #fix uptake
                              lb=obj*0.999,
                              ub=obj)
model.repair()
#
# sink_carbon=[k.replace('EX_','sink_') for k in listcarbonsources(model)]
# sink_carbon=[k.replace('_e','_c') for k in sink_carbon]
# c_source_all_rxns=listcarbonsources_all_rxns(model)
# sink_close=[k for k in sink_diff if k.id in c_source_all_rxns]
#
# expr_sink=symbol_sum([k.flux_expression for k in sink_close])
# model.objective=expr_sink
# model.objective_direction="min"
# sol_sink=model.optimize()
# print("Min flux through sinks is {}".format(sol_sink.objective_value))


# cons1 = model.add_constraint(ModelConstraint,
#                               model,
#                               expr_sink,
#                               id_='pars_constraint',
#                                 #fix uptake
#
#                               ub=(sol_sink.objective_value)*(1.01))





#here maybe also to be tried if only c uptake is minimized is it the same
expr_pars, sol_pars = minimizealluptakes_by_weight(model, rxns_to_consider)
cons1 = model.add_constraint(ModelConstraint,
                              model,
                              expr_pars,
                              id_='pars_constraint',
                                #fix uptake
                              lb=(sol_pars.objective_value)*1.0,
                              ub=(sol_pars.objective_value)*(1.0+1e-4))


model.repair()

my_list,df_all=imeanalysis_all_size(model,rxns_to_consider,getAlternatives=True,max_alternative=50,biomass=biomass,stop_cond=True,criteria=0)
# #
frame = pd.concat(my_list, ignore_index=True)
alternative=frame.alternative.max()

frame_union=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr_updated_010922_TFA.h5')

prev=[k for k in frame_union[frame_union.model==model_name].metabolites.unique()]
print(len(prev)-frame[frame.alternative==1].shape[0])
print([k for k in prev if k not in frame.metabolites.unique() ])
print([k for k in frame.metabolites.unique() if k not in prev])



# """end"""
#
# # frame['uptake_min']=np.ones(frame.shape[0])*sol_pars.objective_value
# # frame['yield_perc']=np.ones(frame.shape[0])*1.0/tolerance_lb
#
# alternative=frame.alternative.max()
# elapsed = time.time() - t
# print('time for ime analysis', elapsed)
#
# from remind.utils.postprocessing import *
# add_size(frame,['alternative'],'alt_size')
# #
from remind.core.decomposition import find_essential_alternates_from_dime_list
a,b=find_essential_alternates_from_dime_list(frame)


# #
# # #
# # #
# # '''Saving data'''
# # #output file to save the results
# # output_file_all='/remind/projects/bee_project/2_Essential_needs_090822/{}/all_vars'.format(model.name.split('tfa_')[-1])
# # if not os.path.exists(output_file_all):
# #     os.makedirs(output_file_all)
# #



# output_file_alts='/remind/projects/bee_project/PRAY_061023_TFA_3_new_Essential_needs/{}/alternatives'.format(model.name.split('tfa_')[-1])
# if not os.path.exists(output_file_alts):
#     os.makedirs(output_file_alts)
#
# # # # #"store them in hdf5 files"
# frame.to_hdf(output_file_alts+'/alternative_mets_for_{}_growth_{}_alt_{}_pars_{}_tol_pars_{}_tol_lb_{}.h5'.format(model.id,limit,alternative,parsimonious,tolerance,tolerance_lb),key='s')
# # # df_all.to_hdf(output_file_all+'/alternatives_all_vars_for_{}_growth_{}_alt_{}_pars_{}_tol_pars_{}_tol_lb_{}.h5'.format(model.id,limit,alternative,parsimonious,tolerance,tolerance_lb),key='s')

#
#
#
#

#
#
# import cobra
# model = load_json_model(allFiles[model_no])
# apply_solver_settings(model)
# model_name=model.name.split('tfa_')[-1]
# data_dir_cobra= '/remind/models/bee_models/core_members/carveme_models'
# model_cobra=cobra.io.read_sbml_model(join(data_dir_cobra,"{}.xml".format(model_name)))
# diff_rxns=[k for k in model.reactions if k not in model_cobra.reactions]
# sink_diff=[k for k in diff_rxns if "sink" in k.id]
# ddf=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_80.0_pars_c_uptakes_only_C_moles_tol_pars_1.2222222222222223_tol_lb_1.1111111111111112.h5")
#
# ddf=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_129.0_pars_c_uptakes_only_C_moles_tol_pars_2.75_tol_lb_2.5.h5")
#
#
#
# ddf=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_533.0_pars_c_uptakes_only_C_moles_tol_pars_5.5_tol_lb_5.0.h5")
# # ddf=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/alternatives/alternative_mets_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_956.0_pars_c_uptakes_only_C_moles_tol_pars_11.0_tol_lb_10.0.h5")
# # ddf=pd.read_hdf("/remind/projects/bee_project/010922_correct_TFA_allowed_True_limited_carbon/Gilliamella_apicola_wkB1GF/all_vars/alternatives_all_vars_for_RelaxedModel_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_956.0_pars_c_uptakes_only_C_moles_tol_pars_11.0_tol_lb_10.0.h5")
# kk=ddf.loc[[k.id for k in sink_diff]]
# d=pd.DataFrame()
# d['min']=kk.min(axis=1)
# d['max']=kk.max(axis=1)
# d['mean']=kk.mean(axis=1)
# d['median']=kk.median(axis=1)
# d['sum']=kk.sum(axis=1)
# kk.max(axis=0)
#
# sink_carbon=[k.replace('EX_','sink_') for k in listcarbonsources(model)]
# sink_carbon=[k.replace('_e','_c') for k in sink_carbon]
# c_source_all_rxns=listcarbonsources_all_rxns(model)
# sink_close=[k for k in sink_diff if k.id in c_source_all_rxns]
# other_sink=[k for k in sink_diff if k not in sink_close]
# kk.sum()
#
# kk.loc[[k.id for k in sink_close]].sum(axis=0)