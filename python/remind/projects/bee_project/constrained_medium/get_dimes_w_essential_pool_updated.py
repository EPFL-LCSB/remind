import os
import glob
import os
import time
from sys import argv
from pytfa.optim.utils import symbol_sum

import numpy as np
import pandas as pd
from pytfa.io.json import load_json_model
from pytfa.optim.constraints import ModelConstraint

from remind.core.exchanges import imeanalysis_all_size
from remind.core.medium import constrain_sources
from remind.core.parsimonious import *


from os.path import join
import cobra

#these code is the one that were used for the latest analysis

_, model_num,tol=argv
model_no=int(model_num)

#
# model_no=0
# tol=1.0


tolerance_lb=1.0/float(tol)
tolerance=tolerance_lb*1.1 # a bit of gap given for the fixed yield_cut

#
# model_no=6
# tolerance_lb=1.0/0.80
# tolerance=tolerance_lb*1.1

max_alternative=5000
growth_limit=0.2
biomass='Growth'
allowed=True
max_consmp=3 #maxmum carbon sources to be used at the same time from the given list
max_flux_lim=1.0
constrain_fluxes=True
#todo put imm and iee scripts in other folders

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = CPLEX


t=time.time()

#TODO check_here
# parsimonious='all_uptakes_mw'
#to make the yield_constraint with respect to the carbon uptake
parsimonious='c_uptakes_only_C_moles' #yield is calc
#parsimonious='c_uptakes'


# path='/remind/models/bee_models/core_members/tfa_structures_corrected'
# path='/remind/models/bee_models/core_members/tfa_real'
path='/remind/models/bee_models/core_members/tfa_real_010922'
path='/remind/models/bee_models/core_members/tfa_real_111023'

allFiles = glob.glob(path + "/*"+"json")
biomass_rxn='Growth'

# model = load_json_model_tmodel_struct(allFiles[model_no])
model = load_json_model(allFiles[model_no])


print('Number of reactions {} and number of genes {} for model {}'.format(len(model.reactions),len(model.genes),model.id))
print('Number of metabolites {}'.format(len(model.metabolites)))


allFolders = os.listdir(path)
model_name=model.name.split('tfa_')[-1]
native_exchanges=len(model.exchanges)
#read experimental data checked for feasibility
#all have only 1 alternative apart from kullabergensis
data_file='/remind/projects/bee_project/stats/'

data_dir_cobra= '/remind/models/bee_models/core_members/carveme_models'
model_cobra=cobra.io.read_sbml_model(join(data_dir_cobra,"{}.xml".format(model_name)))
# model_name=allFolders[model_no].split(".xml")[0]
# t=time.time()


"possible carbon sources data"
# frame_screened=pd.read_csv(data_file+'/frame_screened_w_exp_data_stats_feasibility_checked_removed_duplicates_2.csv',index_col=0)
# frame_possible=pd.read_csv(data_file+'/possible_carbon_sources_bee.csv',index_col=None)
# frame_possible=pd.read_csv(data_file+'/possible_carbon_sources_bee_updated.csv',index_col=None)
# frame_possible=pd.read_csv(data_file+'/stats_updated_bee_sources.csv',index_col=None)
# frame_possible=pd.read_csv(data_file+'/possible_carbon_sources_bee_empty_3.csv',index_col=None)
frame_possible=pd.read_csv(data_file+"possible_c_sources_list_300822_initial.csv",index_col=None)


"union of essential metabolites data"
# frame_union=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr.h5')
# frame_union=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr_updated_090822.h5')
# frame_union=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr_updated_110822.h5')
# frame_union=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr_updated_170822_TFA_correct.h5')
# frame_union=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr_updated_180822_TFA.h5')
# frame_union=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr_updated_300822_TFA.h5')
frame_union1=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr_updated_010922_TFA.h5')
frame_union2=pd.read_hdf('/remind/projects/bee_project/constrained_medium/essential_mets/essential_mets_w_cat_repr_updated_061023_2_TFA.h5')
frame_union2_species = frame_union2[frame_union2.model.isin([model_name])]
frame_union= frame_union1.append(frame_union2_species).reset_index()

# frame_union= frame_union1


if model_name=="Lactobacillus_mellifer_Bin4":
    max_consmp = 4

"if we only want to add the essential ones for each model"
if allowed:
    met_pool = list(frame_union.metabolites.unique())
    # remove=['EX_metsox_R__L_e','EX_metsox_S__L_e','EX_idon__L_e']
    # # remove=['EX_metsox_R__L_e', 'EX_metox__R_e']
    # met_pool=[k for k in met_pool if k not in remove]
else:
    frame_met=frame_union[frame_union.model.isin([model_name])]
    met_pool = list(frame_met.metabolites.unique())

# met_pool=list(frame_union.metabolites.unique())
#for all
met_pool=[k for k in met_pool if k!='EX_coa_e']

# list_all=[k for k in frame_screened.mets.unique()]
list_possible=[k for k in frame_possible.annotation_bigg.unique() if k!='no']

# list_sink=[k for k in model.reactions if "sink" in k.id]

diff_rxns=[k for k in model.reactions if k not in model_cobra.reactions]
sink_diff=[k for k in diff_rxns if "sink" in k.id]

model.remove_reactions(sink_diff)
model.repair()





#to close carbon containing sinks
# for reaction in sink_close:
#     reaction.bounds=(0,0)


'''here put list remove'''

# path_fdp='/remind/projects/bee_project/ime_alternatives_180522_ranking_fdps/{}/'.format(model_name)
# file_ = glob.glob(path_fdp + "/*" + "h5")
# df_fdp = pd.read_hdf(file_[0])

'''if fdp approach is done'''
# df_choose=pd.read_hdf('/remind/projects/bee_project/constrained_medium/fdp_stats/all_fdps_for_bee_core.h5')
# list_fdp_all=[]
# for fdp in df_choose.fdp_list:
#     list_fdp_all+=fdp
# list_fdp_all=list(set(list_fdp_all))
# c_sources_to_include=list(set(list_possible+list_fdp_all+met_pool))

c_sources_to_include=list(set(list_possible+met_pool))

# "only for gapicola make this"
#activate this anyways
"""this was for fdp now with thermo firsttry without fdp"""
if model_name=='Gilliamella_apicola_wkB1GF':
    frame_fdp=pd.read_hdf('/remind/projects/bee_project/constrained_medium/fdp_stats/all_fdps_for_bee_gapicola_050822.h5')
    list_fdp=[]
    for fdp in frame_fdp.fdp_list:
        list_fdp+=fdp
    list_fdp=list(set(list_fdp))
    chosen=['EX_23camp_e', 'EX_cmp_e', 'EX_lac__D_e', 'EX_nac_e','EX_etoh_e']
    remove_gapic=[k for k in list_fdp if k not in chosen]
    #latest removal because they alternate
    remove_gapic_2=[ 'EX_csn_e', 'EX_ura_e', 'EX_cytd_e', 'EX_dad_2_e','EX_mmet_e']
    remove_gapic=['EX_nac_e','EX_gln__L_e']+remove_gapic_2
    c_sources_to_include=[k for k in c_sources_to_include if k not in remove_gapic]


if model_name=="Snodgrassella_alvi_wkB2GF":
    not_accept_list=["EX_etoh_e","EX_meoh_e","EX_hxan_e"]
    c_sources_to_include = [k for k in c_sources_to_include if k  not in not_accept_list]


if model_name=='Gilliamella_apicola_wkB1GF':
    not_accept_list=['EX_4hthr_e','EX_pheme_e','EX_etoh_e','EX_uri_e']
    c_sources_to_include = [k for k in c_sources_to_include if k not in not_accept_list]

if model_name=='Lactobacillus_apis_Hma11GF':
    not_accept_list=['EX_phe__L_e','EX_nac_e']
    c_sources_to_include = [k for k in c_sources_to_include if k not in not_accept_list]



# "for fdp removal"
if model_name=='Bifidobacterium_asteroides_PRL2011GF':
    frame_fdp_combined=pd.read_hdf('/remind/projects/bee_project/constrained_medium/fdp_stats/all_fdps_for_bee_150822.h5')
    frame_f=frame_fdp_combined[frame_fdp_combined.model==model_name]
    list_fdp=[]
    for fdp in frame_f.fdp_list:
        list_fdp+=fdp
    list_fdp = list(set(list_fdp))

    chosen=frame_f.fdp_list.iloc[0]
    remove=[k for k in list_fdp if k not in chosen]
    remove_last=[k for k in remove if k not in list_possible]
    print(model_name,remove_last)
    c_sources_to_include = [k for k in c_sources_to_include if k not in remove_last]


if model_name=='Bifidobacterium_asteroides_PRL2011GF':
    not_accept_list=['EX_lys__L_e']
    # not_accept_list=[ ]
    c_sources_to_include = [k for k in c_sources_to_include if k not in not_accept_list]

if model_name=='Lactobacillus_kullabergensis_Biut2GF':
    not_accept_list=['EX_xyl__D_e']
    c_sources_to_include = [k for k in c_sources_to_include if k not in not_accept_list]


if model_name=='Lactobacillus_mellifer_Bin4':
    not_accept_list=['EX_idon__L_e']
    c_sources_to_include = [k for k in c_sources_to_include if k not in not_accept_list]


"remove them"
list_remove=[k for k in listcarbonsources(model) if k not in c_sources_to_include]

#this is to force o2 uptake
if model_name=="Snodgrassella_alvi_wkB2GF":
    model.reactions.EX_o2_e.bounds=(-25,-1)

#fdp chosen for the model
# df_fdp=df_choose[df_choose.model.isin([model_name])]
# #sum all metabolites
# list_fdp=[]
# for fdp in df_fdp.fdp_list:
#     list_fdp+=fdp
# list_fdp=list(set(list_fdp))
# selected_fdp=[k for k in df_fdp.iloc[0].fdp_list]
# list_remove+=[k for k in list_fdp if k not in selected_fdp]
# if model_name=='Lactobacillus_kullabergensis_Biut2GF':
#     list_remove+=['EX_glu__L_e']


# exchanges = [k.id.replace('EX_', '') for k in model.exchanges]
# not_common = [k for k in exchanges if k not in list_all]
# #remove the exchanges that are not of interest
# rxns_to_rmv = ['EX_{}'.format(k) for k in not_common]
'add the not selected fdp metabolites'

rxns_to_rmv=list_remove

model.remove_reactions(rxns_to_rmv)
after_exchanges=len(model.exchanges)

list_not_main_c_source=[model.reactions.get_by_id(k) for k in listcarbonsources(model) if k in met_pool]

print("For model {} native exchanges is {} after screening is {}".format(model.name,native_exchanges,after_exchanges))


max_flux=25

if constrain_fluxes:
    for k in model.exchanges:
        if k.id!='EX_o2_e':
            k.bounds=(-25,max_flux)
    for r in  list_not_main_c_source:
        r.bounds=(-1*max_flux_lim,max_flux)

#new for inorganics
list_carbon_sources=listcarbonsources(model)
inorganics=[k for k in model.exchanges if k.id not in list_carbon_sources]
inorganics=[k for k in inorganics if k.id!="EX_o2_e"]

for rxn_inor in inorganics:
    rxn_inor.bounds=(-25,50)

if 'EX_glc__D_e' in model.reactions:
    model.reactions.get_by_id('EX_glc__D_e').upper_bound=0

if 'EX_fru_e' in model.reactions:
    model.reactions.get_by_id('EX_fru_e').upper_bound=0

if 'EX_cit_e' in model.reactions:
    model.reactions.get_by_id('EX_cit_e').upper_bound=0

if 'EX_icit_e' in model.reactions:
    model.reactions.get_by_id('EX_icit_e').upper_bound=0

#new additional info
"it is known that it cannot catabolize carbohydrates"
if model_name=="Snodgrassella_alvi_wkB2GF":
    model.reactions.get_by_id('EX_fru_e').lower_bound = 0
    model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0


# Solver settings
def apply_solver_settings(model, solver = solver):
    model.solver = solver
    # model.solver.configuration.verbosity = 1
    model.solver.problem.parameters.reset()
    model.solver.configuration.tolerances.feasibility = 1e-9
    model.solver.configuration.tolerances.optimality = 1e-9
    model.solver.configuration.tolerances.integrality = 1e-9

    if solver == 'optlang_gurobi':
        model.solver.problem.Params.NumericFocus = 3
    model.solver.configuration.presolve = True

apply_solver_settings(model)


limit=growth_limit
model.reactions.get_by_id(biomass).upper_bound=limit
model.reactions.get_by_id(biomass).lower_bound=limit*0.98

#force citrate for snod


"""catabolite repression decide to put 0 or 1"""
list_carbon_sources = listcarbonsources(model)
list_catabolite=[k for k in list_possible if k in list_carbon_sources ]
constrain_sources(model,list_catabolite, max_cnsmp=max_consmp, min_cnsmp=0)


#for sink reaction minimization


#
'2) get list of rxn names for parsimonious uptakes'
if (parsimonious=='c_uptakes'):
    list_carbon_sources=listcarbonsources(model)
    '3)apply parsimonious and get solution and expression'
    'minimize uptake of carbon sources or all uptakes'
    expr_pars,sol_pars=minimizecarbonuptakes(model, list_carbon_sources)
 #todo write another function for C molar balance
if (parsimonious=='c_uptakes_only_C_moles'):
    list_carbon_sources=listcarbonsources(model)
    '3)apply parsimonious and get solution and expression'
    'minimize uptake of carbon sources or all uptakes'
    expr_pars,sol_pars=minimizecarbonuptakes_cmoles_only(model, list_carbon_sources)
if (parsimonious=='c_uptakes_mw'):
    list_carbon_sources=listcarbonsources(model)
    '3)apply parsimonious and get solution and expression'
    'minimize uptake of carbon sources or all uptakes'
    expr_pars,sol_pars=minimizecarbonuptakes_by_weight(model, list_carbon_sources,molecular_weight='formula_weight')
if (parsimonious=='all_uptakes_mw'):
    #list_carbon_sources=listcarbonsources(mytfa)
    '3)apply parsimonious and get solution and expression'
    'minimize uptake of carbon sources or all uptakes'
    # rxn_list=[model.reactions.get_by_id(r) for r in listcarbonsources(model) if r in list_possible]
    expr_pars,sol_pars=minimizealluptakes_by_weight(model,model.exchanges, molecular_weight='formula_weight')
elif (parsimonious=='all_uptakes'):
    expr_pars,sol_pars=minimizealluptakes(model)

print(sol_pars)
'4)add constrain for parsimonious'
'minimize uptake of carbon sources or all uptakes user defined'
cons0 = model.add_constraint(ModelConstraint,
                              model,
                              expr_pars,
                              id_='pars_constraint',
                                #fix uptake
                              lb=(sol_pars.objective_value)*tolerance_lb,
                              ub=(sol_pars.objective_value)*(tolerance))


model.repair()



# rxns_to_consider=[model.reactions.get_by_id(k) for k in list_carbon_sources if k not in list_possible]
rxns_to_consider=[model.reactions.get_by_id(k) for k in list_carbon_sources ]

# if model_name=="Snodgrassella_alvi_wkB2GF":
#     rxns_to_consider+=[model.reactions.EX_o2_e]

# list_exchanges=[k.id for k in mytfa2.exchanges]
# list_c_reverse= [mytfa2.reactions.get_by_id(k).reverse_variable.name for k in list_carbon_sources]

'generate untill max alternative'
#mytfa2.exchanges

# my_list,df_all=imeanalysis_all_size(model,model.exchanges,getAlternatives=True,max_alternative=max_alternative,biomass=biomass,stop_cond=True,criteria=2)
my_list,df_all=imeanalysis_all_size(model,rxns_to_consider,getAlternatives=True,max_alternative=max_alternative,biomass=biomass,stop_cond=True,criteria=2)

'generate untill max alternative'
#my_list,df_all=imeanalysis_all_size_o2(mytfa2,mytfa2.exchanges,getAlternatives=True,max_alternative=10)

frame = pd.concat(my_list, ignore_index=True)
# frame['uptake_min']=np.ones(frame.shape[0])*sol_pars.objective_value
frame['yield_perc']=np.ones(frame.shape[0])*1.0/tolerance_lb

alternative=frame.alternative.max()
elapsed = time.time() - t
print('time for ime analysis', elapsed)

from remind.utils.postprocessing import *
add_size(frame,['alternative'],'alt_size')

# from remind.core.decomposition import find_essential_alternates_from_dime_list
# a,b=find_essential_alternates_from_dime_list(frame)


#output file to save the results
output_file_all='/remind/projects/bee_project/DIMES_correct_TFA_allowed_{}_limited_carbon/{}/all_vars'.format(allowed,model.name.split('tfa_')[-1])
if not os.path.exists(output_file_all):
    os.makedirs(output_file_all)

output_file_alts='/remind/projects/bee_project/DIMES_correct_TFA_allowed_{}_limited_carbon/{}/alternatives'.format(allowed,model.name.split('tfa_')[-1])
if not os.path.exists(output_file_alts):
    os.makedirs(output_file_alts)



#"store them in hdf5 files"
frame.to_hdf(output_file_alts+'/alternative_mets_for_{}_growth_{}_alt_{}_pars_{}_tol_pars_{}_tol_lb_{}.h5'.format(model.id,limit,alternative,parsimonious,tolerance,tolerance_lb),key='s')
df_all.to_hdf(output_file_all+'/alternatives_all_vars_for_{}_growth_{}_alt_{}_pars_{}_tol_pars_{}_tol_lb_{}.h5'.format(model.id,limit,alternative,parsimonious,tolerance,tolerance_lb),key='s')


# d1=pd.read_hdf('/remind/projects/bee_project/ime_alternatives_180522_remove_inorganics/Gilliamella_apicola_wkB1GF/alternatives/alternative_mets_for_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_5000.0_pars_all_uptakes_mw_tol_pars_1.11000001_tol_lb_1.11.h5')
# d2=pd.read_hdf('/remind/projects/bee_project/ime_alternatives_180522_remove_inorganics/Gilliamella_apicola_wkB1GF/alternatives/alternative_mets_for_Gilliamella_apicola_wkB1_xml_growth_0.2_alt_5517.0_pars_all_uptakes_mw_tol_pars_1.11000001_tol_lb_1.11.h5')
#

