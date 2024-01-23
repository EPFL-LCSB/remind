

# from pytfa.thermo.tmodel_struct import ThermoModelStructure
import os
from os.path import join
from pytfa.io import import_matlab_model
from pytfa.optim.constraints import ModelConstraint
from remind.core.exchanges import imeanalysis
from remind.core.medium import constrain_sources
from remind.core.parsimonious import minimize_element_upt, minimizealluptakes_by_weight
import time
import pandas as pd
# import numpy as np
import pytfa
from def_medium_toy import c_sources


CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = GUROBI

analysis_type = 'DIME'

organism = 'Geobacter'
# organism = 'Rhodoferax'
# organism = 'Shewanella'

growth_ids = {
    'Geobacter' : 'agg_GS13m',
    'Rhodoferax': 'BIO_Rfer3',
    'Shewanella': 'Biomass',
    }

GROWTH = growth_ids[organism]

output_file='./ime_alternatives/imes_toy'
#output_file='./milp_3step_011120'
if not os.path.exists(output_file):
    os.makedirs(output_file)

data_dir='../models/Toy_community'
model_name = '{}.mat'.format(organism)

raw_model = import_matlab_model(join(data_dir, model_name))
raw_model.name = organism
# For an unknown reason St is not good ID for a reaction
try:
    raw_model.reactions.St.id = 'STTT'
except AttributeError:
    pass


t=time.time()


# 'here thermomodel is created without thermo constraints but FU BU constraints'
# 'Another class called ThermoModelStructure is created (ask me if required inside pytfa so not pushed)'
# mytfa =ThermoModelStructure ([], salvi)
thermo_data = {'name':'', 'units':'', 'metabolites':{}, 'cues':{}}
mytfa = pytfa.ThermoModel(thermo_data, raw_model)

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

### Fix the growth
growth_limit = 0.2 # to set it to a biologically more releveant value later
# make sure all the exchanges are open
for rxn in mytfa.exchanges:
    rxn.bounds = (0,50)
    
# define the medium
medium = pd.read_csv('ime_alternatives/imes_toy/medium.csv',index_col=0).to_dict(orient='index')
c_source_exchange = []
for met_id in medium.keys():
    try:
        met = mytfa.metabolites.get_by_id(met_id)
        exchange = [rxn for rxn in met.reactions if len(rxn.metabolites)==1]
        if len(exchange) == 0: # no uptakes of this nutrient
            continue
        elif len(exchange) > 1:
            raise Exception('There might be two exchange reactions associated with this metabolite. '
                            'Please remove the redundant exchanges.')
        else:
            exchange = exchange[0]
            # open the reaction
            exchange.bounds = (-50,50)
            if met_id in c_sources:
                c_source_exchange += [exchange.id]
    
    except KeyError: # the model does not have this metabolite
        pass

 
### finding the maximum yield
tol_ub = 1
tol_lb = 1
# Presolve and find alternating groups (only for carbon and oxygen)
mytfa.reactions.get_by_id(GROWTH).lower_bound = growth_limit
mytfa.reactions.get_by_id(GROWTH).upper_bound = growth_limit

# find the minimal uptake by taking parsimonious uptake
# sometimes it is better to couple and sometimes it is better to decouple
'3)apply parsimonious and get solution and expression'
'minimize uptake of carbon sources or all uptakes'
constrain_sources(mytfa, c_source_exchange)
expr_pars, sol_pars=minimizealluptakes_by_weight(mytfa, molecular_weight='formula_weight')
min_uptake = sol_pars.objective_value
mytfa.add_constraint(ModelConstraint,
                            mytfa,
                            expr_pars,
                            id_='pars_constraint',
                              #fix uptake
                            lb = tol_lb*min_uptake,
                            ub = tol_ub*min_uptake)

    

my_list, _ = imeanalysis(mytfa, mytfa.exchanges, 
                                 # secretions=the_rest,
                                  element='all',
                                  enforce_size=False,
                                  getAlternatives=True, 
                                  max_alternative=1000)
    



frame = pd.concat(my_list, ignore_index=True)
frame['yield'] = len(frame) * [1/tol_ub]
alternative=frame.alternative.max()


frame.to_csv(output_file+'/{}_{}_alternative_mets_for_growth_{}_{}.csv'.format(
    analysis_type, organism, growth_limit, tol_ub))
# df_all.to_csv(output_file+'/{}_{}_alternatives_all_vars_for_growth_{}.csv'.format(
#     analysis_type, model_name, growth_limit))
elapsed = time.time() - t
print('time for ime analysis', elapsed)

