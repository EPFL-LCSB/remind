#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:03:46 2022

@author: ReMIND Team
"""
import pandas as pd

from os.path import join

from pytfa.io import import_matlab_model

import pytfa

from remind.core.medium import biotic_mets


CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = GUROBI

organisms = [
            'Geobacter',
            'Rhodoferax',
            # 'Shewanella',
            ]

growth_ids = {
    'Geobacter' : 'agg_GS13m',
    'Rhodoferax': 'BIO_Rfer3',
    'Shewanella': 'Biomass',
    }


def open_exchanges(model, organism):
    
    if organism == 'Geobacter':
        allowed = ['h_e',
             # 'fe2_e',
             # 'co2_e',
             'so4_e',
             'pi_e',
             'nh4_e',
             'mg2_e',
             'k_e',
             'h2o_e',
             'fe3_e',
             # 'ac_e',
             'zn2_e',
             # 'succ_e',
             'ss_e',
             's_e',
             # 'oxa_e',
             'ni2_e',
             'na1_e',
             'n2_e',
             'mobd_e',
             'mn2_e',
             # 'mal-L_e',
             # 'lac-L_e',
             'h2s_e',
             'h2_e',
             # 'fum_e',
             # 'for_e',
             'cu2_e',
             'cobalt2_e',
             'cl_e',
             # 'cit_e',
             'cd2_e',
             'ca2_e']
    elif organism == 'Rhodoferax':
        allowed = [
            # '2ddglcn_e',
             # 'ac_e',
             # 'arab-L_e',
             # 'buts_e',
             # 'bz_e',
             'ca2_e',
             'cd2_e',
             # 'cellb_e',
             # 'cit_e',
             # 'co2_e',
             'cobalt2_e',
             'cro4_e',
             'cu2_e',
             # 'eths_e',
             # 'etoh_e',
             # 'fe2_e',
             'fe3_e',
             # 'fru_e',
             # 'fum_e',
             # 'glc_e',
             # 'glyclt_e',
             'h_e',
             'h2_e',
             'h2o_e',
             # 'hexs_e',
             # 'istnt_e',
             'k_e',
             # 'lac-L_e',
             # 'mal-L_e',
             'mg2_e',
             'mn2_e',
             'mobd_e',
             'na1_e',
             'nh4_e',
             'ni2_e',
             # 'no2_e',
             # 'no3_e',
             # 'o2_e',
             'pi_e',
             # 'ppa_e',
             # # 'ppn_e',
             # 'pyr_e',
             # 'rib-D_e',
             'so4_e',
             # 'succ_e',
             # 'sula_e',
             'tsul_e',
             'zn2_e',
            ]
    elif organism == "Shewanella":
        allowed = ['ac_e',
                # 'akg_e',
                'arsna_e',
                'arsni2_e',
                'ca2_e',
                'cl_e',
                # 'co2_e',
                'cobalt2_e',
                'cobalt3_e',
                'CrOH3_e',
                'cro4_e',
                'cu2_e',
                # 'dms_e',
                # 'dmso_e',
                # 'etoh_e',
                # 'fe2_e',
                'fe3_e',
                # 'for_e',
                # 'fum_e',
                # 'glyclt_e',
                # 'glyc-R_e',
                'h_e',
                'h2_e',
                'h2o_e',
                'h2o2_e',
                'h2s_e',
                # 'hdca_e',
                'hg2_e',
                # 'inoshp_e',
                # 'inospp1_e',
                'k_e',
                # # 'lac-D_e',
                # 'lac-L_e',
                # 'mal-L_e',
                'mg2_e',
                'mn2_e',
                'mn4o_e',
                'mobd_e',
                'na1_e',
                'nh4_e',
                'ni2_e',
                # 'no2_e',
                # 'no3_e',
                # 'o2_e',
                # 'ocdca_e',
                'pi_e',
                # 'ppa_e',
                # 'pyr_e',
                'so3_e',
                'so4_e',
                # 'succ_e',
                'tsul_e',
                # 'ttdca_e',
                'tttnt_e',
                'wo4_e',]

    # we close demand reactions here, as the other functions do not care about the demands
    # close the demand reactions in Geobacter
    if organism == "Geobacter":
        for rxn in model.reactions:
            if 'DM_' in rxn.id: # demand reaction?
                rxn.bounds = (0,0)
                   
                
    return allowed


c_sources = ['ac_e', 'succ_e', 'oxa_e', 'mal-L_e', 'lac-L_e', 'cit_e', 'fum_e',
             'for_e', 'arab-L_e', 'etoh_e', 'eths_e', 'fru_e', 'glc_e', 'hexs_e',
             'glyclt_e', 'istnt_e', 'rib-D_e', 'buts_e', 'ppa_e', 'pyr_e', 'sula_e',
             'dmso_e', 'dms_e', 'glyc-R_e', 'akg_e']  # carbon sources are defined


if __name__ == "__main__":
    data_dir='../models/Toy_community'
    model_list = []
    essential_nutrients = [] # we take the union
    for organism in organisms:
        model_name = '{}.mat'.format(organism)
        
        raw_model = import_matlab_model(join(data_dir, model_name))
        raw_model.name = organism
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
        # opening the essentials other than carbon
        these_mets = open_exchanges(mytfa, organism)
        for met in these_mets:
            if met not in essential_nutrients: # a unique union of all
                essential_nutrients.append(met) 
            
        model_list += [mytfa]
    
    # if we assume catabolite repression exists it should be considered at this stage too
    # so we add carbon sources 1 by 1
    total_biotics = []
    for c_source in c_sources:
        abiotics = essential_nutrients + [c_source]
        biotics,_ = biotic_mets(model_list, growth_ids, abiotics)
        total_biotics += [x for x in biotics if x not in total_biotics]

    medium = {}
    for met in (essential_nutrients + c_sources + total_biotics):
        medium[met] = {}
        if met in (essential_nutrients+c_sources):
            medium[met]['type'] = 'abiotic'
        else:
            medium[met]['type'] = 'biotic'
        
    data = pd.DataFrame.from_dict(medium, orient='index')
    data.to_csv('ime_alternatives/imes_toy/medium.csv')