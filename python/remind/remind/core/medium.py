#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 08:46:58 2022

@author: ReMIND Team
"""
from warnings import warn
from copy import deepcopy
import pandas as pd
from pytfa.optim.variables import BackwardUseVariable, ForwardBackwardUseVariable
from pytfa.optim.constraints import ModelConstraint
from pytfa.optim.utils import symbol_sum
from .parsimonious import minimizealluptakes_by_weight
from .exchanges import imeanalysis, ActiveReaction
from ..utils.postprocessing import get_subs_prods_with_id
from ..optim.constraints import AbioticResourceConstraint,AbioticResourceConstraint2,AbioticResourceConstraint3
from ..optim.variables import AbioticResource, PositiveInteraction, UptakeActivation
import numpy as np
from optlang.exceptions import SolverError
from optlang.interface import INFEASIBLE
from .parsimonious import get_C_number_from_an_exch_rxn
from ..optim.variables import MaximalFluxActiveVariable, MaximalFluxActiveVariableLinearize, \
    MaximalFluxActiveVariableSum
from ..optim.constraints import LinearizeInteger_1, LinearizeInteger_2, LinearizeInteger_3, \
    ConstraintBond, ConstraintEssential, ConstraintMFA, ConstraintMaxFlux, ConstraintSum

big_M=100


def biotic_mets(species_models, growth_ids, abiotics, num_iter=100,
                min_growth=0.1, min_yield_cut=0.1, growth=0.2):
    '''
    An iterative procedure to find all potential byproducts and make them available
    in the next round. The idea is to open these exchanges for DiMEs.

    Parameters
    ----------
    species_models : list of pytfa.ThermoModel
        DESCRIPTION.
    growth_ids : dict
        Keys are model names and values are growth ids.
    abiotics : list
        list of metabolite ids that we expect to be present in the medium.
    num_iter : int, optional
        maximum number of iterations. The default is 100.
    min_growth : float, optional
        fraction of the maximum growth by which we can say an organism can grow. The default is 0.1.
    min_yield_cut : float, optional
        The minimum yield cut that you expect to explore. The default is 0.1.
    growth : float, optional
        The value of growth to be fixed. The default is 0.2.

    Raises
    ------
    Exception
        DESCRIPTION.

    Yields
    ------
    biotics : list
        list of metabolite IDs.
    nongrowing : list
        list of the model names.

    '''
    
    iteration = 0
    converged = False
    
    
    # first open all available abiotic metabolites
    for model in species_models:
        # make sure all uptakes are closed and secretions are open in the beginning
        for rxn in model.exchanges:
            rxn.bounds = (0,1000)
            
        for met_id in abiotics:
            try:
                met = model.metabolites.get_by_id(met_id)
                exchange = [rxn for rxn in met.reactions if len(rxn.metabolites)==1]
                if len(exchange) == 0: # no uptakes of this nutrient
                    continue
                elif len(exchange) > 1:
                    raise Exception('There might be two exchange reactions associated with this metabolite. '
                                    'Please remove the redundant exchanges.')
                else:
                    exchange = exchange[0]
                    # open the reaction
                    exchange.bounds = (-1000,1000)
            
            except KeyError: # the model does not have this metabolite
                continue
            
            
                    
                
    biotics = []
    nongrowing = [model.name for model in species_models] # first assume nobody grows
    while iteration < num_iter and not converged:
        num_current_biot = len(biotics)
        # In each iteration we need to recalculate parsimonious as it might change
        # subject to the available nutrients, the ID of byproducts is saved -> not important which one grows
        for species in species_models:
            model = deepcopy(species)
            growth_id = growth_ids[model.name]
            growth = model.slim_optimize()
            if growth > min_growth: # it grows!
                if model in nongrowing:
                    nongrowing.remove(model)
                
                # if parsimonious constraint was already added, it should be removed
                try:
                    model.remove_constraint(model._cons_dict['MODC_pars_constraint'])
                except KeyError:
                    pass
                
                # we do dimes, but not decomposed based on the yields
                # Since one carbon source is active at anytime, we do not expect co-utilization so low yield solution should not cut high yield solutions
                model.reactions.get_by_id(growth_id).lower_bound = growth
                model.reactions.get_by_id(growth_id).upper_bound = growth
                expr_pars, sol_pars = minimizealluptakes_by_weight(model, molecular_weight='formula_weight')
                min_uptake = sol_pars.objective_value
                model.add_constraint(ModelConstraint,
                                            model,
                                            expr_pars,
                                            id_='pars_constraint',
                                              #fix uptake
                                            lb = min_uptake,
                                            ub = 1/min_yield_cut*min_uptake)
                my_list, _ = imeanalysis(model, model.exchanges, 
                                 # secretions=the_rest,
                                  element='all',
                                  enforce_size=False,
                                  getAlternatives=True, 
                                  max_alternative=100)
                frame = pd.concat(my_list, ignore_index=True)
                _, prods = get_subs_prods_with_id(frame)
                biotics += [m for m in prods if m not in biotics]
            else:
                pass
            
        # check if the length of biotic nutrients not chnaged, then we have converged
        if len(biotics) == num_current_biot:
            converged = True
        
        iteration +=1
        
    if iteration == num_iter:
        warn('The iterations to specify the biotic metabolites did not converge. '
             'Consider increasing the number of iterations.')
    
    return biotics, nongrowing

def constrain_sources(model, exchange_list, max_cnsmp=1, min_cnsmp=1):
    '''
    Constrain the number of C sources that are used at any simulation to be 1, i.e.,
    catabolite repression

    Parameters
    ----------
    model : TYPE
        DESCRIPTION.
    exchange_list : TYPE
        DESCRIPTION.
    max_cnsmp : TYPE, optional
        DESCRIPTION. The default is 1. use at most 1
    min_cnsmp : TYPE, optional
        DESCRIPTION. The default is 1. use at least one

    Returns
    -------
    None.

    '''
    
    # This function might be called before parsimonious or before DiMEs (even before biotic_mets)
    # Since it relies on adding BFUSE variables we need to make sure they exist in the model
    buse_vars = model.get_variables_of_type(BackwardUseVariable)
    
    
     
    constrain_expr = symbol_sum([buse_vars.get_by_id(id_) for id_ in exchange_list])
    model.add_constraint(ModelConstraint,
                        model,
                        constrain_expr,
                        id_='catabolite_repression',
                          #fix uptake
                        lb = min_cnsmp,
                        ub = max_cnsmp)
    
    return




def constrain_abiotics(model, metabolite_list = []):
    '''
    This fundtion to minimize the number of abiotic nutrient uptakes out of a given 
    list metabolite_list. This can be used to minimize abiotic nutrients in a community
    with more than two species or to minimize the uptake of new nutrients in n±1 analysis.

    Parameters
    ----------
    model : InteractionModel
        The community model.
    metabolite_list : list, optional
        The list of metabolites whose uptake we constrained. The default is [].

    Returns
    -------
    None.

    '''

    if len(metabolite_list) == 0:
        metabolite_list = model.metabolites
    
    
    pos_int = model.get_variables_of_type(PositiveInteraction)
    uptakes = model.get_variables_of_type(UptakeActivation)
    
    obj_vars = []
    for met in metabolite_list:
            
        this_met_upts = [var for var in uptakes if var.reaction.metabolite == met]
        if len(this_met_upts) == 0: # this metabolite is not uptaken by any
            continue # do nothing!
        
        nn_var = model.add_variable(AbioticResource,
                            met,
                            lb = 0,
                            ub = 1,
                            )
        obj_vars += [nn_var]
        
        this_pi_var = [var for var in pos_int if var.metabolite == met]

        expr=nn_var-symbol_sum([upt_var for upt_var in this_met_upts])
        model.add_constraint(AbioticResourceConstraint3,
                             hook=met,
                             expr=expr,
                             ub=0)

        if len(this_pi_var) == 0: # no positive interaction is possible for this metabolite
            for upt_var in this_met_upts:
                rxn = upt_var.reaction
                expr = upt_var - nn_var
                model.add_constraint(AbioticResourceConstraint,
                                    hook = rxn,
                                    expr = expr,
                                    ub = 0)
        else:
            this_pi_var = this_pi_var[0]
            for upt_var in this_met_upts:
                rxn = upt_var.reaction
                expr = upt_var - this_pi_var - nn_var
                model.add_constraint(AbioticResourceConstraint,
                                    hook = rxn,
                                    expr = expr,
                                    ub = 0)

                expr2=nn_var+this_pi_var

                model.add_constraint(AbioticResourceConstraint2,
                                    hook = rxn,
                                    expr = expr2,
                                    ub = 1.0)

    model.objective = symbol_sum(obj_vars)
    model.objective_direction = 'min'
    #if not repair it wont find the variable afterwards
    model.repair()
    
    return

"""this functionto be deleted later just to make sure i keep for now"""
def constrain_uptake_old(model,reactions,getAlternatives=True,max_alternative=500,stop_cond=False,block_essentials=True,essential_list=[],criteria=0,biomass='Growth',max_flux_lim=3.0):
    """
    reactions:list of reactions for the constraint to be added
    block essentials : if you dont want to include certain substances
    as a possible C source bcs maybe they are essential anyways prevent them
    by providing a list
    """
    #constrain to make
    max_flux_non_c_source = max_flux_lim
    for rxn in reactions:

        # print(get_C_number_from_an_exch_rxn(rxn))
        # C_NUM=get_C_number_from_an_exch_rxn(rxn)
        max_flux_non_c_source=max_flux_lim
        FU = model.forward_use_variable.get_by_id(rxn.id)
        BU = model.backward_use_variable.get_by_id(rxn.id)
        BFUSE = model.add_variable(ForwardBackwardUseVariable, rxn)
        # BFUSE=model._var_dict['BFUSE_{}'.format(rxn.id)]

        F_ = model.reactions.get_by_id(rxn.id).forward_variable
        B_ = model.reactions.get_by_id(rxn.id).reverse_variable
        MFA = model.add_variable(MaximalFluxActiveVariable,rxn)
        MFAL = model.add_variable(MaximalFluxActiveVariableLinearize, rxn)
        MFAS = model.add_variable(MaximalFluxActiveVariableSum, rxn)

        #first add linearization =>MFAL=MFA*BU

        expression_1 = MFAL-MFA
        expression_2 = MFAL-BU
        expression_3 = MFAL-BU-MFA+1.0
        'if bfuse=1 BU=0 why we have this C --> ? from redGEMx '
        'not required'  # todo remove later
        expression=MFA-BU
        cons0 = model.add_constraint(LinearizeInteger_1,
                                     rxn,
                                     expression_1,
                                     ub = 0)

        cons1 = model.add_constraint(LinearizeInteger_2,
                                     rxn,
                                     expression_2,
                                     ub = 0)


        cons3 = model.add_constraint(LinearizeInteger_3,
                                     rxn,
                                     expression_3,
                                     lb = 0)

        #mfa<=BU if bu==0 mfa ==0
        cons4 = model.add_constraint(ConstraintMFA,
                                 rxn,
                                 expression,
                                 ub=0)

        expr_max_flux=B_-max_flux_non_c_source*MFAL-(1-MFAL)*big_M


        cons5 = model.add_constraint(ConstraintMaxFlux,
                                 rxn,
                                 expr_max_flux,
                                 ub=0)

        expr_new_var=MFAS-(BU+MFA)


        cons6 = model.add_constraint(ConstraintSum,
                                 rxn,
                                 expr_new_var,
                                 lb=0,ub=0)


        #bondinf constraint you have FU bcs it indicates that
        #with this main c ssource uptake there is also little secretion
        #indicating that it doesnt uptake all the little C sources at the same time and secretes
        #products of these little c sources if there is another c source that can do the job w less products
        #if interested the other can be run too
        #then the essential_met_pool will be higher resulting in more alternatives and more little c sources playing
        #with each other
        #1
        expr_bonding=BU-MFAL+BFUSE+FU
        #2
        # expr_bonding=BU-MFAL+BFUSE

        cons7 = model.add_constraint(ConstraintBond,
                                 rxn,
                                 expr_bonding,
                                 lb=0,ub=1.1)
    model.repair()
    #block essentials if you dont want to add essential substrate to a main C source list
    #as this can cause infeasibility of the other models and also lead to utilization of
    #this essential substrate only in big fluxes due to the catabolite repression constraint

    if block_essentials:
        for r in essential_list:
            rxn=model.reactions.get_by_id(r)
            BU_var=model._var_dict['BU_{}'.format(r)]
            MFA_var=model._var_dict['MFA_{}'.format(r)]
            expr_essential=(BU_var+MFA_var)-2.0*BU_var
            cons_ess = model.add_constraint(ConstraintEssential,
                                     rxn,
                                     expr_essential,
                                     lb=0)



        # mytfa_imm = imeAnalysis(mytfa, mytfa.exchanges)
    MFA_VARS = model._var_kinds[MaximalFluxActiveVariable.__name__]
    MFA_var_names = [x.name for x in MFA_VARS]

    MFAL_VARS = model._var_kinds[MaximalFluxActiveVariableLinearize.__name__]
    MFAL_var_names = [x.name for x in MFAL_VARS]

    MFAS_VARS = model._var_kinds[MaximalFluxActiveVariableSum.__name__]
    MFAS_var_names = [x.name for x in MFAS_VARS]

    BFUSE_VARS = model._var_kinds[ForwardBackwardUseVariable.__name__]
    BFUSE_var_names = [x.name for x in BFUSE_VARS]

    expression = sum(BFUSE_VARS) #before sum(MFA/MFAL)
    model.objective = expression
    model.objective_direction = 'max'
    sol = model.optimize()
    list_ = []
    df_all = pd.DataFrame()
    df_all['Alternative_1'] = np.zeros(len(model.variables))
    df_all.index = sol.raw.index

    df_sol = sol.raw
    df_flux = sol.fluxes
    df_bfuse = df_sol.loc[BFUSE_var_names]
    idx = df_bfuse[df_bfuse < 0.5].index
    active_vars = [model._var_dict[k] for k in idx]
    active_rxn_id = [i.reaction.id for i in active_vars]
    min_size = len(active_rxn_id)


    df_mfas=df_sol.loc[MFAS_var_names]
        # todo check integer feasibility <1 but close to 1 or no !
    idx_c_source = df_mfas[(df_mfas < 1.5) &(df_mfas > 0.5)].index
    active_vars_c_source=[model._var_dict[k] for k in idx_c_source]
    active_rxn_id_c_source = [i.reaction.id for i in active_vars_c_source]

    min_size_c=len(active_rxn_id_c_source)
    print('MIN SIZE high flux c_sources is {}'.format(min_size_c))

    'from here to be changed'
    if getAlternatives:
        # add integer cut constraints
        alternative = 1
        while ((model.solver.status != INFEASIBLE) and (alternative <= max_alternative)):

            # todo look for primal from lUMPGEM
            df_sol = sol.raw
            df_flux = sol.fluxes
            df_bfuse = df_sol.loc[BFUSE_var_names]
            # todo check integer feasibility <1 but close to 1 or no !
            idx = df_bfuse[df_bfuse < 0.5].index
            active_vars = [model._var_dict[k] for k in idx]
            df = pd.DataFrame(columns=['metabolites', 'flux', 'FU', 'BU'])
            active_rxn_id = [i.reaction.id for i in active_vars]

            df_mfas = df_sol.loc[MFAS_var_names]
            # todo check integer feasibility <1 but close to 1 or no !
            #check the active ones with MFAS variable
            idx_c_source = df_mfas[(df_mfas < 1.5) & (df_mfas > 0.5)].index
            active_vars_c_source_mfas = [model._var_dict[k] for k in idx_c_source]
            active_rxn_id_c_source_mfas = [i.reaction.id for i in active_vars_c_source_mfas]
            active_vars_c_source = [k for k in model._var_kinds[MaximalFluxActiveVariable.__name__] if k.id in active_rxn_id_c_source_mfas]
            active_vars_c_source_bu = [k for k in model._var_kinds[BackwardUseVariable.__name__] if k.id in active_rxn_id_c_source_mfas]
            print('control is {} '.format(len(active_rxn_id_c_source_mfas)))

            if stop_cond:
                if (abs(len(active_rxn_id_c_source_mfas) - min_size_c) > criteria):
                    print(
                        'Stopping criteria is reached at alternative {} min+{} alternatives are exhausted for model'.format(
                            alternative, criteria, model.id))
                    break
            fu_vars = [model.forward_use_variable.get_by_id(k).name for k in active_rxn_id]
            bu_vars = [model.backward_use_variable.get_by_id(k).name for k in active_rxn_id]

            mfa_vars=[k.name for k in model._var_kinds[MaximalFluxActiveVariable.__name__] if k.id in active_rxn_id]
            mfal_vars=[k.name for k in model._var_kinds[MaximalFluxActiveVariableLinearize.__name__] if k.id in active_rxn_id]
            mfas_vars=[k.name for k in model._var_kinds[MaximalFluxActiveVariableSum.__name__] if k.id in active_rxn_id]

            df['metabolites'] = active_rxn_id
            df['BFUSE'] = df_bfuse[df_bfuse < 0.5].values
            'error was here this var is only for forward take net flux'
            df['flux'] = df_flux[active_rxn_id].values
            df['FU'] = df_sol[fu_vars].values
            df['BU'] = df_sol[bu_vars].values
            df['MFA']=df_sol[mfa_vars].values
            df['MFAL']=df_sol[mfal_vars].values
            df['MFAS']=df_sol[mfas_vars].values
            df['alternative'] = np.ones(len(active_rxn_id)) * alternative
            # add growth also as a column
            df['Growth'] = np.ones(len(active_rxn_id)) * sol.fluxes[biomass]

            df_all['Alternative_{}'.format(alternative)] = df_sol.values

            # expr = [k.reverse_variable.name for k in tmodel.exchanges]
            # df['uptake_yield']=np.ones(len(active_rxn_id)) * df_all.loc[expr].sum()
            df['uptake_yield'] = np.ones(len(active_rxn_id)) * abs(df[df.BU > 0.5].flux.sum())

            print('Generating alternative {} of size {} for model {}'.format(alternative, len(active_rxn_id),
                                                                             model.id))
            list_.append(df)
            #try to break
            if len(active_rxn_id_c_source)<1:
                break
            expression = symbol_sum(BFUSE_VARS)

            #add integer cuts based on the main c source
            expression2 = sum(active_vars_c_source_bu)
            cons1 = model.add_constraint(ModelConstraint,
                                          model,
                                          expression2,
                                          id_='max_flux_constraint_integer_cut_{}'.format(alternative),
                                          ub=0)
            model.repair()

            try:
                # tmodel.objective= sympy.S.Zero
                sol = model.optimize()
            except SolverError:
                print('Number of alternatives generated is ', alternative)
                break

            alternative += 1


    return list_, df_all

    return



def constrain_uptake(model,reactions,getAlternatives=True,max_alternative=500,stop_cond=False,block_essentials=True,essential_list=[],criteria=0,biomass='Growth',max_flux_lim=3.0):
    """
    This is the updated and clearer version other to be deleted later
    linearization for the two integers were not necessary which makes the problem simpler
    reactions:list of reactions for the constraint to be added
    block essentials : if you dont want to include certain substances
    as a possible C source bcs maybe they are essential anyways prevent them
    by providing a list
    """
    #constrain to make
    max_flux_non_c_source = max_flux_lim
    for rxn in reactions:

        # print(get_C_number_from_an_exch_rxn(rxn))
        C_NUM=get_C_number_from_an_exch_rxn(rxn)
        max_flux_non_c_source=max_flux_lim
        FU = model.forward_use_variable.get_by_id(rxn.id)
        BU = model.backward_use_variable.get_by_id(rxn.id)
        BFUSE = model.add_variable(ForwardBackwardUseVariable, rxn)
        # BFUSE=model._var_dict['BFUSE_{}'.format(rxn.id)]

        F_ = model.reactions.get_by_id(rxn.id).forward_variable
        B_ = model.reactions.get_by_id(rxn.id).reverse_variable
        MFA = model.add_variable(MaximalFluxActiveVariable,rxn)
        #this variable is helper not required but used here to screen the results
        #could be also done by only screening BU
        MFAS = model.add_variable(MaximalFluxActiveVariableSum, rxn)

        #first add linearization =>MFAL=MFA*BU


        'if bfuse=1 BU=0 why we have this C --> ? from redGEMx '
        'not required'  # todo remove later
        """coupling constraint MFA and BU"""
        expression=MFA-BU

        cons4 = model.add_constraint(ConstraintMFA,
                                 rxn,
                                 expression,
                                 ub=0)

        expr_max_flux=B_-max_flux_non_c_source*MFA-(1-MFA)*big_M


        cons5 = model.add_constraint(ConstraintMaxFlux,
                                 rxn,
                                 expr_max_flux,
                                 ub=0)

        expr_new_var=MFAS-(BU+MFA)


        cons6 = model.add_constraint(ConstraintSum,
                                 rxn,
                                 expr_new_var,
                                 lb=0,ub=0)


        #bondinf constraint you have FU bcs it indicates that
        #with this main c ssource uptake there is also little secretion
        #indicating that it doesnt uptake all the little C sources at the same time and secretes
        #products of these little c sources if there is another c source that can do the job w less products
        #if interested the other can be run too
        #then the essential_met_pool will be higher resulting in more alternatives and more little c sources playing
        #with each other
        #1
        expr_bonding=BU-MFA+BFUSE+FU
        #2
        # expr_bonding=BU-MFA+BFUSE

        cons7 = model.add_constraint(ConstraintBond,
                                 rxn,
                                 expr_bonding,
                                 lb=0,ub=1.1)
    model.repair()
    #block essentials if you dont want to add essential substrate to a main C source list
    #as this can cause infeasibility of the other models and also lead to utilization of
    #this essential substrate only in big fluxes due to the catabolite repression constraint

    if block_essentials:
        for r in essential_list:
            rxn=model.reactions.get_by_id(r)
            BU_var=model._var_dict['BU_{}'.format(r)]
            MFA_var=model._var_dict['MFA_{}'.format(r)]
            expr_essential=(BU_var+MFA_var)-2.0*BU_var
            cons_ess = model.add_constraint(ConstraintEssential,
                                     rxn,
                                     expr_essential,
                                     lb=0)



        # mytfa_imm = imeAnalysis(mytfa, mytfa.exchanges)
    MFA_VARS = model._var_kinds[MaximalFluxActiveVariable.__name__]
    MFA_var_names = [x.name for x in MFA_VARS]

    # MFAL_VARS = model._var_kinds[MaximalFluxActiveVariableLinearize.__name__]
    # MFAL_var_names = [x.name for x in MFAL_VARS]

    MFAS_VARS = model._var_kinds[MaximalFluxActiveVariableSum.__name__]
    MFAS_var_names = [x.name for x in MFAS_VARS]

    BFUSE_VARS = model._var_kinds[ForwardBackwardUseVariable.__name__]
    BFUSE_var_names = [x.name for x in BFUSE_VARS]

    expression = sum(BFUSE_VARS) #before sum(MFA/MFAL)
    model.objective = expression
    model.objective_direction = 'max'
    sol = model.optimize()
    list_ = []
    df_all = pd.DataFrame()
    df_all['Alternative_1'] = np.zeros(len(model.variables))
    df_all.index = sol.raw.index

    df_sol = sol.raw
    # df_flux = sol.fluxes
    df_bfuse = df_sol.loc[BFUSE_var_names]
    idx = df_bfuse[df_bfuse < 0.5].index
    active_vars = [model._var_dict[k] for k in idx]
    # active_rxn_id = [i.reaction.id for i in active_vars]


    df_mfas=df_sol.loc[MFAS_var_names]
    idx_c_source = df_mfas[(df_mfas < 1.5) &(df_mfas > 0.5)].index
    active_vars_c_source=[model._var_dict[k] for k in idx_c_source]
    active_rxn_id_c_source = [i.reaction.id for i in active_vars_c_source]

    min_size_c=len(active_rxn_id_c_source)
    print('MIN SIZE high flux c_sources is {}'.format(min_size_c))

    'from here to be changed'
    if getAlternatives:
        alternative = 1
        while ((model.solver.status != INFEASIBLE) and (alternative <= max_alternative)):

            # todo look for primal from lUMPGEM
            df_sol = sol.raw
            df_flux = sol.fluxes
            df_bfuse = df_sol.loc[BFUSE_var_names]
            # todo check integer feasibility <1 but close to 1 or no !
            idx = df_bfuse[df_bfuse < 0.5].index
            active_vars = [model._var_dict[k] for k in idx]
            df = pd.DataFrame(columns=['metabolites', 'flux', 'FU', 'BU'])
            active_rxn_id = [i.reaction.id for i in active_vars]

            df_mfas = df_sol.loc[MFAS_var_names]
            # todo check integer feasibility <1 but close to 1 or no !
            #check the active ones with MFAS variable
            #here is where it is used for screening
            idx_c_source = df_mfas[(df_mfas < 1.5) & (df_mfas > 0.5)].index
            active_vars_c_source_mfas = [model._var_dict[k] for k in idx_c_source]
            active_rxn_id_c_source_mfas = [i.reaction.id for i in active_vars_c_source_mfas]
            active_vars_c_source = [k for k in model._var_kinds[MaximalFluxActiveVariable.__name__] if k.id in active_rxn_id_c_source_mfas]
            active_vars_c_source_bu = [k for k in model._var_kinds[BackwardUseVariable.__name__] if k.id in active_rxn_id_c_source_mfas]
            print('control is {} '.format(len(active_rxn_id_c_source_mfas)))

            if stop_cond:
                if (abs(len(active_rxn_id_c_source_mfas) - min_size_c) > criteria):
                    print(
                        'Stopping criteria is reached at alternative {} min+{} alternatives are exhausted for model'.format(
                            alternative, criteria, model.id))
                    break
            fu_vars = [model.forward_use_variable.get_by_id(k).name for k in active_rxn_id]
            bu_vars = [model.backward_use_variable.get_by_id(k).name for k in active_rxn_id]

            mfa_vars=[k.name for k in model._var_kinds[MaximalFluxActiveVariable.__name__] if k.id in active_rxn_id]
            # mfal_vars=[k.name for k in model._var_kinds[MaximalFluxActiveVariableLinearize.__name__] if k.id in active_rxn_id]
            mfas_vars=[k.name for k in model._var_kinds[MaximalFluxActiveVariableSum.__name__] if k.id in active_rxn_id]

            df['metabolites'] = active_rxn_id
            df['BFUSE'] = df_bfuse[df_bfuse < 0.5].values
            'error was here this var is only for forward take net flux'
            df['flux'] = df_flux[active_rxn_id].values
            df['FU'] = df_sol[fu_vars].values
            df['BU'] = df_sol[bu_vars].values
            df['MFA']=df_sol[mfa_vars].values
            # df['MFAL']=df_sol[mfal_vars].values
            df['MFAS']=df_sol[mfas_vars].values
            df['alternative'] = np.ones(len(active_rxn_id)) * alternative
            # add growth also as a column
            df['Growth'] = np.ones(len(active_rxn_id)) * sol.fluxes[biomass]

            df_all['Alternative_{}'.format(alternative)] = df_sol.values

            print('Generating alternative {} of size {} for model {}'.format(alternative, len(active_rxn_id),
                                                                             model.id))
            list_.append(df)
            #try to break
            if len(active_rxn_id_c_source)<1:
                break

            #add integer cuts based on the main c source
            expression2 = sum(active_vars_c_source_bu)
            cons1 = model.add_constraint(ModelConstraint,
                                          model,
                                          expression2,
                                          id_='max_flux_constraint_integer_cut_{}'.format(alternative),
                                          ub=0)
            model.repair()

            try:
                sol = model.optimize()
            except SolverError:
                print('Number of alternatives generated is ', alternative)
                break

            alternative += 1


    return list_, df_all

    return

