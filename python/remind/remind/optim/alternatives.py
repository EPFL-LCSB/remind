#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:56:04 2021

@author: ReMIND Team
"""
import numpy as np
import pandas as pd

from optlang.interface import INFEASIBLE
from pytfa.optim.utils import symbol_sum

from .constraints import ModelConstraint
from .variables import SecretionActivation, UptakeActivation, \
    PositiveInteraction, NegativeInteraction,YieldUse,DIMEVariable
from ..utils.postprocessing import get_active_rxn_yield
from ..core.decomposition import check_feasibility, recover_infeas_sol
from ..optim.variables import AbioticResource

def extract_sol(rec_ids, species_id):
    
    raw = {'RC_' + species_id + '_' + id_.replace('BFUSE_',''):{'d':1} \
           for id_ in rec_ids}
    raw = pd.DataFrame.from_dict(raw, orient='index')
    sol = raw['d']
    
    return sol


def find_alternative_solutions(model, growth_ids=[], growth_limits=[],
                               alternation = 'reaction', max_alts=100,
                               check_feasibility = False,stop_cond=True,criteria=2,\
                               tolerance=1e-7):
    '''
    Finds alternative solution for a model, regardless of the objective function and
    by adding integer cuts. The integer cuts can be defined for the interactions or
    reactions. Also, checks for the feasibility of the solutions. For feasibility
    checking, we should anyway add integer cuts based on reactions. If the solution
    is feasible, the integer cuts based on interactions can be applied
    alternation reaction: considering all DiMES for each species
    alternation interaction: considering the ids of rxns that give the interaction

    Parameters
    ----------
    model : InteractionModel
        DESCRIPTION.
    growth_ids : dict
        DESCRIPTION
    growth_limits : dict
        DESCRIPTION
    alternation : 'reaction' or 'interaction', optional
        DESCRIPTION. The default is 'reaction'.
    max_alts : int, optional
        DESCRIPTION. The default is 100.
    check_feasibility : Boolean
        To check the feasibility of DiMEs? The default is False

    Returns
    -------
    df : pd.DataFrame
        DESCRIPTION.

    '''
    #stores all solution.raw
    sol_list = []
    #stores only active variables like in dimes format
    list_active=[]
    
    # first solution
    alternative = 1 # hashing and counting feasible sols
        
    # these are used to remove the infeasible solutions
    basic_int_cut_vars = [var for var in model.get_variables_of_type(SecretionActivation)] + \
        [var for var in model.get_variables_of_type(UptakeActivation)]
        
    
    # these are used in case that the solution is feasible     
    if alternation == 'reaction':
        int_cut_vars = [var for var in model.get_variables_of_type(SecretionActivation)] + \
        [var for var in model.get_variables_of_type(UptakeActivation)]

    elif alternation == 'reaction_general':
        int_cut_vars = [var for var in model.get_variables_of_type(SecretionActivation)] + \
        [var for var in model.get_variables_of_type(UptakeActivation)]
        
    elif alternation == 'interaction':
        int_cut_vars = [var for var in model.get_variables_of_type(PositiveInteraction)] + \
            [var for var in model.get_variables_of_type(NegativeInteraction)]
            
    elif alternation == 'interaction_directional':
        int_cut_vars = [var for var in model.get_variables_of_type(PositiveInteraction)] + \
            [var for var in model.get_variables_of_type(NegativeInteraction)]

    elif alternation=='interaction_positive_directional_general':
        int_cut_vars = [var for var in model.get_variables_of_type(PositiveInteraction)]
        
    elif alternation=='interaction_negative_directional_general':
        int_cut_vars = [var for var in model.get_variables_of_type(NegativeInteraction)]

    elif alternation == 'interaction_positive':
        int_cut_vars = [var for var in model.get_variables_of_type(PositiveInteraction)]

    elif alternation == 'interaction_negative':
        int_cut_vars = [var for var in model.get_variables_of_type(NegativeInteraction)]

    elif alternation == 'abiotic_sources':
        int_cut_vars = [var for var in model.get_variables_of_type(AbioticResource)]


#new integer cuts withyield_cut also added
#todo can be generalized

    elif alternation == 'interaction_yield':
        int_cut_vars = [var for var in model.get_variables_of_type(PositiveInteraction)] + \
                       [var for var in model.get_variables_of_type(NegativeInteraction)] + \
                       [var for var in model.get_variables_of_type(YieldUse)]
            #yieldvars
            
    elif alternation == 'interaction_directional_yield':
        int_cut_vars = [var for var in model.get_variables_of_type(PositiveInteraction)] + \
                       [var for var in model.get_variables_of_type(NegativeInteraction)] + \
                       [var for var in model.get_variables_of_type(YieldUse)]
            #yieldvars

    elif alternation == 'interaction_positive_yield':
        int_cut_vars = [var for var in model.get_variables_of_type(PositiveInteraction)] + \
                       [var for var in model.get_variables_of_type(YieldUse)]

    elif alternation == 'interaction_positive_directional_general_yield':

        int_cut_vars = [var for var in model.get_variables_of_type(PositiveInteraction)] + \
                       [var for var in model.get_variables_of_type(YieldUse)]

    elif alternation == 'interaction_negative_directional_general_yield':

        int_cut_vars = [var for var in model.get_variables_of_type(NegativeInteraction)] + \
                       [var for var in model.get_variables_of_type(YieldUse)]

    elif alternation == 'interaction_negative_yield':
        int_cut_vars = [var for var in model.get_variables_of_type(NegativeInteraction)] + \
                       [var for var in model.get_variables_of_type(YieldUse)]

    elif alternation == 'abiotic_sources_yield':
        int_cut_vars = [var for var in model.get_variables_of_type(AbioticResource)] + \
                       [var for var in model.get_variables_of_type(YieldUse)]

    elif alternation == 'abiotic_sources_uptake':
        int_cut_vars = [var for var in model.get_variables_of_type(AbioticResource)]

    elif alternation=="dime_variable":
        int_cut_vars = [var for var in model.get_variables_of_type(DIMEVariable)]


    else:
        raise ValueError('The alternations can be based on either interaction or reaction.')
    


    " here alreadyy solve for once to store the ooptimal value if infeasible"
    try:
        sol_opt= model.optimize()
        optimal_obj=sol_opt.objective_value

    except Exception:
        print('model is not feasible')

    # generating alternative solution
    while((model.solver.status != INFEASIBLE) and (alternative <= max_alts)):
        print('Generating alternative', alternative)
        try:
            sol = model.optimize()
            solution = sol.raw
            obj=sol.objective_value

            #aslis data
            df_active=pd.DataFrame(solution[solution.values>=0.5].index.to_list(),columns=['variables'])
            df_active['alternative']=int(alternative)*np.ones(df_active.shape[0])
            print('Objective is {}'.format(obj))
            df_active['objective']=obj*np.ones(df_active.shape[0])


            solution_all=solution.to_frame('Alternative_{}'.format(alternative))
            solution_all.loc['objective'] = obj


            print('Generating alternative {} and objective is {}'.format(alternative,obj))


        except Exception:
            print('Number of alternatives generated is ', alternative)
            break

        if stop_cond:
            #if abs is used direction wont be necessary
            #direction=model.objective.direction
                if (abs(obj-optimal_obj)>(criteria+tolerance)):
                    print('Stopping criteria is reached at alternative {} min+{} alternatives are exhausted for model'.format(alternative, criteria))
                    break

        # check the feasibility for all the species in the community
        if check_feasibility:
            for sp in model.species:
                yield_, substrates, products = get_active_rxn_yield(model, sp.id)
                growth_id = growth_ids[sp.id]
                growth_limit = growth_limits[sp.id]
                state = check_feasibility(sp, growth_id, growth_limit, 
                                          yield_, substrates, products)
                if state == 'infeasible':
                    
                    
                    # TODO: if a solution is once recovered, we do not have to do it again
                    # Add integer cut so this solution is not found again
                    print('An infeasible solution was found for {} and is being recycled.'.format(sp.id))
                    
                    query_rxns = find_query_rxns(model, sp, substrates, products)
                    
                    rec_ids = recover_infeas_sol(sp, query_rxns, growth_id, growth_limit, 
                                          yield_, substrates, products)
                    new_sol = extract_sol(rec_ids, sp.id) # preparing the IDs for recovery to be saved
                    solution = pd.concat([solution,new_sol])

            
        # adding the integer cuts anyway: 
        alternative += 1
        sol_list.append(solution_all)
        list_active.append(df_active)
        # find the inactive/active variables
        #for all integer cuts active is noted with =1 so it works
        #todo keep in mind if we look for smth active=0
        act_vars = [var for var in int_cut_vars if np.isclose(var.variable.primal, 1)]
        inact_vars = [var for var in int_cut_vars if np.isclose(var.variable.primal, 0)]
        # For the interactions add basic integer cuts, but for reactions add non-basic (ours)
        #reasoning dime variable is anyways one per species and it cannot be more for now! so our
        #integer cuts should work
        if alternation in['reaction','abiotic_sources',\
                        'reaction_yield','abiotic_sources_yield','dime_variable']:
            expr = symbol_sum([1-var for var in act_vars])
        elif alternation in ['interaction','interaction_positive','interaction_negative'\
                ,'interaction_yield','interaction_positive_yield','interaction_negative_yield','reaction_general']:
            expr = symbol_sum([var for var in inact_vars]) + \
                symbol_sum([1-var for var in act_vars])


        #this general interaction is a more recent one
        elif alternation in ['interaction_positive_directional_general','interaction_positive_directional_general_yield']:
            #increase active vars
            pos_hook=[var for var in act_vars]
            pos_int_hook = [k.hook for k in pos_hook if k.name.startswith('PI_')]
            upt_act = [k for k in model.get_variables_of_type(UptakeActivation) if k.hook.metabolite in pos_int_hook]
            sec_act = [k for k in model.get_variables_of_type(SecretionActivation) if k.hook.metabolite in pos_int_hook]
            additional=upt_act+sec_act
            additional_act_vars = [var for var in additional if np.isclose(var.variable.primal, 1)]
            act_vars+=additional_act_vars
            additional_inac_vars=[var for var in additional if np.isclose(var.variable.primal, 0)]
            inact_vars+=additional_inac_vars
            expr = symbol_sum([var for var in inact_vars]) + \
                   symbol_sum([1 - var for var in act_vars])
                   
        elif alternation in ['interaction_negative_directional_general','interaction_negative_directional_general_yield']:
            #increase active vars
            neg_hook=[var for var in act_vars]
            neg_int_hook = [k.hook for k in neg_hook if k.name.startswith('NI_')]
            upt_act = [k for k in model.get_variables_of_type(UptakeActivation) if k.hook.metabolite in neg_int_hook]
            sec_act = [k for k in model.get_variables_of_type(SecretionActivation) if k.hook.metabolite in neg_int_hook]
            additional=upt_act+sec_act
            additional_act_vars = [var for var in additional if np.isclose(var.variable.primal, 1)]
            act_vars+=additional_act_vars
            additional_inac_vars=[var for var in additional if np.isclose(var.variable.primal, 0)]
            inact_vars+=additional_inac_vars
            expr = symbol_sum([var for var in inact_vars]) + \
                   symbol_sum([1 - var for var in act_vars])
                   
        elif alternation in ['interaction_directional','interaction_directional_yield']:
            #increase active vars
            active_hooks=[var for var in act_vars]
            int_hook = [k.hook for k in active_hooks if k.name.startswith('NI_') or k.name.startswith('PI_')]
            upt_act = [k for k in model.get_variables_of_type(UptakeActivation) if k.hook.metabolite in int_hook]
            sec_act = [k for k in model.get_variables_of_type(SecretionActivation) if k.hook.metabolite in int_hook]
            additional=upt_act+sec_act
            additional_act_vars = [var for var in additional if np.isclose(var.variable.primal, 1)]
            act_vars+=additional_act_vars
            additional_inac_vars=[var for var in additional if np.isclose(var.variable.primal, 0)]
            inact_vars+=additional_inac_vars
            expr = symbol_sum([var for var in inact_vars]) + \
                   symbol_sum([1 - var for var in act_vars])

    #todo check for 7 member species and run
        elif alternation in ['abiotic_sources_uptake']:
            active_hooks=[var for var in act_vars]
            abiot_hook = [k.hook for k in active_hooks if k.name.startswith('AV_')]
            upt_act = [k for k in model.get_variables_of_type(UptakeActivation) if k.hook.metabolite in abiot_hook]
            additional=upt_act
            additional_act_vars = [var for var in additional if np.isclose(var.variable.primal, 1)]
            act_vars+=additional_act_vars
            additional_inac_vars=[var for var in additional if np.isclose(var.variable.primal, 0)]
            inact_vars+=additional_inac_vars
            #if now another species uptakes you need to capture
            expr = symbol_sum([var for var in inact_vars]) + \
                   symbol_sum([1 - var for var in act_vars])




        #as we stick to min abiotic sources in min +1 a superset should not be possible
        #if we are interested in other
        # if alternation == 'abiotic_sources':
        #     expr = symbol_sum([1-var for var in act_vars])
        
        model.add_constraint(ModelConstraint,
                                 model,
                                 expr,
                                 id_='feas_integer_cut_{}'.format(alternative),
                                 lb = 1)
        # else:
        #     ind += 1
        #     for sp_id, expr in infeasible_integer_cuts.items():
        #         model.add_constraint(ModelConstraint,
        #                                  model,
        #                                  expr,
        #                                  id_='infeas_integer_cut_{}_{}'.format(sp_id,ind),
        #                                  lb = 1)
        

            
    df_all = pd.concat(sol_list, axis=1)
    df_active_vars = pd.concat(list_active, ignore_index=True)
    return df_all,df_active_vars







def find_query_rxns(interaction_model, infeasible_model,
                                   substrates, products):
    '''
    This function finds the f exchange reactions that can be 
    opened to recover an infeasible solution (the solution is generated by assembling
    DiMEs). This function ensures that opening the new exchange does not affect 
    the interactions. Hence, the current products, substrates, and the reactions
    that contibute to the interactions are removed.

    Parameters
    ----------
    interaction_model : InteractionModel
        A solved interaction model.
    infeasible_model : CommunitySpecies
        An infeasible model for an organism.
    substrates : list
        subtrates for this solution.
    products : list
        products for this solution.

    Returns
    -------
    query_rxns : list
        Reaction IDs of exchanges that can be opened to recover the solution.

    '''
    
    active_interactions = [var for var in infeasible_model.get_variables_of_type(PositiveInteraction) \
                           if np.isclose(var.variable.primal,1)] + \
            [var for var in infeasible_model.get_variables_of_type(NegativeInteraction) \
             if np.isclose(var.variable.primal,1)]
                
    # find the exchange IDs that participate in the interactions
    exch_int_act = [var.name.replace('PI__','EX_').replace('NI__','EX_') \
                    for var in active_interactions]
    
    exclude_rxns = substrates + products + exch_int_act # these reactions cannot be open to recover the solution
    
    query_rxns = [rxn.id for rxn in infeasible_model.exchanges if rxn.id not in exclude_rxns]
    
    
    return query_rxns
