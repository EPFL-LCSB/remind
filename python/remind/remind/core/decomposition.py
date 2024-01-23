#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 09:41:45 2021

@author: ReMIND Team
"""
import pandas as pd
import numpy as np
from cobra import Reaction
from copy import deepcopy
import re

from pytfa.optim.constraints import ModelConstraint, ReactionConstraint
from pytfa.optim.variables import ForwardBackwardUseVariable
from pytfa.optim.utils import symbol_sum

from .exchanges import imeanalysis
from ..utils.postprocessing import get_indispensability_score
from .parsimonious import minimizealluptakes_by_weight


class ForwardUse(ReactionConstraint):
    """
    Class to represent a constraint the implements
    """
    prefix = 'FwU_'
    
class BackwardUse(ReactionConstraint):
    """
    Class to represent a constraint the implements
    """
    prefix = 'BwU_'
    

element_symbols = {
        'Fe': 'Iron',
        'Ca': 'Calcium',
        'Co': 'Cobalt',
        'Cu': 'Copper',
        'Cr': 'Chromium',
        'Mg': 'Magnesium',
        'Mn': 'Manganese',
        'Cl': 'Chlorine',
        'Zn': 'Zinc',
        'Na': 'Sodium',
        'Ni': 'Nickel',
        'Se': 'Selenium',
        'Br': 'Bromine',
        'Mo': 'Molybdenum',
        'Ce': 'Cerium',
        'Li': 'Lithium',
        'La': 'Lanthanum',
        'Hg': 'Mercury',
        'Ag': 'Silver',
        'Au': 'Gold',
        'As': 'Arsenic',
        'Cd': 'Cadmium',
        'C' : 'Carbon',
        'S' : 'Sulfur',
        'P' : 'Phosphorus',
        'O' : 'Oxygen',
        'N' : 'Nitrogen',
        'K' : 'Potassium',
        'F' : 'Flourine',
        'I' : 'Iodine',
        'U' : 'Uranium',
        'W' : 'Tungsten',
        'H' : 'Hydrogen',
        }

def find_elements(formula):
    
    """
    A fundction to find the elements in a metabolite
    
    """
    
    # element_symbols.update(additional_elements) # If new elements should be considered
    
    # first we should make sure the formula is clean and understandable
    try:
        if ';' in formula:
            formula.split(';')[0] # sometimes the formula is repeated; we just need one copy of it
        
        form = re.findall('[A-Z][^A-Z]*', formula) # split the string from capital letters
        element_syms = dict()
        for x in form:
            sym = re.findall(r'\D+', x)[0] # take the non-numbers
            try:
                num = re.findall(r'\d+', x)[0] # try to take numbers
            except IndexError:
                num = 1 # the coefficient is one and not written
            element_syms[sym] = num
        
        element_list = {}
        for sym, num in element_syms.items():
            the_element  = element_symbols[sym]
            element_list[the_element] = num
    except KeyError: # the formula is not know or it is a genereic formuila containing R or X
        generic_formula = {'Carbon':1, 'Nitrogen':1,
                           'Oxygen':1, 'Hydrogen':1, 'Phosphorus':1,
                           'Sulfur':1} # we can define a general formula containing main elements so that the elemental composition is reflected
        element_list = generic_formula
    
    return element_list


def find_elemental_exchange(model, element, reactions=None):
    
    """
    A fundtion to return all exchanges in the model that are associated with a specified element
    element is the name of the element in string, e.g, "Iron"
    """
    exchange_list = []
    if reactions is None:
        reactions = model.exchanges
    for rxn in reactions:
        for met in rxn.metabolites:
            element_list = find_elements(met.formula)
            if element in element_list.keys():
                exchange_list.append(rxn)
    
    return exchange_list

def decompose_biomass_block(model, growth_id, bbb, growth_limit):
    '''
    This adds an exchange for the specified bbb and also fixes its flux based
    on the stoichiometric coefficient and fixed growth.

    Parameters
    ----------
    model : pyTFA
    growth_id : str
    bbb : Metabolite
        The biomas building block.
    growth_limit : float
        The value of fixed growth.

    Returns
    -------
    ex_rxn : Reaction
        The exchange reaction that is added for the bbb.

    '''
    
    ex_rxn = Reaction('DCMP_{}'.format(bbb.id))
    ex_rxn.add_metabolites({bbb:-1})
    growth = model.reactions.get_by_id(growth_id)
    st_coef = growth.get_coefficient(bbb.id)
    bound = -st_coef * growth_limit
    ex_rxn.bounds = (bound,bound)
    model.add_reactions([ex_rxn])
    
    return ex_rxn


def find_essential_alternates(model_in, exchange_list, presolve_iter=100):
    '''
    Presolves for a specified number of MiExes and then, based on results, tries
    to find essential alternating metabolites:
        1) find the Exchanges that repeat in almost half of the cases
        2) find their couple that satisfies two condition: i) the intersection is
           empty and ii) the union is all.

    Parameters
    ----------
    model_in : Model
    presolve_iter : int, optional
        The number of solutions in presolve step. The default is 100.

    Returns
    -------
    df_list : list
        Similar to output of iMEs, instead of element, the hook is an alternating group.
    rxn_ids : list
        The IDs of the alternating reactions which should be removed.

    '''
    
    model = deepcopy(model_in) 
    # because we copy here the reactions should be re-specified
    new_exchage = [model.reactions.get_by_id(rxn.id) for rxn in exchange_list]
    
    # Presolve
    pre_list, _ = imeanalysis(model, new_exchage, element=None, 
                                  enforce_size=False,
                                  getAlternatives=True, 
                                  max_alternative=presolve_iter)
    frame = pd.concat(pre_list, ignore_index=True)
    alternative = frame.alternative.nunique()
    # Heuristically, the exchanges that appear in more than 1/5 of the total are suspicious
    threshold = (1/5) * alternative
    score = get_indispensability_score(frame)
    
    # find in which alternatoves the are repeated?
    rxn_dict = dict(tuple(frame.groupby('metabolites'))) # a dict where data grouped for each reaction (metabolite)
    susp_dict = dict()
    for rxn, rxn_data in rxn_dict.items():
        upt = score[rxn]['uptake']
        sec = score[rxn]['secretion']
        if upt > threshold and upt < alternative: # should not appear in all cases
            uptakes = \
                (rxn_data[np.isclose(rxn_data['BU'],1)])['alternative'].tolist()
            susp_dict[rxn] = uptakes
        elif sec > threshold and sec < alternative:
            secretions = \
                (rxn_data[np.isclose(rxn_data['FU'],1)])['alternative'].tolist()
            susp_dict[rxn] = secretions
    
    # check which couples have empty intersections (do not appear together)
    empty_intersect = dict()
    for key,val in susp_dict.items():
        empty_intersect[key] = [k for k,v in susp_dict.items() \
                                if k != key and len(set(val).intersection(v))==0]
    
    # we are looking for the essential ones, so the union should be all
    # since the intersect is empty, union is the sum
    remove = []
    for k,v in empty_intersect.items():
        union = len(susp_dict[k]) + sum([len(susp_dict[x]) for x in v])
        if union != alternative:
            remove.append(k)
    for k in remove:
        empty_intersect.pop(k)
        
    # finalizing the lists
    final_groups = []
    for k,v in empty_intersect.items():
        pot_group = (v + [k])
        pot_group.sort()
        if pot_group not in final_groups: # not to add redundants
            final_groups.append(pot_group)
    
    # create a DataFrame emulating the MiEx list
    df_list = []
    rxn_ids = []
    for ind, group in enumerate(final_groups):
        for alt, rxn_id in enumerate(group):
            df = pd.DataFrame(columns=['metabolites', 'flux', 'FU', 'BU'])
            df['metabolites'] = [rxn_id]
            df['flux'] = [np.nan] # not important
            df['FU'] = [1] if score[rxn_id]['secretion'] else [0]
            df['BU'] = [1] if score[rxn_id]['uptake'] else [0]
            df['BFUSE'] = [0]
            df['element'] = ['Group_{}'.format(ind)]
            df['alternative'] = [alt + 1]
            df_list.append(df)
            rxn_ids.append(rxn_id)

    df = pd.concat(df_list, ignore_index=True)

    return rxn_ids, df


def find_essential_alternates_from_dime_list(frame,thres=1/5):
    #it needs to be unique alternatives not max! if i am passing a part of aa dataframe
    alternative = frame.alternative.nunique()
    # Heuristically, the exchanges that appear in more than 1/5 of the total are suspicious
    threshold = (thres) * alternative
    score = get_indispensability_score(frame)

    # find in which alternatoves the are repeated?
    rxn_dict = dict(tuple(frame.groupby('metabolites')))  # a dict where data grouped for each reaction (metabolite)
    susp_dict = dict()
    for rxn, rxn_data in rxn_dict.items():
        upt = score[rxn]['uptake']
        sec = score[rxn]['secretion']
        if upt > threshold and upt < alternative:  # should not appear in all cases
            uptakes = \
                (rxn_data[np.isclose(rxn_data['BU'], 1)])['alternative'].tolist()
            susp_dict[rxn] = uptakes
        elif sec > threshold and sec < alternative:
            secretions = \
                (rxn_data[np.isclose(rxn_data['FU'], 1)])['alternative'].tolist()
            susp_dict[rxn] = secretions

    # check which couples have empty intersections (do not appear together)
    empty_intersect = dict()
    for key, val in susp_dict.items():
        empty_intersect[key] = [k for k, v in susp_dict.items() \
                                if k != key and len(set(val).intersection(v)) == 0]

    # we are looking for the essential ones, so the union should be all
    # since the intersect is empty, union is the sum
    "trial to comment"
    remove = []
    for k, v in empty_intersect.items():
        union = len(susp_dict[k]) + sum([len(susp_dict[x]) for x in v])
        if union != alternative:
            remove.append(k)
    for k in remove:
        empty_intersect.pop(k)
    "trial comment end"
    # finalizing the lists
    final_groups = []
    for k, v in empty_intersect.items():
        pot_group = (v + [k])
        pot_group.sort()
        if pot_group not in final_groups:  # not to add redundants
            final_groups.append(pot_group)

    # create a DataFrame emulating the MiEx list
    df_list = []
    rxn_ids = []
    for ind, group in enumerate(final_groups):
        for alt, rxn_id in enumerate(group):
            df = pd.DataFrame(columns=['metabolites', 'flux', 'FU', 'BU'])
            df['metabolites'] = [rxn_id]
            df['flux'] = [np.nan]  # not important
            df['FU'] = [1] if score[rxn_id]['secretion'] else [0]
            df['BU'] = [1] if score[rxn_id]['uptake'] else [0]
            df['BFUSE'] = [0]
            df['element'] = ['Group_{}'.format(ind)]
            df['alternative'] = [alt + 1]
            df_list.append(df)
            rxn_ids.append(rxn_id)
    df = pd.concat(df_list, ignore_index=True)

    return rxn_ids, df


def process_alternating_metabolites_data(frame,groups=['model','element']):
    fr = frame.groupby(groups).metabolites.unique().reset_index()
    fr['mets'] = [tuple(k) for k in fr.metabolites]
    return fr




def choose_removal_alternating_metabolites(fr,frame_all,model_col='id'):
    '''
    fr: dataframe containing groups to be removed generated from
    find_essential_alternates_from_dime_list function
    frame_all: dataframe containing all exchanges considering all models
    in the community
    model_col= column name which indicates the model
    selection is based on the frequency of appearance of metabolites among all models
    the ones that appears the least is considered to be removed heuristically!
    Though it still needs visual inspection of the results to avoid conflicting removals!
    '''
    to_be_removed=[]
    for k in fr.index:
        f=fr.iloc[k]
        choose_dict={}
        for idx in range(len(f.metabolites)):
            met=f.metabolites[idx]
            cnts=frame_all[frame_all.mets==met.replace('EX_','')][model_col].nunique()
            choose_dict[met]=cnts
        to_be_removed.append(min(choose_dict, key=choose_dict. get))
    fr['to_be_removed']=to_be_removed
    return fr

def rank_fdps(frame_fdp,frame_presolve):
    """frame_fdp:dataframe containing alternative fdps which is generated after solving presolve
 and then validating feasible fdps
 frame_presolve=this is the presolve dataframe from ime

returns:
f(dataframe) containing fdps unique_metabolites and its size ranked
    """
    f=frame_fdp.groupby('alternative').metabolites.unique().reset_index()
    fdp_alt_list={}
    fdp_list=[]
    for fdp in f.metabolites:

        union_alts={}
        combined_alternates=[]
        size_fdp=len(fdp)
        for met in fdp:
            alts=[k for k in frame_presolve[frame_presolve.metabolites==met].alternative.unique()]
            combined_alternates+=alts
            union_alts[met]=alts
        combined=pd.Series(combined_alternates).value_counts()
        combined=combined[combined==size_fdp].index.to_list()
        fdp_name = str(fdp)
        replacements={"'":"",
                      " ": "+",
                      "[": "",
                      "]": ""}
        #change fdp name accordingly
        for k,v in replacements.items():
            fdp_name=fdp_name.replace(k,v)
        #get fdp_name as the key and put the alternative numbers
        fdp_alt_list[fdp_name]=combined
        fdp_list.append(list(fdp))
    "here it is the choosing criteria for now it is based on # of unique metabolites"
    unique_mets={}
    for k,v in fdp_alt_list.items():
        unique_mets[k]=list(frame_presolve[frame_presolve.alternative.isin(v)].metabolites.unique())

    f['unique_mets']=[v for  k,v in unique_mets.items()]
    f['len_unique_mets']=[len(k) for k in f.unique_mets]
    # f['metabolites']=tuple(f.metabolites)
    f=f.sort_values('len_unique_mets', ascending=False)
    list_metabolites_union=[]
    for k in range(len(f)):
        list_metabolites_union+=f.iloc[k].unique_mets
    #check union to the selected
    list_metabolites_union=list(set(list_metabolites_union))
    selected_list=f.iloc[0].unique_mets
    inersection_list=set(selected_list).intersection(list_metabolites_union)
    not_common=[k for k in list_metabolites_union if k not in selected_list]
    not_common_of_interest=[k for k in not_common if k not in frame_fdp.metabolites.unique()]
    #there are uncommon mets if they are already fdp ones u dont care
    print('Selected fdp which has max unique metabolites has {} common mets and {} uncommon mets which is {}'.format(len(inersection_list),len(not_common_of_interest),not_common_of_interest))
    return f

def rank_fdps_from_find_alternates(frame_fdp,frame_presolve):
    """frame_fdp:dataframe taken directly from the function
    find_essential_alternates_from_dime_list
    frame_presolve=this is the presolve dataframe from ime

    returns:
    f(dataframe) containing fdps unique_metabolites and its size ranked
        """
    f=frame_fdp.groupby('alternative').metabolites.unique().reset_index()
    fdp_alt_list={}
    fdp_list=[]
    for fdp in f.metabolites:

        union_alts={}
        combined_alternates=[]
        size_fdp=len(fdp)
        for met in fdp:
            alts=[k for k in frame_presolve[frame_presolve.metabolites==met].alternative.unique()]
            combined_alternates+=alts
            union_alts[met]=alts
        combined=pd.Series(combined_alternates).value_counts()
        combined=combined[combined==size_fdp].index.to_list()
        fdp_name = str(fdp)
        replacements={"'":"",
                      " ": "+",
                      "[": "",
                      "]": ""}
        #change fdp name accordingly
        for k,v in replacements.items():
            fdp_name=fdp_name.replace(k,v)
        #get fdp_name as the key and put the alternative numbers
        fdp_alt_list[fdp_name]=combined
        fdp_list.append(list(fdp))
    "here it is the choosing criteria for now it is based on # of unique metabolites"
    unique_mets={}
    for k,v in fdp_alt_list.items():
        unique_mets[k]=list(frame_presolve[frame_presolve.alternative.isin(v)].metabolites.unique())

    f['unique_mets']=[v for  k,v in unique_mets.items()]
    f['len_unique_mets']=[len(k) for k in f.unique_mets]
    # f['metabolites']=tuple(f.metabolites)
    f=f.sort_values('len_unique_mets', ascending=False)
    list_metabolites_union=[]
    for k in range(len(f)):
        list_metabolites_union+=f.iloc[k].unique_mets
    #check union to the selected
    list_metabolites_union=list(set(list_metabolites_union))
    selected_list=f.iloc[0].unique_mets
    inersection_list=set(selected_list).intersection(list_metabolites_union)
    not_common=[k for k in list_metabolites_union if k not in selected_list]
    not_common_of_interest=[k for k in not_common if k not in frame_fdp.metabolites.unique()]
    #there are uncommon mets if they are already fdp ones u dont care
    print('Selected fdp which has max unique metabolites has {} common mets and {} uncommon mets which is {}'.format(len(inersection_list),len(not_common_of_interest),not_common_of_interest))
    return f

def check_feasibility(org_model, growth_id, growth_limit, yield_cut, subs, prods,
                      water_exch = 'EX_h2o_e', proton_exch = 'EX_h_e'):
    '''
    This function checks if the assembled DiME is feasible.

    Parameters
    ----------
    org_model : pytfa.Model
        DESCRIPTION.
    growth_id : str
        DESCRIPTION
    growth_limit : float
        DESCRIPTION
    yield_cut : float
        DESCRIPTION
    subs : list
        substrates for this solution.
    prod : list
        products for this solution.

    Returns
    -------
    state : str
        "feasible" or "infeasible".

    '''
    
    with org_model as temp_model:
        '4)add constrain for parsimonious'
        'minimize uptake of carbon sources or all uptakes user defined'
        # tolerance defined for each given yield cut 
        temp_model.reactions.get_by_id(growth_id).lower_bound = growth_limit
        temp_model.reactions.get_by_id(growth_id).upper_bound = growth_limit
        
        # find the minimal uptake by taking parsimonious uptake
        # sometimes it is better to couple and sometimes it is better to decouple
        '3)apply parsimonious and get solution and expression'
        'minimize uptake of carbon sources or all uptakes'
        expr_pars, sol_pars=minimizealluptakes_by_weight(temp_model, molecular_weight='formula_weight')
        min_uptake = sol_pars.objective_value
        
        tolerance = 1 / yield_cut
        temp_model.add_constraint(ModelConstraint,
                                          temp_model,
                                          expr_pars,
                                          id_='pars_constraint',
                                          # fix uptake
                                          lb=min_uptake * tolerance,
                                          ub=min_uptake * tolerance)

        temp_model.repair()


        for rxn in temp_model.exchanges:
            rxn.bounds = (0, 0)

        for id_ in subs:
            temp_model.reactions.get_by_id(id_).bounds = (-1000, 0)

        for id_ in prods:
            temp_model.reactions.get_by_id(id_).bounds = (0, 1000)
            
        # water and H+ exchenge is not necessarily accounted for, so open them anyway
        temp_model.reactions.get_by_id(water_exch).bounds = (-1000, 1000)
        temp_model.reactions.get_by_id(proton_exch).bounds = (-1000, 1000)

        growth = temp_model.slim_optimize()
        
        state = 'infeasible' if np.isnan(growth) else 'feasible'
    
    
    
    return state


def recover_infeas_sol(org_model, query_rxns, growth_id, growth_limit, yield_cut, subs, prods,
                      water_exch = 'EX_h2o_e', proton_exch = 'EX_h_e'):
    '''
    This function checks if the assembled DiME is feasible.

    Parameters
    ----------
    org_model : pytfa.Model
        DESCRIPTION.
    query_rxns : list
        DESCRIPTION.
    growth_id : str
        DESCRIPTION
    growth_limit : float
        DESCRIPTION
    yield_cut : float
        DESCRIPTION
    subs : list
        substrates for this solution.
    prod : list
        products for this solution.

    Returns
    -------
    state : str
        "feasible" or "infeasible".

    '''
    
    with org_model as new_model:
        '4)add constrain for parsimonious'
        'minimize uptake of carbon sources or all uptakes user defined'
        # tolerance defined for each given yield cut 
        new_model.reactions.get_by_id(growth_id).lower_bound = growth_limit
        new_model.reactions.get_by_id(growth_id).upper_bound = growth_limit
        
        # find the minimal uptake by taking parsimonious uptake
        # sometimes it is better to couple and sometimes it is better to decouple
        '3)apply parsimonious and get solution and expression'
        'minimize uptake of carbon sources or all uptakes'
        expr_pars, sol_pars=minimizealluptakes_by_weight(new_model, molecular_weight='formula_weight')
        min_uptake = sol_pars.objective_value
        
        tolerance = 1 / yield_cut
        new_model.add_constraint(ModelConstraint,
                                          new_model,
                                          expr_pars,
                                          id_='pars_constraint',
                                          # fix uptake
                                          lb=min_uptake * tolerance,
                                          ub=min_uptake * tolerance)

        new_model.repair()


        # First close all exchanges
        for rxn in new_model.exchanges:
            rxn.bounds = (0, 0)

        for id_ in subs:
            new_model.reactions.get_by_id(id_).bounds = (-1000, 0)

        for id_ in prods:
            new_model.reactions.get_by_id(id_).bounds = (0, 1000)
            
        bfuse_vars = []
        bfuse_cons = []
        bigM = 1000
        for id_ in query_rxns:
            rxn = new_model.reactions.get_by_id(id_)
            rxn.bounds = (-1000,1000)
            bfuse = new_model.add_variable(ForwardBackwardUseVariable,rxn)
            bfuse_vars += [bfuse]
            cons1 = new_model.add_constraint(ForwardUse,
                                     rxn,
                                     rxn.flux_expression - bigM*bfuse,
                                     ub = 0)
            cons2 = new_model.add_constraint(BackwardUse,
                                     rxn,
                                     rxn.flux_expression + bigM*bfuse,
                                     lb = 0)
            bfuse_cons += [cons1, cons2]
            
        # water and H+ exchenge is not necessarily accounted for, so open them anyway
        new_model.reactions.get_by_id(water_exch).bounds = (-1000, 1000)
        new_model.reactions.get_by_id(proton_exch).bounds = (-1000, 1000)

        new_model.objective = symbol_sum(bfuse_vars)
        new_model.objective_direction = 'min'
        new_model.slim_optimize()
        
        active_vars = [var.name for var in bfuse_vars if np.isclose(var.variable.primal,1)]
        
        # context specific coding does not remove added variables and constraints properly
        for var in bfuse_vars:
            new_model.remove_variable(var)
        for cons in bfuse_cons:
            new_model.remove_constraint(cons)
    
    
    return active_vars

    