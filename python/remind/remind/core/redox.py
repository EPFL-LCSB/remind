#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 19:05:45 2021

@author: ReMIND Team
"""
import re
from itertools import combinations

ATP = 'ATP'
ADP = 'ADP'
COA = 'CoA'
CO2 = 'CO2'
NH4 = 'NH4'
NO2 = 'NO2'
NO3 = 'NO3'
SO4 = 'SO4'
SO3 = 'SO3'
CO = 'CO'
CH4 = 'CH4'
PO4 = 'PO4'
PO3 = 'PO3'
PPI = 'PPi'
H2O2 = 'H2O2'
H2O = 'H2O'
H2 = 'H2'
O2 = 'O2'
N2 = 'N2'
H = 'H'


def dict_compare(d1, d2):
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    added = d1_keys - d2_keys
    removed = d2_keys - d1_keys
    modified = {o : (d1[o], d2[o]) for o in shared_keys if d1[o] != d2[o]}
    same = set(o for o in shared_keys if d1[o] == d2[o])
    return added, removed, modified, same

def check_redox_def(potent_pair):
    '''
    It tries to find a redox pair based on the definition and by comparing the formulae
    By definition, three changes are observed in redox pairs:
        1) the formula is exactly the same (then the charge is different)
        2) The difference is only in the number of hydrogens
        3) The difference is only in the number of oxygens
    Since I am not able to exclude transport between compartments, I did not
    check the first definition.
    
    Parameters
    ----------
    potent_pair : TYPE
        DESCRIPTION.

    Returns
    -------
    found_pair

    '''
    found_pair = []
    for x in combinations(potent_pair, 2): # taking all pairs
        met_1, met_2 = x[0], x[1]
        formula_1 = find_element_coeff(met_1.formula.split(';')[0]) # CarveMe models!
        formula_2 = find_element_coeff(met_2.formula.split(';')[0])
        added, removed, modified, _ = dict_compare(formula_1, formula_2)
        changed = list(added) + list(removed) + list(modified)
        if len(changed) == 1 and ('H' in changed or 'O' in changed):
            found_pair += [met_1, met_2]
    
    return found_pair

def electron_source_sink(model, redox_pairs, extracellular = '_e', 
                         compartments = ['_c','_p']):
    '''
    This function finds which redox pairs are present in extracellular.
    These are the potential external electron acceptors and donors.
    This function is written for BiGG IDs. For other annotations, may need to be modified.

    Parameters
    ----------
    model : TYPE
        DESCRIPTION.
    reox_pair : list of tuples
        This is the output of find_redox_pairs(model, rxns, oxidized, reduced)

    Returns
    -------
    sources_sinks : TYPE
        DESCRIPTION.

    '''
    sources_sinks = []
    
    for pair in redox_pairs:
        for met_id in pair:
            for comps in compartments:
                if extracellular in met_id:
                    excell_id = met_id
                elif comps in met_id:
                    excell_id = met_id.replace(comps, extracellular)
            try:
                met = model.metabolites.get_by_id(excell_id)
                if met.id not in sources_sinks:
                    sources_sinks.append(met.id)
            except KeyError: # it is only internal
                pass
    
    return sources_sinks


def find_element_coeff(met_formula):
    '''
    

    Parameters
    ----------
    met_formula : str
        metabolic chemical formula.

    Returns
    -------
    element_dict : dict
        keys are the name of element and values are their number

    '''
    if ';' in met_formula:
        met_formula.split(';')[0] # sometimes the formula is repeated; we just need one copy of it
    
    # the coefficient 1 is not written!
    el_spl = re.findall('[A-Z][^A-Z]*',met_formula) # split based on capital letters
    for ind, sym in enumerate(el_spl):
        coeff = [s for s in sym if s.isdigit()]
        if len(coeff) == 0: # no coefficient for this element
            el_spl[ind] = sym + '1'
    
    formula = ''.join(el_spl)
    
    splitted = re.split('(\d+)',formula)
    # remove nonfunctional characters
    splitted = [s for s in splitted if s!='']
                
    
    # the odd indices are symbols and the even indices are numbers
    length = int(len(splitted)/2)
    assert int(len(splitted)/2) == len(splitted)/2
        
    element_dict = {splitted[2*i]:int(splitted[2*i+1]) for i in range(0, length)}
                                  
    
    # remove Hydrogen if exists
    # element_dict.pop('H',1)
    
    return element_dict

def find_electron_cycle(model, growth_id, nad='nad_c', nadh='nadh_c'):
    '''
    This function starts from NAD/NADH to find redox reactions. By comparing metabolite
    formulae, it tries to find redox pairs in the model. 

    Parameters
    ----------
    model : pytfa.model
    growth_id : str
    nad : str, optional
        The metabolite id for NAD. The default is 'nad_c'.
    nadh : str, optional
        The metabolite id for NADH. The default is 'nadh_c'.

    Returns
    -------
    redox_mets : list of tuples
        A dictionary of redox pairs.

    '''
    redox_mets = [] # starts with NAD/NADH
    
    check_list = [(nad,nadh)] # the set of active redox pairs
    checked_rxns = [] # If we already checked a reaction we don't have to see it the second time
    
    not_found_rxns = []
    while len(check_list) !=0: # as long as there is something to check
        query = check_list[0]
        oxidized = model.metabolites.get_by_id(query[0])
        reduced = model.metabolites.get_by_id(query[1])
        if query not in redox_mets:
            redox_mets += [query]
        else:
            check_list.remove(query)
            continue
        redox_rxns = [rxn for rxn in oxidized.reactions if rxn in \
                      reduced.reactions and rxn.id != growth_id]
        # make sure that the redox pair appears on different sides
        redox_rxns = [r for r in redox_rxns if \
              r.get_coefficient(oxidized.id)*r.get_coefficient(reduced.id) < 0]
        # remove the reactions that have been already checked
        redox_rxns = [r for r in redox_rxns if r not in checked_rxns]
        
        new_pairs, not_found = find_redox_pairs(model, redox_rxns, oxidized, reduced)
        check_list += new_pairs
        not_found_rxns += not_found
        
        # update the checked reactions
        checked_rxns += redox_rxns
    
    return redox_mets, not_found_rxns


def find_redox_pairs(model, rxns, oxidized, reduced):
    '''
    This function finds new redoc pairs based on redox reactions
    and already found pair (oxidized, reduced)

    Parameters
    ----------
    rxns : list of reaction objects
        Redox reactions that include oxidized and reduced.
    oxidized : Metabolite
    reduced : Metabolite

    Returns
    -------
    new_pairs : list of tuples

    '''
    not_found = []
    
    new_pairs = []
    current_pair = [oxidized, reduced]
    for rxn in rxns:
        potent_pair = [met for met in rxn.metabolites if met not in current_pair]
        max_stoc_coeff = max([abs(rxn.get_coefficient(m.id)) for m in rxn.metabolites])
        if len(rxn.metabolites) == 3 and max_stoc_coeff > 1:
            # it is possible that one compound serves as electron donor and acceptor
            # but it's only possible if the stoichiometric coefficient is more than 1
            # e.g.: 2 H2O2 --> 2H2O + O2
            assert len(potent_pair) == 1 
            curr_pot = potent_pair[0]
            # take the metabolite that is on the other side of the potent_pair metabolite
            new_pot = [met for met in rxn.metabolites if \
                       rxn.get_coefficient(met.id)*rxn.get_coefficient(curr_pot.id)<0]
            potent_pair += new_pot
            new_ox = [met.id for met in potent_pair if \
                      rxn.get_coefficient(met.id)*rxn.get_coefficient(oxidized.id)<0][0]
            new_red = [met.id for met in potent_pair if \
                       rxn.get_coefficient(met.id)*rxn.get_coefficient(reduced.id)<0][0]
            new_pairs += [(new_ox,new_red)]
        # it is possible that H2 is oxidized to H+ (this is captured in the first)
        elif len(rxn.metabolites) == 4: # the easiest case
            # the other two metabolites are redox pair
            # the new oxidized form is on different side of the current oxidized
            try:
                new_ox = [met.id for met in rxn.metabolites if met not in current_pair \
                      and rxn.get_coefficient(met.id)*rxn.get_coefficient(oxidized.id)<0][0]
                new_red = [met.id for met in rxn.metabolites if met not in current_pair \
                      and rxn.get_coefficient(met.id)*rxn.get_coefficient(reduced.id)<0][0]
                new_pairs += [(new_ox,new_red)]
            except IndexError: # it is possible we have non-oxidative (de)carboxylation
                continue
       
        else: # more complicated cases
            spec_new_pairs = find_special_redox_pairs(rxn,potent_pair)
            if len(spec_new_pairs) != 2:
                not_found += [rxn]
                continue
            new_ox = [met.id for met in spec_new_pairs if \
                      rxn.get_coefficient(met.id)*rxn.get_coefficient(oxidized.id)<0][0]
            new_red = [met.id for met in spec_new_pairs if \
                       rxn.get_coefficient(met.id)*rxn.get_coefficient(reduced.id)<0][0]
            new_pairs += [(new_ox,new_red)]
    
    return new_pairs, not_found


def find_special_redox_pairs(rxn, potent_pair):
    
    assert len(potent_pair) > 2 # make sure it is really a complicated redox reaction
    
   
    # first try to find the well-known metabolites, as the rules are defined based on them
    known_mets = {}
    metalion = [] # a list to keep potential metal ions
    
    for met in potent_pair:
        met_formula = find_element_coeff(met.formula)
        # ATP
        try:
            if met_formula['C'] == 10 and met_formula['H'] > 10 and met_formula['N'] == 5 \
            and met_formula['O'] == 13 and met_formula['P'] == 3:
                known_mets[ATP] = met
        except KeyError:
            pass
        # ADP
        try:
            if met_formula['C'] == 10 and met_formula['H'] > 10 and met_formula['N'] == 5 \
            and met_formula['O'] == 10 and met_formula['P'] == 2:
                known_mets[ADP] = met
        except KeyError:
            pass
        # CoA
        try:
            if met_formula['C'] == 21 and met_formula['H'] == 32 and met_formula['N'] == 7 \
            and met_formula['O'] == 16 and met_formula['P'] == 3 and met_formula['S'] == 1:
                known_mets[COA] = met
        except KeyError:
            pass
        # Ammonium
        try:
            if 'C' not in met_formula and met_formula['H'] == 4 and met_formula['N'] == 1:
                known_mets[NH4] = met
        except KeyError:
            pass
        # Nitrite
        try:
            if 'C' not in met_formula and met_formula['O'] == 2 and met_formula['N'] == 1:
                known_mets[NO2] = met
        except KeyError:
            pass
        # Nitrate
        try:
            if 'C' not in met_formula and met_formula['O'] == 3 and met_formula['N'] == 1:
                known_mets[NO3] = met
        except KeyError:
            pass
        # CO2
        try:
            if len(met_formula)==2 and met_formula['O'] == 2 and met_formula['C'] == 1:
                known_mets[CO2] = met
        except KeyError:
            pass
        # CO
        try:
            if len(met_formula)==2 and met_formula['O'] == 1 and met_formula['C'] == 1:
                known_mets[CO] = met
        except KeyError:
            pass
        # Methane
        try:
            if 'O' not in met_formula and met_formula['H'] == 4 and met_formula['C'] == 1:
                known_mets[CH4] = met
        except KeyError:
            pass
        # Sulfite
        try:
            if 'C' not in met_formula and met_formula['O'] == 3 and met_formula['S'] == 1:
                known_mets[SO3] = met
        except KeyError:
            pass
        # Sulfate
        try:
            if 'C' not in met_formula and met_formula['O'] == 4 and met_formula['S'] == 1:
                known_mets[SO4] = met
        except KeyError:
            pass
        # Phosphite
        try:
            if 'C' not in met_formula and met_formula['O'] == 3 and met_formula['P'] == 1:
                known_mets[PO3] = met
        except KeyError:
            pass
        # Phosphate
        try:
            if 'C' not in met_formula and met_formula['O'] == 4 and met_formula['P'] == 1:
                known_mets[PO4] = met
        except KeyError:
            pass
        # Pyrophosphate
        try:
            if 'C' not in met_formula and met_formula['O'] == 7 and met_formula['P'] == 2:
                known_mets[PPI] = met
        except KeyError:
            pass
        # Water
        try:
            if len(met_formula)==2 and met_formula['O'] == 1 and met_formula['H'] == 2:
                known_mets[H2O] = met
        except KeyError:
            pass
        # Hydrogen peroxide
        try:
            if len(met_formula)==2 and met_formula['O'] == 2 and met_formula['H'] == 2:
                known_mets[H2O2] = met
        except KeyError:
            pass
        # Oxygen
        try:
            if len(met_formula)==1 and met_formula['O'] == 2:
                known_mets[O2] = met
        except KeyError:
            pass
        # Hydroygen
        try:
            if len(met_formula)==1 and met_formula['H'] == 2:
                known_mets[H2] = met
        except KeyError:
            pass
        
        # Nitrogen
        try:
            if len(met_formula)==1 and met_formula['N'] == 2:
                known_mets[N2] = met
        except KeyError:
            pass
        # Proton
        try:
            if len(met_formula)==1 and met_formula['H'] == 1:
                if H in known_mets: # protons in different compartments
                    known_mets[H] += [met]
                else:
                    known_mets[H] = [met]
        except KeyError:
            pass            
        # Metal ions (complexes)
        if len(set(met_formula.keys()).intersection(['Co','Fe','U','Cu',
                'Zn','Mo','Ni','Cr'])):
            metalion.append(met)
    
    
    
    # since here, we define some generic assumptions to find the potential redox pair
    if H in known_mets and H2 not in known_mets:
        'only H2 oxidation results in H+'
        for proton in known_mets[H]:
            potent_pair.remove(proton)
    if H in known_mets and H2 in known_mets:
        # in some cases, H+ is reduced but at the same time we may have proton transfer, e.g.:
        # 'h2_c + 2.0 h_c + mqn7_c --> 2.0 h_e + mql7_c'
        # In this reaction H2 is oxidized to H+, but the other H+ is just transferred
        for proton in known_mets[H]:
            if rxn.get_coefficient(known_mets[H2]) * rxn.get_coefficient(proton) > 0:
                potent_pair.remove(proton)
    if ATP in known_mets and ADP in known_mets:
        'assumption: Phosphorylation od ADP to ATP or vice versa'
        potent_pair.remove(known_mets[ATP])
        potent_pair.remove(known_mets[ADP])
    if COA in known_mets:
        'assumption: CoA is not reduced or oxidized, but just translocated'
        potent_pair.remove(known_mets[COA])
    if PPI in known_mets:
        'assumption: PPi is not reduced or oxidized, but just translocated'
        potent_pair.remove(known_mets[PPI])
    if NH4 in known_mets and (NO2 not in known_mets and NO3 not in known_mets \
                              and N2 not in known_mets):
        'If N2, NO2 and NO3 are not in the reaction, it is a (de)amination'
        potent_pair.remove(known_mets[NH4])
    if PO4 in known_mets and PO3 not in known_mets:
        'If PO3 are not in the reaction, it is a (de)phosphorylation'
        potent_pair.remove(known_mets[PO4])
    if H2O in known_mets and (H2 not in known_mets and O2 not in known_mets \
                              and H2O2 not in known_mets):
        'If neither O2 nor H2 in the reaction, it is a (de)hydration'
        potent_pair.remove(known_mets[H2O])
    if CO2 in known_mets and (CO not in known_mets and CH4 not in known_mets):
        'If neither CO nor CH4 in the reaction, it is an oxidative (de)carboxylation'
        'CO2 is the new oxidized and anything else on the same side of the reaction should be excluded'
        rem_met = [met for met in potent_pair if met.id != known_mets[CO2].id and \
               rxn.get_coefficient(met.id) * rxn.get_coefficient(known_mets[CO2].id) > 0]
        for met in rem_met:
            potent_pair.remove(met)
    
    
    # if still the match is not found, try to find the redox pair based on some well known paris
    if len(potent_pair) > 2:
        # these are positive rules in contrast to the previous rules
        if len(metalion) > 0:
            'If a metal ion is present, probably it acts as the electron donor (acceptor)'
            return metalion
        if CO2 in known_mets and CO in known_mets:
            if rxn.get_coefficient(known_mets[CO].id) * \
                rxn.get_coefficient(known_mets[CO2].id) < 0:
                    return [known_mets[CO2], known_mets[CO]]
        if CO2 in known_mets and CH4 in known_mets:
            if rxn.get_coefficient(known_mets[CH4].id) * \
                rxn.get_coefficient(known_mets[CO2].id) < 0:
                    return [known_mets[CO2], known_mets[CH4]]
        if NH4 in known_mets and N2 in known_mets:
            if rxn.get_coefficient(known_mets[NH4].id) * \
                rxn.get_coefficient(known_mets[N2].id) < 0:
                    return [known_mets[NH4], known_mets[N2]]
        if O2 in known_mets and H2O in known_mets:
            if rxn.get_coefficient(known_mets[H2O].id) * \
                rxn.get_coefficient(known_mets[O2].id) < 0:
                    return [known_mets[O2], known_mets[H2O]]
        if O2 in known_mets and H2O2 in known_mets:
            if rxn.get_coefficient(known_mets[H2O2].id) * \
                rxn.get_coefficient(known_mets[O2].id) < 0:
                    return [known_mets[O2], known_mets[H2O2]]
        if H2O2 in known_mets and H2O in known_mets:
            if rxn.get_coefficient(known_mets[H2O2].id) * \
                rxn.get_coefficient(known_mets[H2O].id) < 0:
                    return [known_mets[H2O2], known_mets[H2O]]
        if NO2 in known_mets:
            reduced = [met for met in potent_pair if 'N' in find_element_coeff(met.formula)\
                       and rxn.get_coefficient(met.id) * rxn.get_coefficient(known_mets[NO2].id) < 0]
            if len(reduced) == 1:
                return [reduced[0], known_mets[NO2]]
        if NO3 in known_mets:
            reduced = [met for met in potent_pair if 'N' in find_element_coeff(met.formula)\
                       and rxn.get_coefficient(met.id) * rxn.get_coefficient(known_mets[NO3].id) < 0]
            if len(reduced) == 1:
                return [reduced[0], known_mets[NO3]]
        if SO3 in known_mets:
            reduced = [met for met in potent_pair if 'S' in find_element_coeff(met.formula)\
                       and rxn.get_coefficient(met.id) * rxn.get_coefficient(known_mets[SO3].id) < 0]
            if len(reduced) == 1:
                return [reduced[0], known_mets[SO3]]
        if SO4 in known_mets:
            reduced = [met for met in potent_pair if 'S' in find_element_coeff(met.formula)\
                       and rxn.get_coefficient(met.id) * rxn.get_coefficient(known_mets[SO4].id) < 0]
            if len(reduced) == 1:
                return [reduced[0], known_mets[SO4]]
    
    
    # Somtimes, it's easy! we just need to check which pair satisfies the redox definition
    # I've put it at the end, as it does not for the loss of H+
    found_pair = check_redox_def(potent_pair)
    if len(found_pair) == 2:
        return found_pair
        
    
    return potent_pair