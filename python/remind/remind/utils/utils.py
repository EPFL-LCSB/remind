import numpy as np
import itertools
from tqdm import tqdm

from ..optim.variables import UptakeActivation, PositiveInteraction


def batch_growth_from_abund(abundance_dict, community_growth, time):
    '''
    This function tries to capture the growth rate of single species based on their
    relative abundance and the growth rate of the community.
    Two main assumtions:
        1) the initial abundance of all species is the same
        2) the growth rate of the community is the rate of increase in total biomass.

    If the relative abundance of the first species is a, the abundance of the other
    will be 1-a (for two organisms). We have the following two equations based on
    the two abovementined assumptions:
        1) log((1-a)/a) = (mu_2 - mu_1) * time
        2) mu = a*mu_1 + (1-a)*mu_2

    Parameters
    ----------
    abundance_dict : dict
        A dictionary where the keys are organism names and values are their rel. abundance 
    community_growth : float
        A value for the growth rate of the community
    time: float (h)
        A value to show at which time point the abundances have been obtained

    Returns
    -------
    growth_dict : dict
        A dictionary where the keys are organism names and values are their growth rate.

    '''
    
    
    tot_rel_abund = sum([k for k in abundance_dict.values()])
    if not np.isclose(tot_rel_abund, 1):
        raise ValueError('The values for the relative abundance are not acceptable.'
                         'Make sure that no member of the community is misssing.')
        
    # it is a system of n variables and n equations
    rhs = [0] * len(abundance_dict) # initiating the right hand side vector    
    # the first equation is: mu = a*mu_1 + (1-a)*mu_2
    com_eq = [0] * len(abundance_dict)
    for ind, val in enumerate(abundance_dict.items()):
        com_eq[ind] = val[1]
    rhs[0] = community_growth
    
    # log((1-a)/a)/t = mu_2 - mu_1 --> n-1 equations    
    for ind, val in enumerate(abundance_dict.items()):
        if ind == 0: 
            ref_abund = val # define all equations wrt the first organism
    
    coeff_matrix = []
    for ind, val in enumerate(abundance_dict.items()):
        if ind != 0:
          recurring_coeff = [0] * len(abundance_dict)
          recurring_coeff[0] = -1 # define all equations wrt the first organism
          recurring_coeff[ind] = 1
          rhs[ind] = np.log(val[1] / ref_abund[1]) / time
          coeff_matrix += [recurring_coeff]
          
    coeff_matrix = np.array(coeff_matrix)
    
    
    
    # concatening the last equations to the rest
    coeff_matrix = np.concatenate((np.array([com_eq]), coeff_matrix), axis=0)
    
    # the RHS vector is zero for all n-1 equations, except the last one
    rhs_vec = np.array(rhs)
    
    growth_vec = np.linalg.solve(coeff_matrix, rhs_vec)
    
    growth_dict = {}
    for ind, val in enumerate(abundance_dict.items()):
        growth_dict[val[0]] = growth_vec[ind]
    
    return growth_dict

def chemostat_growth_from_abund(species_list, community_growth):
    '''
    This function tries to capture the growth rate of single species based on their
    relative abundance and the growth rate of the community assuming that system
    has converged to a stable steady-state. For this, the composition of the system
    is independent of time and hence, the initial abundances is similar to the current
    abundances.
    Taking these assumtpions, it can be deduced that the growth rate of all species
    are identical and equal to the growth rate of the community.

    Parameters
    ----------
    species_list: list
        A list of organism names to be used to create the dictionary of the output 
    community_growth : float
        A value for the growth rate of the community

    Returns
    -------
    growth_dict : dict
        A dictionary where the keys are organism names and values are their growth rate.

    '''
    
    # simple to implement as the growth rates are the same
    growth_dict = {species:community_growth for species in species_list}
    
    return growth_dict


def find_chemical_similarity(query, base, sigma=1):
    
    '''
    A fundction to find similarity in chemical formula between a compund and a list of compunds
    
    Parameters
    ----------
    query: dict
        key is the metabolite ID and value is decomposed fomula (see find_element_coeff(met_formula))
    base: dict
        key is the metabolite ID and value is decomposed fomula (see find_element_coeff(met_formula))
    sigma: int
        number of difference in oxygens that is equivalent to unit difference in carbons
    
    Returns
    -------
    score: float
        
    '''
    
    # we want to check what fraction of elements is common
    list1 = [x for x in query.keys()]
    list2 = [x for x in base.keys()]
    N = len(set(list1).intersection(list2))
    D = len(set(list1).union(list2))
    
    if N == 0:
        return 100
    else:
         elementwise = D-N
         
    # specifically the coefficents of O and C are important
    try:
        C_score = abs(query['C'] - base['C'])
        O_score = abs(query['O'] - base['O'])/(sigma+1)
    except KeyError:
        return 10
    
    # for the other elements compare the numbers if possible
    other_el_score = 0
    for e in set(list1).union(list2):
        if e =='C' or e == 'O': # these already considered
            pass
        else:
            try:
                other_el_score += abs(query[e] - base[e])/(sigma)
            except KeyError:
                pass
    
    # difference in Carbon is more important than Oxygen
    score = elementwise + C_score + O_score + other_el_score
    
    
    return score

def fix_dir_by_thermo(model, tfa_model, growth_lb, growth_id, bigM=1000):
    '''
    

    Parameters
    ----------
    model : pytfa.Model
        a model without thermodynamic constraints.
    tfa_model : pytfa.Model
        a model with thermodyanamic constraints.
    growth_lb : float
        the lower bound for the growth.
    growth_id : str

    Returns
    -------
    None.

    '''
    
    tfa_model.reactions.get_by_id(growth_id).lower_bound = growth_lb
    
    for rxn in tqdm(tfa_model.reactions, desc='Doing TVA to find directionalities:'):
        tfa_model.objective = rxn
        tfa_model.objective_direction = 'min'
        min_flux = tfa_model.slim_optimize()
        tfa_model.objective_direction = 'max'
        max_flux = tfa_model.slim_optimize()
        q_rxn = model.reactions.get_by_id(rxn.id)
        if q_rxn.lower_bound == -bigM and min_flux >= 0:
            q_rxn.lower_bound = 0
        if q_rxn.upper_bound == bigM and max_flux <= 0:
            q_rxn.upper_bound = 0
    
    
    return


def find_environment(model, solution):
    '''
    A function to find abiotic metabolic IDs.

    Parameters
    ----------
    model : TYPE
        DESCRIPTION.
    solution : TYPE
        DESCRIPTION.

    Returns
    -------
    met_ids : TYPE
        DESCRIPTION.

    '''
    
    uptakes = model.get_variables_of_type(UptakeActivation)
    pos_int = model.get_variables_of_type(PositiveInteraction)
    
    # first, find all the metabolites for which positive interactions exist
    pos_int_mets = [var.metabolite for var in pos_int if \
                    np.isclose(solution.loc[var.name][0],1)]
        
    # find the uptakes for which we do not have positive interaction, as aboitic resources
    met_ids = []
    for var in uptakes:
       if np.isclose(solution.loc[var.name][0],1):
            if var.reaction.metabolite not in pos_int_mets:
                if var.reaction.metabolite.id not in met_ids: # not to add replicates for competitions
                    met_ids.append(var.reaction.metabolite.id)
    
    return met_ids
    