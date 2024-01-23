#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from ..core.medium import constrain_abiotics

def add_rem_member(new_comm, prev_medium):
    '''
    Adding constraints and variables to the n±1

    Parameters
    ----------
    new_comm : TYPE
        DESCRIPTION.
    prev_medium : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    
    old_medium = []
    for met_id in prev_medium:
        try:
            old_medium += [new_comm.metabolites.get_by_id(met_id)]
        except KeyError: # maybe the new community is smaller
            pass
    
    new_nutrients = []
    for met in new_comm.metabolites:
        if met not in old_medium:
            new_nutrients.append(met)
                
    constrain_abiotics(new_comm, new_nutrients)
    
    return
