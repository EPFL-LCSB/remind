#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 17:18:31 2022

@author: ReMIND Team
"""

import pandas as pd
from remind.utils.postprocessing import remove_redundant_sols

organism = 'Geobacter'
# organism = 'Rhodoferax'
# organism = 'Shewanella'

yield_cuts = [1, 1.111, 1.25, 1.42, 1.66, 2, 2.5, 3.333, 5, 10]

input_path = '../projects/ime_alternatives/imes_toy/DIME_{}_alternative_mets_for_growth_0.2_{}.csv'
output_path = '../projects/ime_alternatives/imes_toy/DIME_{}_alternative_mets.csv'

# importing and concatenating the DiMEs for different yields
dime_list = [] 
for y in yield_cuts:
    dime_list += [pd.read_csv(input_path.format(organism,y))]
    
df = pd.concat(dime_list)

data  = remove_redundant_sols(df)
data.to_csv(output_path.format(organism))
    
    