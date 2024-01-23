#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 15:14:21 2023

@author: ReMIND Team
"""
import pandas as pd
import numpy as np
import os
from os.path import join as pjoin
from pandas import ExcelWriter



load_dir = 'ilp_solutions_Geobacter_Rhodoferax/_engineering/directional'

rep = ' // '
NH4 = 'nh4_e'
Fe3 = 'fe3_e'

# columns = ['Negative_Interaction',	'Positive_Interaction',	'Geobacter_Uptakes',
#            'Rhodoferax_Uptakes',	'Geobacter_Secretions',	'Rhodoferax_Secretions']

not_needed_int = ['so4_e', 'pi_e', 
                   'h_e', 'h2o_e'
                  ]
not_needed_exh = ['EX_so4_e' , 'EX_pi_e' , 'EX_mg2_e' , 'EX_ca2_e' , 'EX_k_e', 'EX_h_e', 'EX_h2o_e']

results = {}
writer = ExcelWriter('ilp_solutions_Geobacter_Rhodoferax/_engineering/engineering.xlsx')


for path in os.listdir(load_dir):
    if '.csv' in path: # only csv files to be read
        results = {}
        path_dir = pjoin(load_dir, path)
        data = pd.read_csv(path_dir, index_col = 0)
        for ind, df in data.iterrows():
            group_1 = rep.join([x for x in df['Negative_Interaction'].split(rep) if \
                                x not in not_needed_int])# competitions
            group_2 = float(df['Geobacter_Yield'])
            group_3 = float(df['Rhodoferax_Yield'])
            
            try:
                group_4 = rep.join([x for x in df['Positive_Interaction'].split(rep) if \
                                    'EX_'+x in df['Rhodoferax_Uptakes'].split(rep) \
                                        and 'EX_'+x not in not_needed_exh]) # From Geobacter to Rhodoferax
                group_5 = rep.join([x for x in df['Positive_Interaction'].split(rep) if \
                                    'EX_'+x in df['Geobacter_Uptakes'].split(rep) \
                                    and 'EX_'+x not in not_needed_exh]) # From Rhodoferax  to Geobacter
                
            except AttributeError: # it is empty
                group_4 = ''
                group_5 = ''
            
            order = 0
            # further ordering 
            c_comp = group_1.split(rep)
            if NH4 in c_comp: c_comp.remove(NH4)
            if Fe3 in c_comp: c_comp.remove(Fe3)
            if len(c_comp) > 0: # competeing over a carbon source
                if 'lac-L_e' in c_comp:
                    order -= 0.7 # further penalty
                elif 'mal-L_e' in c_comp:
                    order -= 0.6 # further penalty
                elif 'cit_e' in c_comp:
                    order -= 0.5 # further penalty
                    
            c_coop_1 = group_5.split(rep)
            # if NH4 in c_comp: c_coop_1.remove(NH4) # no cooperation for amonium
            if Fe3 in c_coop_1: c_coop_1.remove(Fe3)
            if len(c_coop_1) > 0: # competeing over a carbon source
                if 'mal-L_e' in c_coop_1:
                    order -= 0.4 # further penalty
                elif 'ac_e' in c_coop_1:
                    order -= 0.3 # further penalty
                    
            c_coop_2 = group_4.split(rep)
            # if NH4 in c_comp: c_coop_1.remove(NH4) # no cooperation for amonium
            if Fe3 in c_coop_2: c_coop_2.remove(Fe3)
            if len(c_coop_2) > 0: # competeing over a carbon source
                if 'mal-L_e' in c_coop_2:
                    order -= 0.2 # further penalty
                elif 'ac_e' in c_coop_2:
                    order -= 0.1 # further penalty
            
            results[ind] = {'Comp': group_1,
                            'G->R': group_4,
                            'R->G': group_5,
                            'Geo_Y' : group_2,
                            'Rhod_Y': group_3,
                            'order': order
                            }
        processed = pd.DataFrame.from_dict(results, orient='index')
        processed = processed.sort_values(['order'])
        processed = processed.sort_values(['Geo_Y'], ascending=False)
        processed.to_excel(writer, sheet_name=path.replace('_interaction_directional.csv', '').replace('x_Geobacter_Rhodoferax_',''))
writer.save()               