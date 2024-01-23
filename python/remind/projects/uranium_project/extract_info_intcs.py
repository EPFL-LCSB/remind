#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
from os.path import join as pjoin


load_dir = 'ilp_solutions'


organisms = ['Geobacter','Rhodoferax',
             'Shewanella'
             ]

neg_int_mod = 'NI__'
pos_int_mod = 'PI__'
upt_mod     = 'UA_{}_'
sec_mod     = 'SA_{}_'
yield_mod   = 'YU_{}_'


for path in os.listdir(load_dir):
    if '.csv' in path: # only csv files to be read
    
        path_dir = pjoin(load_dir, path)
        data = pd.read_csv(path_dir, index_col = 0)
        data_dict = data.to_dict()
        
        
        store = []
        for _, alt in data_dict.items():
            
            df = pd.DataFrame(columns=['Negative_Interaction', 'Positive_Interaction'] + \
                          ['{}_Uptakes'.format(x) for x in organisms] + \
                          ['{}_Secretions'.format(x) for x in organisms] + \
                          ['{}_Yield'.format(x) for x in organisms] + \
                              ['Objective_value'])
            
            df['Objective_value'] = [alt['objective']]
                
            df['Negative_Interaction'] = [' // '.join([k.replace(neg_int_mod,'') for k,v \
                                                      in alt.items() \
                                                      if neg_int_mod in k and np.isclose(v, 1)])]
            df['Positive_Interaction'] = [' // '.join([k.replace(pos_int_mod,'') for k,v \
                                                      in alt.items() \
                                                      if pos_int_mod in k and np.isclose(v, 1)])]
            for org in organisms:
                org_upt_mod = upt_mod.format(org)
                org_sec_mod = sec_mod.format(org)
                org_yld_mod = yield_mod.format(org)
                
                df['{}_Uptakes'.format(org)] = [' // '.join([k.replace(org_upt_mod,'') for k,v \
                                                      in alt.items() \
                                                      if org_upt_mod in k and np.isclose(v, 1)])]
                df['{}_Secretions'.format(org)] = [' // '.join([k.replace(org_sec_mod,'') for k,v \
                                                      in alt.items() \
                                                      if org_sec_mod in k and np.isclose(v, 1)])]
                df['{}_Yield'.format(org)] = [' // '.join([k.replace(org_yld_mod,'') for k,v \
                                                      in alt.items() \
                                                      if org_yld_mod in k and np.isclose(v, 1)])]
            store += [df]
        
        frame = pd.concat(store, ignore_index=True)
        frame.to_csv(path_dir.replace('ilp_solutions/','ilp_solutions/x_')) # add an x to the beginning of the name
