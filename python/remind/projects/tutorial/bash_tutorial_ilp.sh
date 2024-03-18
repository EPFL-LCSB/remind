#!/bin/bash


#to run for the tutorial script with 2 models from the bee gut community

#for positive interaction cross-feeding patterns
#first argument is the objective function and the second is the alternation with respect
#to which the integer cuts are added
#max cooperation
python run_ilp_tutorial_community_model.py 2 2 #patterns for metabolites
python run_ilp_tutorial_community_model.py  2 11 #patterns with directional
python run_ilp_tutorial_community_model.py  2 12 #patterns directional and according to biomass yields


#min competition
python run_ilp_tutorial_community_model.py 1 3 #patterns for metabolites
python run_ilp_tutorial_community_model.py  1 13 #patterns with directional
python run_ilp_tutorial_community_model.py  1 14 #patterns directional and according to biomass yields

#min nutritional requirements
python run_ilp_tutorial_community_model.py 6 4 #patterns for metabolites in min environment
python run_ilp_tutorial_community_model.py  6 9 #patterns w metabolites and yield

