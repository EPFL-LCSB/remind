"""

@author: ReMIND Team
"""

from remind.io.json import load_json_model
from remind.optim.variables import PositiveInteraction, NegativeInteraction, \
    UptakeActivation, SecretionActivation, YieldUse
from remind.optim.alternatives import find_alternative_solutions

from pytfa.optim.utils import symbol_sum
#
import os
from sys import argv
from remind.core.medium import constrain_abiotics
from pathlib import Path


def get_base_path():
    """
    Returns the base path of the installed remind package if importable.
    Otherwise returns an empty string (or None).
    """
    try:
        import remind
        base = Path(remind.__file__).resolve().parent
        return str(base.parents[1])
    except ImportError:
        # Not installed (e.g., inside docker where local tree is used)
        return ""   # or return None
base=get_base_path()

#default it is 0 if only wanted to stay in optimal solutions

comm_name = '2_member_bee'

#different alternation options basically different yield cuts
alternation_opts=['reaction','interaction','interaction_positive','interaction_negative','abiotic_sources',\
                  'reaction_yield','interaction_yield','interaction_positive_yield','interaction_negative_yield','abiotic_sources_yield',
                  'interaction_positive_directional_2','interaction_positive_directional_general','interaction_positive_directional_general_yield',
                  'interaction_negative_directional_general','interaction_negative_directional_general_yield']

max_alternative=5000

_, objective_no,alt=argv
obj_num=int(objective_no)
alternation=alternation_opts[int(alt)]


filepath=(base+"/remind/projects/tutorial/ilp_model_tutorial/model_2_member_180324.json")
model= load_json_model(filepath)



### setting different objective functions for the community
# the number of competitions and cooperations in the community
pos_int = model.get_variables_of_type(PositiveInteraction)
neg_int = model.get_variables_of_type(NegativeInteraction)
competitions = symbol_sum([v for v in neg_int])
cooperations = symbol_sum([v for v in pos_int])

# total number of reactions
num_upt = symbol_sum(model.get_variables_of_type(UptakeActivation))
num_sec = symbol_sum(model.get_variables_of_type(SecretionActivation))


#minimize competitive interactions
if obj_num == 1:
    model.objective = competitions
    model.objective_direction = 'min'

#maximize cooperative interactions
if obj_num == 2:
    model.objective = cooperations
    model.objective_direction = 'max'

#minimal number of metabolites as nutritional requirements
if obj_num==6: #('minimize abiotic nutrients for two species')

    constrain_abiotics(model)


sol = model.optimize()
solution = sol.objective_value


#this is the stopping condition if it should only generate optimal alternatives
#or should also look for suboptimal solutions
#if you put zero it will only look for optimal solutions
stop_crit = int(solution - 1) #in case you want to get it until all cross-feeding
stop_crit=0
print("MAX COOP NUMBER IS {}".format(solution))


# output file to save the results
data_path = base+'/remind/projects/tutorial/bee_tutorial_ilp_solutions_{}/obj_num_{}_alternatives_{}'.format(comm_name, obj_num,
                                                                                                 alternation)
if not os.path.exists(data_path):
    os.makedirs(data_path)


growth_ids = {
    'Gapicola': 'Growth',
    'Salvi': 'Growth',

}

growth_limits = {
    'Gapicola': 0.2,
    'Salvi': 0.2,
}


sols_all,sols_active = find_alternative_solutions(model, growth_ids, growth_limits,
                                  alternation=alternation, max_alts=max_alternative,check_feasibility = False,stop_cond=True,criteria=stop_crit)

altsno=sols_active.alternative.nunique()

sols_all.to_hdf(data_path+'/sols_all_alternative_for_diff_obj_{}_alternation_{}_altno_{}_.h5'.format(obj_num,alternation,altsno),key='s')


sols_active.to_hdf(data_path+'/sols_active_alternative_for_diff_obj_{}_alternation_{}_altno_{}_.h5'.format(obj_num,alternation,altsno),key='s')

#

