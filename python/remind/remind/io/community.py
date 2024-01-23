# -*- coding: utf-8 -*-
"""
.. module:: etfl
   :platform: Unix, Windows
   :synopsis: Expression and thermodynamics-based models

.. moduleauthor:: Pierre Salvy

Make the model serializable
"""
from collections import OrderedDict

from tqdm import tqdm

from cobra import DictList
from cobra.exceptions import SolverNotFound
from optlang.util import  parse_expr
from pytfa.io.dict import get_solver_string, var_to_dict, cons_to_dict, \
    obj_to_dict, rebuild_obj_from_dict
from pytfa.thermo.tmodel import ThermoModel
import cobra.io.dict as cbd


try:
    from pytfa.thermo.tmodel_struct import ThermoModelStructure 
    from pytfa.io.dict import model_from_dict_tfa_struct, model_to_dict_tfa_struct
except ImportError:
    pass
from pytfa.optim.variables import ReactionVariable, ModelVariable
from pytfa.optim.constraints import ReactionConstraint, ModelConstraint
from pytfa.optim.utils import get_all_subclasses
from pytfa.io.dict import model_from_dict as tfa_model_from_dict 
from pytfa.io.dict import model_to_dict as tfa_model_to_dict



from ..core.community import CommModel
from ..optim import constraints, variables



SOLVER_DICT = {
    'optlang.gurobi_interface':'optlang-gurobi',
    'optlang.cplex_interface':'optlang-cplex',
    'optlang.glpk_interface':'optlang-glpk',
}

def make_subclasses_dict(cls):
    the_dict = {x.__name__:x for x in get_all_subclasses(cls)}
    the_dict[cls.__name__] = cls
    return the_dict

REACTION_VARIABLE_SUBCLASSES    = make_subclasses_dict(ReactionVariable)
REACTION_CONSTRAINT_SUBCLASSES  = make_subclasses_dict(ReactionConstraint)
MODEL_VARIABLE_SUBCLASSES       = make_subclasses_dict(ModelVariable)
MODEL_CONSTRAINT_SUBCLASSES     = make_subclasses_dict(ModelConstraint)



def archive_variables(var_dict):
    obj = OrderedDict()

    obj['variables'] = []
    for classname,var in var_dict.items():
        obj[classname] = list(map(var_to_dict, var))

    return obj

def archive_constraints(cons_dict):
    obj = OrderedDict()

    for classname,cons in cons_dict.items():
        obj[classname] = list(map(cons_to_dict, cons))

    return obj

def v_to_dict(variable):
    obj = var_to_dict(variable)
    return obj

def c_to_dict(constraint):
    obj = cons_to_dict(constraint)
    return obj

def species_to_dict(species):
    # this is like a normal pytfa model

    try:
        if isinstance(species,ThermoModelStructure):
            obj = model_to_dict_tfa_struct(species)
            obj['type'] = 'ThermoModelStructure'
    except NameError:
        pass
    if isinstance(species, ThermoModel):
        obj = tfa_model_to_dict(species)
        obj['type'] = 'ThermoModel'


    return obj



def model_to_dict(model):
    """

    :param model:
    :return:
    """

    # Take advantage of cobra's dict serialization for metabolites and
    # reactions
    obj = cbd.model_to_dict(model)
    # obj = OrderedDict()
    obj['name'] = 'Model' if model.name is None else model.name 

    obj['solver'] = get_solver_string(model)
    obj['objective'] = obj_to_dict(model)

    # Copy variables, constraints
    # obj['var_dict'] = archive_variables(model._var_kinds)
    # obj['cons_dict'] = archive_constraints(model._cons_kinds)
    obj['variables'] = list(map(v_to_dict, model._var_dict.values()))
    obj['constraints'] = list(map(c_to_dict, model._cons_dict.values()))

    is_thermo = False
    if isinstance(model, ThermoModel):
        obj['kind'] = 'ThermoModel'
        obj['thermo_data'] = model.thermo_data #it's a dict
        obj['name'] = model.name
        obj['temperature'] = model.TEMPERATURE
        obj['min_ph'] = model.MIN_pH
        obj['max_ph'] = model.MAX_pH
        is_thermo = True

        # Relaxation info
        try:
            obj['relaxation'] = model.relaxation
        except AttributeError:
            pass    

    if isinstance(model, CommModel):

        # Convenience attributes
        obj['kind'] = 'CommModel'
        
        
        # Species
        obj['species'] = list(map(species_to_dict, model.species)) 
        
        obj['environment'] = model.environment
        obj['cooperation'] = model.cooperation

    # Metabolite and Reaction-level cleanup
    for rxn_dict in obj['reactions']:
        rxn = model.reactions.get_by_id(rxn_dict['id'])

        if is_thermo:
            _add_thermo_reaction_info(rxn, rxn_dict)
            
    # Peptides and Thermo
    for met_dict in obj['metabolites']:
        the_met_id = met_dict['id']
        is_peptide = False

        if is_thermo and not is_peptide: # peptides have no thermo
            the_met = model.metabolites.get_by_id(the_met_id)
            _add_thermo_metabolite_info(the_met, rxn_dict)
            met_dict['kind'] = 'Metabolite'

    return obj

def _add_thermo_reaction_info(rxn, rxn_dict):
    if hasattr(rxn, 'thermo'):
        rxn_dict['thermo'] = rxn.thermo

def _add_thermo_metabolite_info(met, met_dict):
    if hasattr(met, 'thermo'):
        met_dict['thermo'] = metabolite_thermo_to_dict(met)
        
def metabolite_thermo_to_dict(metthermo):
    return metthermo.thermo.__dict__


def model_from_dict(obj, solver=None):
    
    # Take advantage of cobra's serialization of mets and reactions
    new = cbd.model_from_dict(obj)

    if solver is not None:
        try:
            new.solver = solver
        except SolverNotFound as snf:
            raise snf
    else:
        try:
            new.solver = obj['solver']
        except KeyError:
            pass
        
    if obj['kind'] == 'CommModel':
        new = CommModel(model=new,
                          name=obj['name'],
                          )
        new.environment = obj['environment']
        new.cooperation = obj['cooperation']
        
    # populate the model
    species_from_dict(new, obj)
    
    
    # Populate variables and constraints
    new._push_queue()
    
    # Force update GPR info
    for rxn in new.reactions:
        rxn.gene_reaction_rule = rxn.gene_reaction_rule

    
    for the_var_dict in tqdm(obj['variables'], desc='rebuilding variables'):
        this_id = the_var_dict['id']
        classname = the_var_dict['kind']
        lb = the_var_dict['lb']
        ub = the_var_dict['ub']
        try: #Backward compat
            scaling_factor = the_var_dict['scaling_factor']
        except KeyError:
            scaling_factor = 1

        if classname in REACTION_VARIABLE_SUBCLASSES:
            hook = new.reactions.get_by_id(this_id)
            this_class = REACTION_VARIABLE_SUBCLASSES[classname]
            nv = new.add_variable(kind=this_class,
                                    hook=hook,
                                    ub=ub,
                                    lb=lb,
                                    scaling_factor=scaling_factor,
                                    queue=True)
        
            
        elif classname in MODEL_VARIABLE_SUBCLASSES:
            hook = new
            this_class = MODEL_VARIABLE_SUBCLASSES[classname]
            nv = new.add_variable(kind=this_class,
                                    hook=hook,
                                    id_=this_id,
                                    ub=ub,
                                    lb=lb,
                                    scaling_factor=scaling_factor,
                                    queue=True)
           
                
        else:
            raise TypeError(
                'Class {} serialization not handled yet' \
                    .format(classname))

    new._push_queue()

    variable_parse_dict = {x.name:x for x in new.variables}

    for the_cons_dict in tqdm(obj['constraints'], desc='rebuilding constraints'):
        this_id = the_cons_dict['id']
        classname = the_cons_dict['kind']
        new_expr = parse_expr(the_cons_dict['expression'],
                              local_dict = variable_parse_dict)
        
        lb = the_cons_dict['lb']
        ub = the_cons_dict['ub']

        if classname in REACTION_CONSTRAINT_SUBCLASSES:
            hook = new.reactions.get_by_id(this_id)
            this_class = REACTION_CONSTRAINT_SUBCLASSES[classname]
            nc = new.add_constraint(kind=this_class, hook=hook,
                                      expr=new_expr,
                                      ub=ub,
                                      lb=lb,
                                      queue=True)
        elif classname in MODEL_CONSTRAINT_SUBCLASSES:
            hook = new
            this_class = MODEL_CONSTRAINT_SUBCLASSES[classname]
            nc = new.add_constraint(kind=this_class, hook=hook,
                                      expr=new_expr, id_=this_id,
                                      ub=ub,
                                      lb=lb,
                                      queue=True)
        else:
            raise TypeError('Class {} serialization not handled yet' \
                            .format(classname))

    try:
        rebuild_obj_from_dict(new, obj['objective'])
    except KeyError:
        pass

    
    new.repair()
    
    
    return new


def species_from_dict(new, obj):
    
    for sp in obj['species']:
        try:
            if sp['type'] == 'ThermoModelStructure':
                raw = model_from_dict_tfa_struct(sp)
        except NameError:
            pass    
            
        if sp['type'] == 'ThermoModel':
            raw = tfa_model_from_dict(sp)

        # raw = tfa_model_from_dict(sp)
        # if isinstance()
        species = raw
        new.species += [species]
        species._model = new
        
    return



import json

from pytfa.io.json import MyEncoder, check_json_extension

def save_json_model(model, filepath):
    """
    Saves the model as a JSON file

    :param model:
    :param filepath:
    :return:
    """

    filepath = check_json_extension(filepath)
    obj = model_to_dict(model)

    with open(filepath, 'w') as fid:
        json.dump(obj, fid, cls=MyEncoder)


def load_json_model(filepath, solver = None):
    """
    Loads a model from a JSON file

    :param filepath:
    :param solver:
    :return:
    """

    filepath = check_json_extension(filepath)
    with open(filepath, 'r') as fid:
        obj = json.load(fid)

    model = model_from_dict(obj, solver=solver)
    return model