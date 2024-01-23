# -*- coding: utf-8 -*-
"""
.. module:: etfl
   :platform: Unix, Windows
   :synopsis: Expression and thermodynamics-based models

.. moduleauthor:: Pierre Salvy

Make the model serializable
"""
from collections import OrderedDict, defaultdict

from tqdm import tqdm

from cobra import DictList
from cobra.exceptions import SolverNotFound
from optlang.util import expr_to_json, parse_expr
from pytfa.io.dict import get_solver_string, var_to_dict, cons_to_dict, \
    obj_to_dict, rebuild_obj_from_dict
from pytfa.thermo.tmodel import ThermoModel

try:
    from pytfa.thermo.tmodel_struct import ThermoModelStructure 
    from pytfa.io.dict import model_from_dict_tfa_struct, model_to_dict_tfa_struct
except ImportError:
    pass
from pytfa.optim.variables import ReactionVariable, ModelVariable, MetaboliteVariable
from pytfa.optim.constraints import ReactionConstraint, ModelConstraint, MetaboliteConstraint
from pytfa.optim.utils import get_all_subclasses
from pytfa.io.dict import model_from_dict as tfa_model_from_dict 
from pytfa.io.dict import model_to_dict as tfa_model_to_dict

from ..core.objects import CommunityMetabolite, CommunityReaction, CommunitySpecies, \
            DIME

from ..core.interaction import InteractionModel
from ..optim import constraints, variables

from ..optim.variables import SecretionActivation, UptakeActivation


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
METABOLITE_VARIABLE_SUBCLASS    = make_subclasses_dict(MetaboliteVariable)
METABOLITE_CONSTRAINT_SUBCLASS  = make_subclasses_dict(MetaboliteConstraint)



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
    try:
        obj['species'] = variable.species
    except AttributeError:
        obj['species'] = ''
    return obj

def c_to_dict(constraint):
    obj = cons_to_dict(constraint)
    try:
        obj['species'] = constraint.species
    except AttributeError:
        obj['species'] = ''
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


    # obj['elements'] = species.elements
    obj['yields'] = [str(x) for x in species.yields]
    obj['dimes'] = [x.id for x in species.dimes]
    obj['substrates'] = [x.id for x in species.substrates]
    obj['products'] = [x.id for x in species.products]
    obj['rxns'] = [x.id for x in species.rxns]

    return obj

def dime_to_dict(dime):
    obj = OrderedDict()
    obj['id'] = dime.id
    obj['species'] = dime.species
    # obj['element'] = dime.element
    obj['yield'] = str(dime._yield)
    obj['substrates'] = [x.id for x in dime.substrates]
    obj['products'] = [x.id for x in dime.products]
    obj['rxns'] = [x.id for x in dime.rxns]
    obj['varname'] = dime.variable.name
    return obj

def metabolite_to_dict(metabolite):
    obj = OrderedDict()
    obj['id'] = metabolite.id
    obj['name'] = metabolite.name
    obj['species'] = [x for x in metabolite.species]
    obj['rxns'] = [r.id for r in metabolite.rxns]
    return obj

def reaction_to_dict(reaction):
    obj = OrderedDict()
    obj['id'] = reaction.id
    obj['name'] = reaction.name
    obj['species'] = reaction.species
    obj['metabolite'] = reaction.metabolite.id
    obj['dimes'] = [x.id for x in reaction.dimes]
    try:
        obj['upt_varname'] = reaction.upt_variable.name
    except AttributeError:
        pass
    try:
        obj['sec_varname'] = reaction.sec_variable.name
    except AttributeError:
        pass
    return obj

def model_to_dict(model):
    """

    :param model:
    :return:
    """

    # Take advantage of cobra's dict serialization for metabolites and
    # reactions
    # obj = cbd.model_to_dict(model)
    obj = OrderedDict()

    obj['solver'] = get_solver_string(model)
    obj['objective'] = obj_to_dict(model)

    # Copy variables, constraints
    # obj['var_dict'] = archive_variables(model._var_kinds)
    # obj['cons_dict'] = archive_constraints(model._cons_kinds)
    obj['variables'] = list(map(v_to_dict, model._var_dict.values()))
    obj['constraints'] = list(map(c_to_dict, model._cons_dict.values()))


    if isinstance(model, InteractionModel):

        # Convenience attributes
        obj['kind'] = 'InteractionModel'
        
        # Metabolites
        obj['metabolites'] = list(map(metabolite_to_dict, model.metabolites))

        # Reactions
        obj['reactions'] = list(map(reaction_to_dict, model.reactions))
        
        # Species
        obj['species'] = list(map(species_to_dict, model.species))

        # DiMEs
        obj['dimes'] = list(map(dime_to_dict, model.dimes))
        

    return obj




def model_from_dict(obj, solver=None):
    
    new = InteractionModel()

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

    # populate the model
    species_from_dict(new, obj)
    dimes_from_dict(new, obj)
    reactions_from_dict(new, obj)
    metabolites_from_dict(new, obj)
    reconstruct_connections(new, obj)
    
    # Populate variables and constraints
    new._push_queue()

    
    for the_var_dict in tqdm(obj['variables'], desc='rebuilding variables'):
        this_id = the_var_dict['id']
        classname = the_var_dict['kind']
        lb = the_var_dict['lb']
        ub = the_var_dict['ub']
        species = the_var_dict['species']
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
        elif classname in METABOLITE_VARIABLE_SUBCLASS:
            hook = new.metabolites.get_by_id(this_id)
            this_class = METABOLITE_VARIABLE_SUBCLASS[classname]
            nv = new.add_variable(kind=this_class,
                                    hook=hook,
                                    ub=ub,
                                    lb=lb,
                                    queue=True)
            
        elif classname in MODEL_VARIABLE_SUBCLASSES:
            hook = new
            this_class = MODEL_VARIABLE_SUBCLASSES[classname]
            try:
                nv = new.add_variable(kind=this_class,
                                    hook=hook,
                                    id_=this_id,
                                    ub=ub,
                                    lb=lb,
                                    scaling_factor=scaling_factor,
                                    queue=True)
            except TypeError:
                nv = new.add_variable(kind=this_class,
                                    hook=hook,
                                    id_=this_id,
                                    ub=ub,
                                    lb=lb,
                                    species = species,
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
        species = the_cons_dict['species']

        if classname in REACTION_CONSTRAINT_SUBCLASSES:
            # in case that we have ReactionActivation constraint the id is not the id of the reaction anymore
            if classname == 'ReactionActivation': # the id looks like this 'Shewanella_EX_hdca_e_Shewanella_0.2_Carbon_1'
                dime_id = species + this_id.split('_' + species)[1]
                dime = new.dimes.get_by_id(dime_id)
                this_id = this_id.split('_' + species)[0]
                hook = new.reactions.get_by_id(this_id)
                this_class = REACTION_CONSTRAINT_SUBCLASSES[classname]
                nc = new.add_constraint(kind=this_class, hook=hook,
                                          dime = dime,
                                          expr=new_expr,
                                          ub=ub,
                                          lb=lb,
                                          queue=True)
            else:
                hook = new.reactions.get_by_id(this_id)
                this_class = REACTION_CONSTRAINT_SUBCLASSES[classname]
                nc = new.add_constraint(kind=this_class, hook=hook,
                                          expr=new_expr,
                                          ub=ub,
                                          lb=lb,
                                          queue=True)
        elif classname in METABOLITE_CONSTRAINT_SUBCLASS:
            hook = new.metabolites.get_by_id(this_id)
            this_class = METABOLITE_CONSTRAINT_SUBCLASS[classname]
            nc = new.add_constraint(kind=this_class, hook=hook,
                                      expr=new_expr,
                                      ub=ub,
                                      lb=lb,
                                      queue=True)
            
        elif classname in MODEL_CONSTRAINT_SUBCLASSES:
            hook = new
            this_class = MODEL_CONSTRAINT_SUBCLASSES[classname]
            try:
                nc = new.add_constraint(kind=this_class, hook=hook,
                                      expr=new_expr, id_=this_id,
                                      ub=ub,
                                      lb=lb,
                                      queue=True)
            except TypeError:
                nc = new.add_constraint(kind=this_class, hook=hook,
                                      expr=new_expr, id_=this_id,
                                      ub=ub,
                                      lb=lb,
                                      species = species,
                                      queue=True)
        else:
            raise TypeError('Class {} serialization not handled yet' \
                            .format(classname))

    try:
        rebuild_obj_from_dict(new, obj['objective'])
    except KeyError:
        pass

    
    new.repair()
    
    for dime in new.dimes:
        dime.init_variable()
        
    sec_act = new.get_variables_of_type(SecretionActivation)
    upt_act = new.get_variables_of_type(UptakeActivation)
    
    for rxn in new.reactions:
        try:
            var = sec_act.get_by_id(rxn.id) # does this variable exists
            rxn.init_sec_variable()
        except KeyError:
            pass
        try:
            var = upt_act.get_by_id(rxn.id) # does this variable exists
            rxn.init_upt_variable()
        except KeyError:
            pass
            
    
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
        species = CommunitySpecies(raw)
        new.species += [species]
        species._model = new

        
    return

def dimes_from_dict(new, obj):
    
    for dime in obj['dimes']:
        d = DIME(id = dime['id'],
                          species_id = dime['species'],
                          rxn_set = set(dime['rxns']), #ids for now
                          # element = dime['element'],
                          yield_ = dime['yield'])
        new.dimes += [d]
        d._model = new
        new.species.get_by_any(dime['species'])[0].dimes += [d]
        
    return

def reactions_from_dict(new, obj):
    
    for rx in obj['reactions']:
        rxn = CommunityReaction(id = rx['id'],
                                name = rx['name'],
                                metabolite = rx['metabolite'], # id for now
                                species_id = rx['species'],
                                dimes = DictList(),
                                )
        # DiMEs are already added
        for di in rx['dimes']:
            rxn.dimes += [new.dimes.get_by_id(di)]
            
        new.reactions += [rxn]
        rxn._model = new
        
    return

def metabolites_from_dict(new, obj):
    
    for m in obj['metabolites']:
        met = CommunityMetabolite(id = m['id'],
                                  name = m['name'],
                                  species = set(m['species']),
                                  rxns = set(
                                      [new.reactions.get_by_id(x) for x in m['rxns']]
                                       ))
        met._model = new
        new.metabolites += [met]
    
    return

def reconstruct_connections(new, obj):
    
    for rxn in new.reactions:
        met = new.metabolites.get_by_id(rxn.metabolite)
        rxn.metabolite = met
        
    for d in obj['dimes']:
        dime = new.dimes.get_by_id(d['id'])
        rxn_set = set([new.reactions.get_by_id(x) for x in dime.rxns])
        dime.rxns = rxn_set
        dime.substrates += [new.metabolites.get_by_id(x) for x in d['substrates']]
        dime.products += [new.metabolites.get_by_id(x) for x in d['products']]
        
    for sp in obj['species']:
        species = new.species.get_by_any(sp['id'])[0]
        species.rxns += [new.reactions.get_by_id(x) for x in sp['rxns']]
        species.substrates += [new.metabolites.get_by_id(x) for x in sp['substrates']]
        species.products += [new.metabolites.get_by_id(x) for x in sp['products']]
        
    
    return

