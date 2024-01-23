#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .model import Model
# from ..io.dict import model_from_dict, model_to_dict
from ..optim.constraints import DIMEUse, ReactionActivation, UptakeInactivation, SecretionInactivation, \
     ActiveNegInteraction, ActivePosInteraction, InactiveNegInteraction, InactiveProdUse,\
         InactiveSubsUse, ActiveProdUse, ActiveSubsUse, SimultUse, YieldChoose, BioticProduction
from ..optim.variables import UptakeActivation, SecretionActivation, PositiveInteraction, \
    NegativeInteraction, SubstrateUse, ProductUse, YieldUse
from .objects import CommunityMetabolite

from pytfa.optim.utils import copy_solver_configuration, symbol_sum
from pytfa.utils.logger import get_bistream_logger
from cobra.core.dictlist import DictList
from tqdm import tqdm
from pytfa.optim.utils import get_all_subclasses
from ..optim.constraints import SpeciesConstraint,CommunityReactionConstraint
from ..optim.variables import SpeciesVariable, CommunityReactionVariable

UPTAKES = 'uptakes'
SECRETIONS = 'secretion'


def get_met_id(rxn_id, species_id):
    '''
    The reactions are exchange and each are associated with a single metabolite.
    The reaction ID is : 'species id' + '_' + 'EX_' + 'metabolite id'
    Extract metabolite id out of it!

    Parameters
    ----------
    rxn_id : str
    species_id : str

    Returns
    -------
    met_id : str

    '''
    
    if 'EX_' in rxn_id:
        met_id = rxn_id.replace('EX_', '')
        
    met_id = met_id.replace(species_id, '')
    
    return met_id

class InteractionModel(Model):
    """ A class representing the community model"""
    
    def __init__(self, model=None, name=None):
        Model.__init__(self, model, name)
        
        self.dimes = DictList()
        self.species = DictList()
        self.logger = get_bistream_logger('Community model' + str(self.name))
        
    def add_species(self, species_list, **kwargs):
        """
        The function to define community constraints and variables. 
        This adds the community members as well as their nutritional 
        requiremnet, i.e. their DIMEs. Then, it adds the reactions for
        each DIME. Finally, it adds the metabolites for each reaction.

        Parameters
        ----------
        species_list : iterable
            list of core.CommunitySpecies.

        Returns
        -------
        None.

        """
        #todo check for generalization for element again
        self.species += [x for x in species_list]
        self.repair()
        
        for species in species_list:
            species._model = self
            self.dimes += species.dimes
            rxn_set = species.rxns
            
            # Also add the constraint for only one active DIME per each element per each yield
            var_dict = dict()

            for dime in tqdm(species.dimes,desc='Adding DIMEs and reactions for {}'.format(species)):
                yield_ = dime._yield
                
                dime._model = self
                dime.init_variable()
                try:
                    var_dict['{}_{}'.format(species.id,yield_)] += [dime.variable]
                except KeyError: # to see this ID for the first time
                    var_dict['{}_{}'.format(species.id,yield_)] = [dime.variable]

                # add the Actiavation constraint
                self._add_rxns_from_dime(dime)
                
            yield_use_vars = dict()
            for id_, this_dime_vars in var_dict.items():
                expr = symbol_sum(this_dime_vars)
                # RHS is a binary variable associated to the yield
                # extract the yield for this set of iMEs
                this_yield = id_.split('_')[1] # id_ is of the following format species_id_yield_element
                try:
                    var = yield_use_vars[this_yield]
                except KeyError: # This is the first time that we see this yield
                    var = self.add_variable(YieldUse,
                                    hook = self,
                                    species = species.id,
                                    id_ = '{}_{}'.format(species.id,this_yield),
                                    lb = 0,
                                    ub = 1)
                                    
                    yield_use_vars[this_yield] = var
                #check this
                self.add_constraint(DIMEUse, 
                                    hook = self,
                                    expr = expr - var,
                                    id_ = id_,
                                    species = species.id,
                                    lb = 0,
                                    ub = 0
                                    )
                
            # add a constraint to make sure only one yield is active
            expr = symbol_sum([var for var in yield_use_vars.values()])
            self.add_constraint(YieldChoose, 
                                hook = self,
                                expr = expr,
                                id_ = species.id,
                                species = species.id,
                                lb = 1,
                                ub = 1
                                )
    
            for rxn, dime_list in species.dime_coupling_dict.items(): # for this species, add the Inactivation constraint
                try:
                    expr = symbol_sum([dime.variable for dime in dime_list[UPTAKES]]) - rxn.upt_variable
                    self.add_constraint(UptakeInactivation,
                                    hook = rxn,
                                    expr = expr,
                                    lb = 0)
                except KeyError:
                    pass
            for rxn, dime_list in species.dime_coupling_dict.items(): # for this species, add the Inactivation constraint
                try:
                    expr = symbol_sum([dime.variable for dime in dime_list[SECRETIONS]]) - rxn.sec_variable
                    self.add_constraint(SecretionInactivation,
                                    hook = rxn,
                                    expr = expr,
                                    lb = 0)
                except KeyError:
                    pass
                
            # a metabolite cannot be uptaken and secreted at the same time by a species
            for rxn in species.rxns:
                try: # make sure that there is both uptake and secretion for this reaction
                    _ = species.dime_coupling_dict[rxn][UPTAKES] + \
                        species.dime_coupling_dict[rxn][SECRETIONS]
                    expr = rxn.upt_variable + rxn.sec_variable
                    self.add_constraint(SimultUse, 
                                        hook =rxn, 
                                        expr = expr,
                                        ub = 1)
                except KeyError: # if not, don't do anything
                    pass
                
            self.repair()
                
        for rxn in tqdm(self.reactions, desc='Adding the metabolites'):
            if 'met_type' in kwargs:
                met_id = get_met_id(rxn.id, rxn.species)
                try:
                    met_type = kwargs['met_type'][met_id]
                except KeyError:
                    met_type = 'biotic' # if the metabolite is not in the original environment
                self._add_mets_from_rxn(rxn, met_type)
            else:
                self._add_mets_from_rxn(rxn)
            
        # now that the metabolites are defined, we can define substrates and products
        for species in self.species:
            species.update_subs_prod()
        
        for met in tqdm(self.metabolites, desc='Adding the interactions'):
            self._add_interaction_from_met(met, **kwargs)
            
        self.repair()
        pos_int = self.get_variables_of_type(PositiveInteraction)
        upt_var = self.get_variables_of_type(UptakeActivation)
        for met in tqdm(self.metabolites, desc='Checking metabolites with biotic source'):
            # if the metabolite is biotic, a new constraint is needed to make sure it is produced
            if met.type == 'biotic': # if it is biotic it should be produced by a species and PI_ variable should exist
                try:
                    this_pi = pos_int.get_by_id(met.id)
                    for rxn in met.rxns:
                        this_ua = upt_var.get_by_id(rxn.id)
                        self.add_constraint(BioticProduction, 
                                                hook = rxn, 
                                                expr = this_ua - this_pi,
                                                ub = 0)
                except KeyError: # it is possible that this organism only produces this metabolite, so no need to do anything
                    pass
            
    def _add_interaction_from_met(self, met, sense='any', queue=False, **kwargs):
        """
        Will add the interaction related constraints and variables for each metabolite.

        Parameters
        ----------
        met : A core.CommunityMetabolite object
        sense : str, optional
            The direction of optimization. The acceptable values are 'min', 'max',
            and 'any'. The default is 'any' (undetermined direction).
        diff_interaction : Bool, optional
            To differentiate between positive and negative interactions? The default is True.

        Returns
        -------
        None.

        """
        bigM = 50
        
        producing = [rxn.sec_variable for rxn in met.rxns if met.id in \
                     self.species.get_by_id(rxn.species).products] # the producing reactions
            
        consuming = [rxn.upt_variable for rxn in met.rxns if met.id in \
                     self.species.get_by_id(rxn.species).substrates] # the consuming reactions

        consuming_species = [rxn.species for rxn in met.rxns if met.id in \
                     self.species.get_by_id(rxn.species).substrates]  # the consuming reactions

        producing_species = [rxn.species for rxn in met.rxns if met.id in \
                     self.species.get_by_id(rxn.species).products]  # the producing reactions

        union_species = consuming_species + producing_species

        
        if len(consuming) == 0: # if it cannot be consumed by any, no interaction exists and hence, no constarint is needed
            return
        
        # this should be added anyway
        #this is to ensure that if at least two species are consuming it add a neg int if only one doesnt matter!!
        if len(consuming) >1 :
            xn_j = self.add_variable(NegativeInteraction,
                                    met,
                                    lb = 0,
                                    ub = 1,
                                    queue=queue)

        #at least someone produces and at least someone consumes this needs to hold anyways
        if (len(producing) > 0) and (len(set(union_species))>1) :
            xp_j = self.add_variable(PositiveInteraction,
                                    met,
                                    lb = 0,
                                    ub = 1,
                                    queue=queue)
            
        
        if sense == 'any' or sense == 'min': 
            # these constarints force xi=1, so they should be applied if the objective is to minimize the interactions
            # todo normally we can differentiate between cooperation and competition
            #  todo but in case someone wants not to they can add less variables and constraints
            # todo for now we dont see a point in keeping this option but fyi
            if len(producing)>0: # cooperation can exist only if at least one can produce it
                # cooperation
                if (len(set(union_species)) > 1):
                    a_j = self.add_variable(SubstrateUse,
                                    met,
                                    lb = 0,
                                    ub = 1,
                                    queue=queue)
                    b_j = self.add_variable(ProductUse,
                                    met,
                                    lb = 0,
                                    ub = 1,
                                    queue=queue)


                    self.add_constraint(ActivePosInteraction,
                                        met,
                                        expr = xp_j+a_j+b_j,
                                        lb = 1)
                    self.add_constraint(ActiveSubsUse,
                                        met,
                                        expr = symbol_sum(consuming) + bigM*a_j,
                                        ub = bigM)
                    self.add_constraint(ActiveProdUse,
                                        met,
                                        expr = symbol_sum(producing) + bigM*b_j,
                                        ub = bigM)
                # competition
                if len(consuming) > 1:
                    self.add_constraint(ActiveNegInteraction,
                                   hook = met, 
                                   expr = -bigM*xn_j + symbol_sum(consuming),
                                   ub = 1)
            else:
                # both can be taken together
                all_rxns = list(set(consuming + producing)) # to make sure the list is unique
                if len(consuming) > 1:
                    self.add_constraint(ActiveNegInteraction,
                                   hook = met, 
                                   expr = -bigM*xn_j + symbol_sum(all_rxns),
                                   ub = 1)

        if sense == 'any' or sense == 'max': 
            # these constarints force xi=0, so they should be applied if the objective is to maximize the interactions
            if len(producing)>0: # cooperation can exist only if at least one can produce it
                # cooperation
                if (len(set(union_species)) > 1):
                    self.add_constraint(InactiveSubsUse,
                                        met,
                                        expr = xp_j - symbol_sum(consuming),
                                        ub = 0)
                    self.add_constraint(InactiveProdUse,
                                        met,
                                        expr = xp_j - symbol_sum(producing),
                                        ub = 0)
                #competition
                if len(consuming) > 1:
                    self.add_constraint(InactiveNegInteraction,
                                       hook = met,
                                       expr = bigM*xn_j - symbol_sum(consuming) + 2,
                                       ub = bigM)
            else:
                # both can be taken together
                all_rxns = list(set(consuming + producing)) # to make sure the list is unique

                if len(consuming) > 1:
                    self.add_constraint(InactiveNegInteraction,
                                       hook = met,
                                       expr = bigM*xn_j - symbol_sum(all_rxns) + 2,
                                       ub = bigM)

            
        self.repair()
        return
    
    def _add_mets_from_rxn(self, rxn, met_type='abiotic'):
        """
        Will add a metabolite to the model object.

        Parameters
        ----------
        rxn : A core.CommunityReaction object
        
        Returns
        -------
        None.

        """
        
        # met = rxn.metabolite # the metabolite object is not defined yet
        met_id = get_met_id(rxn.id, rxn.species)
        
        if met_id in self.metabolites:
            self.logger.warning('The metabolite {} was already in the model and hence,'
                 'neglected.'.format(met_id))
            met = self.metabolites.get_by_id(met_id)
        else:
            met = CommunityMetabolite(id = met_id, species=set(), rxns=set(), type_=met_type)
            self.metabolites += [met]
            # add the corresponding binary variable
            met._model = self
            
            
        met.species.add(rxn.species)
        met.rxns.add(rxn)
        rxn.metabolite = met
        
        self.repair()

    def remove_species(self,species_to_be_removed_list):
        """species_list = a list containing the id of the species to be removed as string"""
        sp_cons = get_all_subclasses(SpeciesConstraint)
        r_cons = get_all_subclasses(CommunityReactionConstraint)
        all_cons = sp_cons + r_cons
        sp_var = get_all_subclasses(SpeciesVariable)
        r_var = get_all_subclasses(CommunityReactionVariable)
        all_vars = sp_var + r_var
        for species_to_be_removed in species_to_be_removed_list:
            to_be_removed_cons = []
            to_be_removed_var = []
            #first store the constraints and variables to be removed for the species
            for cons in all_cons:
                try:
                    to_be_removed_cons += [k for k in self._cons_kinds[cons.__name__] if
                                           k.species == species_to_be_removed]
                except KeyError as e:
                    pass


            for var in all_vars:
                try:
                    to_be_removed_var += [k for k in self._var_kinds[var.__name__] if
                                          k.species == species_to_be_removed]
                except KeyError as e:
                    pass

            #remove them
            for cons in to_be_removed_cons:
                try:
                    self.remove_constraint(cons)
                except KeyError as e:
                    pass


            for var in to_be_removed_var:
                try:
                    self.remove_variable(var)
                except KeyError as e:
                    pass

            self.species.remove(self.species.get_by_id(species_to_be_removed))

        self.repair()

    def _add_rxns_from_dime(self, dime):
        """Add reactions to the model.

        Reactions with identifiers identical to a reaction already in the
        model are ignored.


        """
        this_species = self.species.get_by_id(dime.species)
        reaction_list = dime.rxns
        
        zm = dime.variable
        
        upt_act = self.get_variables_of_type(UptakeActivation)
        sec_act = self.get_variables_of_type(SecretionActivation)
        
        for rxn in reaction_list:
            if (rxn.id.split('EX')[1] in dime.substrates):
                try:
                    yj = upt_act.get_by_id(rxn.id)
                except KeyError:
                    if rxn.id not in self.reactions:
                        self.reactions += [rxn]
                        rxn._model = self
                    rxn.init_upt_variable()
                    yj = rxn.upt_variable
                try:
                    this_species.dime_coupling_dict[rxn][UPTAKES] += [dime]
                except KeyError:
                    try: # the first time as secretion
                        this_species.dime_coupling_dict[rxn][UPTAKES] = [dime]
                    except KeyError: # really the first time
                        this_species.dime_coupling_dict[rxn] = {}
                        this_species.dime_coupling_dict[rxn][UPTAKES] = [dime]
            
            elif (rxn.id.split('EX')[1] in dime.products):
                try:
                    yj = sec_act.get_by_id(rxn.id)
                except KeyError:
                    if rxn.id not in self.reactions:
                        self.reactions += [rxn]
                        rxn._model = self
                    rxn.init_sec_variable()
                    yj = rxn.sec_variable
                try:
                    this_species.dime_coupling_dict[rxn][SECRETIONS] += [dime]
                except KeyError:
                    try: # the first time as secretion
                        this_species.dime_coupling_dict[rxn][SECRETIONS] = [dime]
                    except KeyError: # really the first time
                        this_species.dime_coupling_dict[rxn] = {}
                        this_species.dime_coupling_dict[rxn][SECRETIONS] = [dime]
                
            self.add_constraint(ReactionActivation, 
                                hook = rxn, 
                                expr = yj - zm,
                                dime = dime,
                                lb = 0
                                )

        self.repair()
        
    def _add_simul_use_cstr(self, met):
        
        ''' make sure one metabolite cannot be uptaken and secreted by the same organism.'''
        
    
    
    def __deepcopy__(self, memo):
        """

        :param memo:
        :return:
        """

        return self.copy()    
    
    def copy(self):
        """Provides a parttial deepcopy of the model."""

        from ..io.dict import model_from_dict, model_to_dict
        from pytfa.optim.utils import copy_solver_configuration
        
        dictmodel = model_to_dict(self)
        new = model_from_dict(dictmodel)

        copy_solver_configuration(self, new)

        return new