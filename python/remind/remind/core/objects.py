#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from ..optim.variables import DIMEVariable, UptakeActivation, SecretionActivation,\
    InteractionVariable

from cobra import Metabolite, Reaction, Species, DictList
from pytfa.thermo.tmodel import ThermoModel
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)

# delattr(Reaction, 'forward_variable') # we do not need to keep these properties from cobra
# delattr(Reaction, 'reverse_variable') # we do not need to keep these properties from cobra

def throw_nomodel_error(id_):
    logger.warning('''{} has no model attached - variable attribute
     is not available'''.format(id_))
     
def sanitize_id(id_):
    
    if isinstance(id_, float):
        id_ = '_' + str(int(id_))
    return id_

class CommunityMetabolite(Metabolite):
    
    def __init__(self, species=set(), rxns=set(), type_='abiotic', **kwargs):
        Metabolite.__init__(self, **kwargs)
        self.species = species # a set to represent which species can consume or produce this compound
        self.rxns = rxns
        self.type = type_
        
    @property
    def constraint(self):
        # TODO: returns a specific constraint?
        raise NotImplementedError()
        

class CommunityReaction(Reaction):
    
    def __init__(self, metabolite=None, species_id=None, dimes=DictList(), **kwargs):
        
        Reaction.__init__(self, **kwargs)
        self.species = species_id # a string to represent the species that reaction belongs to
        self.dimes = dimes # DictList
        self._metabolite = metabolite # a single object
        
    @property
    def metabolite(self):
        return self._metabolite
    
    @metabolite.setter
    def metabolite(self, value):
        if isinstance(value, CommunityMetabolite):
            self._metabolite = value
        else:
            raise ValueError('The metabolite should be a CommunityMetabolite object.')
            
    @property
    def dimes(self):
        return self._dimes
    
    @dimes.setter
    def dimes(self, value):
        self._dimes = value
            
    def init_upt_variable(self, queue=False):
        """
        Attach an ExchangeActivation object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._internal_upt_variable = self.model.add_variable(UptakeActivation,
                                    self,
                                    lb = 0,
                                    ub = 1,
                                    queue=queue)
    @property
    def upt_variable(self):
        """
        For convenience in the equations of the constraints

        :return:
        """
        try:
            return self._internal_upt_variable.variable
        except AttributeError:
            throw_nomodel_error(self.id)
            
    def init_sec_variable(self, queue=False):
        """
        Attach an ExchangeActivation object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._internal_sec_variable = self.model.add_variable(SecretionActivation,
                                    self,
                                    lb = 0,
                                    ub = 1,
                                    queue=queue)
    @property
    def sec_variable(self):
        """
        For convenience in the equations of the constraints

        :return:
        """
        try:
            return self._internal_sec_variable.variable
        except AttributeError:
            throw_nomodel_error(self.id)

try:
    from pytfa.thermo.tmodel_struct import ThermoModelStructure
    inherit_class = ThermoModelStructure
except ImportError:
    inherit_class = ThermoModel

class CommunitySpecies(inherit_class):
    
    def __init__(self, model,
                 inplace=True,
                 *args, **kwargs):
        
        if not inplace:
            new = model.copy()
            self.__dict__ = new.__dict__
        else:
            # new = me_model
            self.__dict__ = model.__dict__
            
        self._dimes = DictList()
        self._rxns = DictList() # the boundary reactions
        self._substrates = DictList()
        self._products = DictList()

        self._yields = list() # the list of all yields in the biomass
        
        self.dime_coupling_dict = dict()
        
        self._model = None
        
    @property
    def dimes(self):
        return self._dimes
    
    @dimes.setter
    def dimes(self, value):
        self._dimes = value
        
    @property
    def rxns(self):
        return self._rxns
    
    @rxns.setter
    def rxns(self, value):
        self._rxns = value
        
    @property
    def substrates(self):
        return self._substrates
    
    @substrates.setter
    def substrates(self, value):
        self._substrates = value
        
    @property
    def products(self):
        return self._products
    
    @products.setter
    def products(self, value):
        self._products = value

    # @property
    # def elements(self):
    #     return self._elements
    
    # @elements.setter
    # def elements(self, value):
    #     self._elements = value

    @property
    def yields(self):
        return self._yields
    
    @yields.setter
    def yields(self, value):
        self._yields = value
        
    def add_dimes(self, coupling_dict):
        """
        This function initiates DIME objects, associates them with the species.
        Also, it initiates CommunityReaction objects and associates them to DIME and species.

        Parameters
        ----------
        coupling_dict : dict
            A dict where the key is DIME ID and the value is a list of reaction IDs

        Returns
        -------
        None.

        """
        
        for key, info in tqdm(coupling_dict.items(), desc='Adding the dime'):
            id_ = '{}_{}'.format(self.id , sanitize_id(key))
            if id_ in self.dimes:
                continue # not to add repeated objects

            the_yield = info['yield']
            this_dime = DIME(id = id_, species_id = self.id,
                             yield_=the_yield, rxn_set=set())
            self.dimes += [this_dime]
            

            if the_yield not in self.yields:
                self.yields += [the_yield]
                
            this_dime.substrates = [s.replace('EX','') for s in info['substrates']] # for now it's IDs
            this_dime.products = [s.replace('EX','') for s in info['products']] # for now it's IDs
            # add the reactions of this DIME
            rxn_list = info['reactions']
            for rxn_id in rxn_list:
                id_ = '{}_{}'.format(self.id,rxn_id) # the final ID is species_id + rxn_id
                try:
                    reaction = self.rxns.get_by_id(id_) # check if the reaction already exists
                except KeyError:
                    reaction = CommunityReaction(id = id_, species_id = self.id, dimes = DictList())
                    self.rxns += [reaction]

                # We should add metabolite to the reaction later, as it's defined per the community
                # add DIME to reaction
                reaction.dimes += [this_dime]
                # add reaction to dime
                this_dime.add_reactions(reaction)
                    
        return
    
    def update_subs_prod(self):
        '''
        The dime are already defined, the species are added to the model. And
        the metabolites are defined.
        '''
        
        try:
            model = self._model
        except AttributeError:
            throw_nomodel_error(self.id)
            
        for dime in tqdm(self.dimes, desc='Updating the list of substrates and products'):
            # Substrates
            new_subs = []
            for met_id in dime.substrates:
                met = model.metabolites.get_by_id(met_id) # replace ids by objects
                new_subs += [met]
                if met_id not in self.substrates:
                    self.substrates += [met]
            dime.substrates = new_subs
            # Products
            new_prods = []
            for met_id in dime.products:
                met = model.metabolites.get_by_id(met_id) # replace ids by objects
                new_prods += [met]
                if met_id not in self.products:
                    self.products += [met]
            dime.products = new_prods
        return         


class DIME(Species):

    def __init__(self, id=None, rxn_set=set(), species_id=None,
                 yield_='', *args, **kwargs):
        Species.__init__(self, id = id, *args, **kwargs)

        self.species = species_id
        self.rxns = rxn_set

        self._yield = yield_
        
        self._substrates = DictList()
        self._products = DictList()
        
    def add_reactions(self, rxns):
        """
        Adds reactions to an DIME.

        Parameters
        ----------
        rxns : List of CommunityReaction objects or a CommunityReaction

        Returns
        -------
        None.

        """
        e = 0 # is error?
        if hasattr(rxns, '__iter__'):
            if not isinstance(rxns[0], CommunityReaction):
                e = 1
        elif isinstance(rxns, CommunityReaction):
            rxns = [rxns]
        else:
            e = 1
            
        if e:
            raise ValueError('The input argument must be either a list or an'
                             'object of class CommunityReaction')
            
        # TODO: if the cstr or var exists, it is not that easy
        for x in rxns:
            self.rxns.add(x)
            
    @property
    def substrates(self):
        return self._substrates
    
    @substrates.setter
    def substrates(self, value):
        self._substrates = value
        
    @property
    def products(self):
        return self._products
    
    @products.setter
    def products(self, value):
        self._products = value
    
    def init_variable(self, queue=False):
        """
        Attach an DIMEVariable object to the Species. Needs to have the object
        attached to a model

        :return:
        """
        self._internal_variable = self.model.add_variable(DIMEVariable,
                                    hook = self.model,
                                    species = self.species,
                                    id_ = self.id,
                                    lb = 0,
                                    ub = 1,
                                    queue=queue)
    @property
    def variable(self):
        """
        For convenience in the equations of the constraints

        :return:
        """
        try:
            return self._internal_variable.variable
        except AttributeError:
            throw_nomodel_error(self.id)
