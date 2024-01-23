#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pytfa.optim.variables import ModelVariable, BinaryVariable, \
    MetaboliteVariable, ReactionVariable, get_binary_type
    

class SpeciesVariable(ModelVariable):
    """
    Class to represent a variables attached to the model, per species
    """

    def __init__(self, model, id_, species, **kwargs):
        ModelVariable.__init__(self,
                                   id_= id_,
                                   model=model,
                                   **kwargs)
        self._species = species
        
    @property
    def species(self):
        return self._species

    prefix = 'SV_'


class CommunityReactionVariable(ReactionVariable):
    '''A class associated to exchange reactions in the community model'''
    
    def __init__(self, reaction, **kwargs):
        species = reaction.species

        ReactionVariable.__init__(self,
                                   reaction,
                                   **kwargs)
        self._species = species
    
    @property
    def species(self):
        return self._species

    prefix = 'CRV_'
    
    
class DIMEVariable(SpeciesVariable, BinaryVariable):
    'A class representing the DIME associated binary variables (zm,k).'
    
    def __init__(self, model, id_, species, **kwargs):
        SpeciesVariable.__init__(self, 
                                model = model,
                                id_ = id_,
                                species = species,
                                type=get_binary_type(),
                                **kwargs)

    prefix = 'MXV_'
    
class ExchangeActivation(CommunityReactionVariable, BinaryVariable):
    'A class representing the exchange reaction binary variables (yj,k).'
    
    def __init__(self, reaction, **kwargs):
        CommunityReactionVariable.__init__(self, 
                                           reaction = reaction,
                                           type=get_binary_type(),
                                           **kwargs)

    prefix = 'XA_'
    
class UptakeActivation(ExchangeActivation):
    
    prefix = 'UA_'
    
class SecretionActivation(ExchangeActivation):
    
    prefix = 'SA_'
    
class InteractionVariable(MetaboliteVariable, BinaryVariable):
    'A class representing the interaction binary variables (xi).'
    
    def __init__(self, metabolite, **kwargs):
        MetaboliteVariable.__init__(self, 
                                    metabolite = metabolite,
                                    type=get_binary_type(),
                                    **kwargs)

    prefix = 'InV_'

class PositiveInteraction(InteractionVariable):
    'A class of binary variables to account for cooperation (xp,i).'
    
    prefix = 'PI_'
    
class NegativeInteraction(InteractionVariable):
    'A class of binary variables to account for competition (xn,i).'
    
    prefix = 'NI_'

class SubstrateUse(InteractionVariable):
    'A class of binary variables to account for the use of the substrate (ai).'
    
    prefix = 'SbU_'
    
class ProductUse(InteractionVariable):
    'A class of binary variables to account for the use of the product (bi).'
    
    prefix = 'PrU_'

class YieldUse(SpeciesVariable, BinaryVariable):
    'A class representing the yield associated binary variables (zm,k).'
    
    def __init__(self, model, id_, species, **kwargs):
        SpeciesVariable.__init__(self, 
                                model = model,
                                id_ = id_,
                                species = species,
                                type=get_binary_type(),
                                **kwargs)

    prefix = 'YU_'
    
class AbioticResource(InteractionVariable):
    'A class of binary variables to account for new metabolites to be added to the environment for n+(-)1.'
    
    prefix = 'AV_'
    
    
class MaximalFluxActiveVariable(ReactionVariable, BinaryVariable):
    """
    Class to represent a type of binary variable used to tell whether the
    max_flux constraint is active reaction is active or not such that:
    similar to the basal fluxconstraint
        B_+max_flux*BU*MFA<=(1-MFA)M
        MFA<=BU --> if BU not active MFA not active by definition
        linearization of BU*MFA is required linearize two integers
    """

    def __init__(self, reaction, **kwargs):
        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = 1

        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

    prefix = 'MFA_'


class MaximalFluxActiveVariableLinearize(ReactionVariable, BinaryVariable):
    """
    Class to represent the binary  multiplication variable
    MFA*BU=MFAL

        MFA<=BU --> if BU not active MFA not active by definition
        linearization of BU*MFA is required linearize two integers
    """

    def __init__(self, reaction, **kwargs):
        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = 1

        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

    prefix = 'MFAL_'

class MaximalFluxActiveVariableSum(ReactionVariable, BinaryVariable):
    """
    Class to represent the binary  multiplication variable
    MFA*BU=MFAL
    represents a variable of type BU+MFA
    to postprocess solutions like a helper variable
    BU+MFA=0 backward flux is not active
    BU+MFA=1 required more than min flux so potential main C source

    """

    def __init__(self, reaction, **kwargs):
        ReactionVariable.__init__(self, reaction,
                                  type='integer',
                                  **kwargs)

    prefix = 'MFAS_'