#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pytfa.optim.constraints import ModelConstraint, MetaboliteConstraint, ReactionConstraint
    

class SpeciesConstraint(ModelConstraint):
    """
    Class to represent a constraint attached to the model, per species
    """

    def __init__(self, model, expr, id_, species, **kwargs):
        ModelConstraint.__init__(self,
                                   id_= id_,
                                   expr=expr,
                                   model=model,
                                   **kwargs)
        self._species = species
        
    @property
    def species(self):
        return self._species

    prefix = 'SC_'


class CommunityReactionConstraint(ReactionConstraint):
    '''A class associated to exchange reactions in the community model'''
    
    def __init__(self, reaction, expr, **kwargs):
        species = reaction.species

        ReactionConstraint.__init__(self,
                                   expr=expr,
                                   reaction=reaction,
                                   **kwargs)
        self._species = species
        
    @property
    def species(self):
        return self._species

    prefix = 'CRC_'
    
class ReactionActivation(CommunityReactionConstraint):
    '''A class to enforce a reaction to be active if one of its DIMES is active: 
        zm,k <= yj,k'''
    
    def __init__(self, reaction, expr, dime, **kwargs):
        self._dime = dime
        CommunityReactionConstraint.__init__(self, reaction, expr, **kwargs)
        
    
    @property
    def id(self):
        return self.reaction.id + '_' + self.dime.id
    
    @property
    def dime(self):
        return self._dime
        
    prefix = 'RAc_'
    
class UptakeInactivation(CommunityReactionConstraint):
    '''A class to enforce a reaction to be inactive if all of its DIMES are inactive: 
        sum(zm,k) >= yj,k'''
        
    prefix = 'UIAc_'
    
class SecretionInactivation(CommunityReactionConstraint):
    '''A class to enforce a reaction to be inactive if all of its DIMES are inactive: 
        sum(zm,k) >= yj,k'''
        
    prefix = 'SIAc_'
    
class DIMEUse(SpeciesConstraint):
    '''A class to enforce only one DIME to be active:
        sum(zm,k) = 1'''
        
    prefix = 'MU_'
    
class ActiveNegInteraction(MetaboliteConstraint):
    '''A class to check if an interaction exists: 
        sum(yj,k) -1 <= M*xi'''
        
    prefix = 'ANI_'
    
class InactiveNegInteraction(MetaboliteConstraint):
    '''A class to check if an interaction exists: 
        -sum(yj,k) + 2 <= M*(1-xi)'''
        
    prefix = 'INI_'
    
class ActivePosInteraction(MetaboliteConstraint):
    '''A class to check if an interaction exists: 
        1 - xp,i <= ai + bi'''
        
    prefix = 'API_'
    
class ActiveSubsUse(MetaboliteConstraint):
    '''A class to check if an interaction exists: 
        sum(yj,k) <= M(1 - ai)'''
        
    prefix = 'ASU_'
    
class ActiveProdUse(MetaboliteConstraint):
    '''A class to check if an interaction exists: 
        sum(yj,k) <= M(1 - bi)'''
        
    prefix = 'APU_'
    
class InactiveSubsUse(MetaboliteConstraint):
    '''A class to check if an interaction exists: 
        sum(yj,k) >= xp,i for substrates'''
        
    prefix = 'ISU_'
    
class InactiveProdUse(MetaboliteConstraint):
    '''A class to check if an interaction exists: 
        sum(yj,k) >= xp,i for products'''
        
    prefix = 'IPU_'
    
class SimultUse(CommunityReactionConstraint):
    '''A class of constraints to make sure a compound is not product 
    and substrate at the same time'''
    
    prefix = 'SUse_'
    
class YieldChoose(SpeciesConstraint):
    '''a class of constraints to make sure only one yield is active'''
    
    prefix = 'YC_'
    
class AbioticResourceConstraint(CommunityReactionConstraint):
     '''a class of constraints to check if a new metabolite is added to the medium'''
    
     prefix = 'AB_'


class AbioticResourceConstraint2(CommunityReactionConstraint):
    '''a class of constraints to check if a new metabolite is added to the medium'''

    prefix = 'AB2_'

class AbioticResourceConstraint3(CommunityReactionConstraint):
    '''a class of constraints to check if a new metabolite is added to the medium'''

    prefix = 'AB3_'
    
class BioticProduction(CommunityReactionConstraint):
    '''A class to ensure that biotic metabolites are produced by another organism'''
        
    prefix = 'BP_'
    



"THREE REACTION CONSTRAINTS TO LINEARIZE THE BINARY MULTIP"
class LinearizeInteger_1(ReactionConstraint):
    """
    Class to represent a constraint the implements
            z=x*y
    z<=x
    """

    prefix = 'LI1_'

class LinearizeInteger_2(ReactionConstraint):
    """
    Class to represent a constraint the implements
        z=x*y
    z<=y
    """

    prefix = 'LI2_'


class LinearizeInteger_3(ReactionConstraint):
    """
    Class to represent a constraint the implements
    z=x*y
    z>=x+y-1
    """

    prefix = 'LI3_'


class ConstraintMFA(ReactionConstraint):
    """
    Class to represent a constraint the implements
    MFA<=BU
    """

    prefix = 'C_MFA_'


class ConstraintMaxFlux(ReactionConstraint):
    """
    Class to represent a constraint the implements
    real constraint
    """

    prefix = 'C_MAXF'

class ConstraintSum(ReactionConstraint):
    """
    Class to represent a constraint the implements
    real constraint
    """

    prefix = 'C_MAXFSUM'


class ConstraintBond(ReactionConstraint):
    """
    Class to represent a constraint the implements
    real constraint
    """

    prefix = 'C_BOND'


class ConstraintEssential(ReactionConstraint):
    """
    Class to represent a constraint the implements
    real constraint
    """

    prefix = 'C_ESS'
    
class ExtraExchange(ModelConstraint):
    ''' A set of constraints to control the balance for the exchange of metabolites 
     with the environment'''
    prefix = 'EE_'
    
class IntraExchange(ModelConstraint):
    ''' A set of constraints to control the balance for the exchange of metabolites 
    between different species'''
    prefix = 'IE_'