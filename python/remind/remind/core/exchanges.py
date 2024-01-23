
import pandas as pd
import numpy as np
from optlang.exceptions import SolverError
from optlang.interface import INFEASIBLE, UNBOUNDED
from pytfa.optim.constraints import ReactionConstraint, ModelConstraint, LinearizationConstraint
from pytfa.optim.variables import ForwardBackwardUseVariable, ModelVariable, ReactionVariable
from pytfa.optim.utils import symbol_sum
import sympy
from cobra.util.solver import OptimizationError

class ActiveReaction(ReactionConstraint):
    """
    Class to represent a constraint the implements
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'AR_'

class ActiveReaction2(ReactionConstraint):
    """
    Class to represent a constraint the implements
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'AR2_'
    
class ActiveReaction3(ReactionConstraint):
    """
    Class to represent a constraint the implements
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'AR3_'
    
class FractionalConv(ModelVariable):
    'A scalar varaible to represent the denominator in the objective function'
    
    prefix = 'FrC_'
    
class PeterLinearizeForward(ReactionVariable):
    """
    """
    prefix = 'LZF_'
    
class PeterLinearizeReverse(ReactionVariable):
    """
    """
    prefix = 'LZR_'
    
class UpperBoundConv(ReactionConstraint):
    """
    """
    prefix = 'UBC_'
    
class LowerBoundConv(ReactionConstraint):
    """
    """
    prefix = 'LBC_'



"Some comments"
"if no thermo we should ensure that when FU or BU=1 flux is not Zero"
"When looking for iME's of higer sizes change size parameter"
#todo we should verify this I think if we go higher size we need to add additional constraints

def imeanalysis(tmodel, reactions=[], uptakes=[], secretions=[],
                element=None, bbb=None, getAlternatives=True,
                enforce_size=True, size=0, max_alternative=1000):
    """"
    tmodel:tfa model 
    reactions: reaction objects to consider
    enforce_size: Boolean, There are two possibilities:
       1) either no constrain and optimize --> enforce_size=False
       2) fix size and generate alternatives of this size --> enforce_size=True
    size: a parameter the extent of suboptimality. If size=0 --> find optimal solutions
    """
    #tmodel=_tmodel.copy()
    
    # make sure the reactions, uptakes and secretions are dijoint sets
    if len(set(reactions).intersection(uptakes)) or \
       len(set(reactions).intersection(secretions)) or \
       len(set(secretions).intersection(uptakes)):
           
        raise ValueError('The lists reactions, uptakes, and secretions must be mutually disjoint.')
    
    
    bfuse_variables = []
    'for ime minimize active uptakes by maximizing non active ones'
    for rxn in reactions: # constrain both uptakes and secretions
        f_use = tmodel.forward_use_variable.get_by_id(rxn.id)
        b_use = tmodel.backward_use_variable.get_by_id(rxn.id)


        #its a reaction variable needs to be rxn type the name
        bfuse = tmodel.add_variable(ForwardBackwardUseVariable,rxn)
        bfuse_variables += [bfuse]


        #expression based on use variables
        expression = f_use + b_use + bfuse

        cons0 = tmodel.add_constraint(ActiveReaction,
                                     rxn,
                                     expression,
                                     lb = 0,ub = 1)

        #todo why is this here
    for rxn in uptakes: # constrain only uptakes
        b_use = tmodel.backward_use_variable.get_by_id(rxn.id)


        #its a reaction variable needs to be rxn type the name
        bfuse = tmodel.add_variable(ForwardBackwardUseVariable,rxn)
        bfuse_variables += [bfuse]


        #expression based on use variables
        expression = b_use + bfuse

        cons0 = tmodel.add_constraint(ActiveReaction,
                                     rxn,
                                     expression,
                                     lb = 0,ub = 1)
        
    for rxn in secretions: # constrain only secretions
        f_use = tmodel.forward_use_variable.get_by_id(rxn.id)


        #its a reaction variable needs to be rxn type the name
        bfuse = tmodel.add_variable(ForwardBackwardUseVariable,rxn)
        bfuse_variables += [bfuse]


        #expression based on use variables
        expression = f_use + bfuse

        cons0 = tmodel.add_constraint(ActiveReaction,
                                     rxn,
                                     expression,
                                     lb = 0,ub = 1)


    tmodel.repair()
    
    # setting the objective function
    # BFUSE_VARS = tmodel._var_kinds[ForwardBackwardUseVariable.__name__]
    bfuse_var_names = [x.name for x in bfuse_variables]

    # tmodel.objective = symbol_sum(bfuse_variables) #- subs_penalty*symbol_sum(b_use_variables)
    # if the integer cut variables and non-integer cut variables contribute equally
    # to the objective, it is possible that the iMEs vary with the order of found solutions.
    # To avoid this different weights are used for integer cut variables and non-integer cut variables
    # This means we would prefer not to secret the same element as what are are uptaking.
    obj_expression = symbol_sum(bfuse_variables) 
    tmodel.objective = obj_expression
    tmodel.objective_direction = 'max'
    # todo define objective direction
    sol = tmodel.optimize()
    list_=[]
    df_all=pd.DataFrame()
    df_all['Alternative_1']=np.zeros(len(tmodel.variables))
    df_all.index=sol.raw.index
    # minimal_set = tmodel.objective.value
    if getAlternatives:
        #add integer cut constraints
        alternative=1
        while((tmodel.solver.status != INFEASIBLE) and (alternative<=max_alternative)):
            print('Generating alternative', alternative)
            #todo look for primal from lUMPGEM
            df_sol = sol.raw
            df_flux=sol.fluxes
            df_bfuse = df_sol.loc[bfuse_var_names]
            idx = df_bfuse[np.isclose(df_bfuse,0)].index

            active_vars=[tmodel._var_dict[k] for k in idx]
            df = pd.DataFrame(columns=['metabolites', 'flux', 'FU', 'BU'])
            
            active_rxn_id = [i.reaction.id for i in active_vars]
            fu_vars = [tmodel.forward_use_variable.get_by_id(k).name for k in active_rxn_id]
            bu_vars = [tmodel.backward_use_variable.get_by_id(k).name for k in active_rxn_id]
            df['metabolites'] = active_rxn_id
            df['BFUSE']=df_bfuse[np.isclose(df_bfuse,0)].values
            'error was here this var is only for forward take net flux'
            df['flux'] = df_flux[active_rxn_id].values
            df['FU'] = df_sol[fu_vars].values
            df['BU'] = df_sol[bu_vars].values
            df['alternative'] = np.ones(len(active_rxn_id)) * alternative
            if element is not None:
                df['element'] = [element]*len(active_rxn_id)
            if bbb is not None:
                df['BBB'] = [bbb]*len(active_rxn_id)
            df_all['Alternative_{}'.format(alternative)]=df_sol.values
            list_.append(df)

            'forbidden profile constraint from pytfa can also be used !'
            'integer cut constraint'
            expression2 = symbol_sum([var for var in active_vars \
                                      # if var in integer_cut_vars \
                                      ])
            cons1 = tmodel.add_constraint(ModelConstraint,
                                                 tmodel,
                                                 expression2,
                                                 id_='dime_constraint_integer_cut_{}'.format(alternative),
                                                #maybe later set a real upperbound
                                                 # I think we don't need an upper bound
                                                 lb = 1)
            
            if enforce_size:
                # The idea is to fix the number of active exchanges, and solve the problem for a dummy objective
                upper_bound = df_bfuse.sum()
                expression = symbol_sum(bfuse_variables)
                eps=1e-5
                'THIS CONSTRAINT WE CAN ADD OOUTSIDE OF THE LOOP'
                cons0 = tmodel.add_constraint(ModelConstraint,
                                                     tmodel,
                                                     expression,
                                                     id_='imm_constraint',
                                              #not to fix it there
                                                    lb=upper_bound-size-eps,
                                                    ub=upper_bound-size+eps)
                tmodel.objective= sympy.S.Zero
                
            try:
                sol=tmodel.optimize()
            except: # warning: this catches any error. In CPLEX, it's SolverError but in gurobi it's AttributeError
                print('Number of alternatives generated is ', alternative)
                break
                
            #check if its so slow or not
            tmodel.repair()

            alternative += 1

    return list_,df_all


def immanalysis(tmodel,reactions, eps=1e-5, getAlternatives=True, enforce_size=True,
                size=0, max_alternative=1000):
	
	
    consume = []
    product = []
    for rxn in reactions:
        f_use = tmodel.forward_use_variable.get_by_id(rxn.id)
        b_use = tmodel.backward_use_variable.get_by_id(rxn.id)


        f_flux = tmodel.reactions.get_by_id(rxn.id).forward_variable
        b_flux = tmodel.reactions.get_by_id(rxn.id).reverse_variable
        #its a reaction variable needs to be rxn type the name
        # No need to add new variables
	
        expression = b_use*eps - b_flux
        tmodel.add_constraint(ActiveReaction,
                             rxn,
                             expression,
                             ub = 0)
        
        consume += [b_use]
        product += [f_use] # also to get secretions for comparison with the rest
	
	        
    tmodel.repair()
    
    tmodel.objective = symbol_sum(consume)
    tmodel.objective_direction = 'min'
    # todo define objective direction
    sol = tmodel.optimize()
    list_=[]
    df_all=pd.DataFrame()
    df_all['Alternative_1']=np.zeros(len(tmodel.variables))
    df_all.index=sol.raw.index
    # minimal_set = tmodel.objective.value
    if getAlternatives:
        #add integer cut constraints
        alternative=1
        while((tmodel.solver.status != INFEASIBLE) and (alternative<=max_alternative)):
            print('Generating alternative', alternative)
            #todo look for primal from lUMPGEM
            df_sol = sol.raw
            df_flux=sol.fluxes
            df_buse = df_sol.loc[[var.name for var in consume]]
            idx_b = df_buse[np.isclose(df_buse,1)].index
            df_fuse = df_sol.loc[[var.name for var in product]]
            # df_f_flux = df_flux[[var.reaction.id for var in product]] # there is no guarntee if b_f is 1, the flux is nonzero, so we take the flux directly
            # idx_f = df_fuse[df_f_flux>0].index

            active_vars=[tmodel._var_dict[k] for k in idx_b] + \
                [tmodel._var_dict[k] for k in df_fuse.index if tmodel._var_dict[k].reaction.flux>0]
            df = pd.DataFrame(columns=['metabolites', 'flux', 'FU', 'BU'])
            active_rxn_id = [i.reaction.id for i in active_vars]
            fu_vars = [tmodel.forward_use_variable.get_by_id(k).name for k in active_rxn_id]
            bu_vars = [tmodel.backward_use_variable.get_by_id(k).name for k in active_rxn_id]
            df['metabolites'] = active_rxn_id
            'error was here this var is only for forward take net flux'
            df['flux'] = df_flux[active_rxn_id].values
            df['FU'] = df_sol[fu_vars].values
            df['BU'] = df_sol[bu_vars].values
            df['alternative'] = np.ones(len(active_rxn_id)) * alternative
            df_all['Alternative_{}'.format(alternative)]=df_sol.values
            list_.append(df)

            'forbidden profile constraint from pytfa can also be used !'
            'integer cut constraint'
            expression2 = symbol_sum(active_vars)
            tmodel.add_constraint(ModelConstraint,
                                                 tmodel,
                                                 expression2,
                                                 id_='dime_constraint_integer_cut_{}'.format(alternative),
                                                #maybe later set a real upperbound
                                                 # I think we don't need an upper bound
                                                 lb = 1)
            
            if enforce_size:
                # The idea is to fix the number of active exchanges, and solve the problem for a dummy objective
                upper_bound = tmodel.objective.value
                expression = tmodel.objective.expression
                tol=1e-5
                'THIS CONSTRAINT WE CAN ADD OOUTSIDE OF THE LOOP'
                tmodel.add_constraint(ModelConstraint,
                                                     tmodel,
                                                     expression,
                                                     id_='imm_constraint',
                                              #not to fix it there
                                                    lb=upper_bound-size-tol,
                                                    ub=upper_bound-size+tol)
                tmodel.objective= sympy.S.Zero
                
            try:
                sol=tmodel.optimize()
            except: # warning: this catches any error. In CPLEX, it's SolverError but in gurobi it's AttributeError
                print('Number of alternatives generated is ', alternative)
                break
                
            #check if its so slow or not
            tmodel.repair()

            alternative += 1

    return list_,df_all


# this version involves also the stop cndition for now stop_cond is built in and set for min+2 size
def imeanalysis_all_size(tmodel,reactions,getAlternatives=True,max_alternative=1000,biomass='Growth',stop_cond=True,criteria=2):
    "tmodel:tfa model reactions: reaction objects to consider"
    #tmodel=_tmodel.copy()
    'for imm minimize active uptakes by maximizing non active ones'
    for rxn in reactions:
        FU = tmodel.forward_use_variable.get_by_id(rxn.id)
        BU = tmodel.backward_use_variable.get_by_id(rxn.id)


        F_ = tmodel.reactions.get_by_id(rxn.id).forward_variable
        B_ = tmodel.reactions.get_by_id(rxn.id).reverse_variable
        #its a reaction variable needs to be rxn type the name
        BFUSE = tmodel.add_variable(ForwardBackwardUseVariable,rxn)


        #expression based on use variables
        expression = FU+BU+50*BFUSE
        #whynot FU+BU+BFUSE=1
        #expression2= BFUSE+B_

        'if bfuse=1 BU=0 why we have this C --> ? from redGEMx '
        'not required' #todo remove later
        cons0 = tmodel.add_constraint(ActiveReaction,
                                     rxn,
                                     expression,
                                     lb = 0,ub = 50)




        #to ensure that when BFUSE=0 R_rxn carries a flux R_rxn!=0
        # eps=1e-7
        # cons1 = tmodel.add_constraint(ActiveReaction2,
        #                              rxn,
        #                              expression2,
        #                              lb = eps,ub = 50)

    tmodel.repair()
    #mytfa_imm = imeAnalysis(mytfa, mytfa.exchanges)
    BFUSE_VARS = tmodel._var_kinds[ForwardBackwardUseVariable.__name__]
    BFUSE_var_names = [x.name for x in BFUSE_VARS]

    expression = sum(BFUSE_VARS)

    tmodel.objective = expression
    tmodel.objective_direction = 'max'
    # todo define objective direction
    sol = tmodel.optimize()
    list_=[]
    df_all=pd.DataFrame()
    df_all['Alternative_1']=np.zeros(len(tmodel.variables))
    df_all.index=sol.raw.index

    #to get min_Size
    df_sol = sol.raw
    df_flux = sol.fluxes
    df_bfuse = df_sol.loc[BFUSE_var_names]
    # todo check integer feasibility <1 but close to 1 or no !
    idx = df_bfuse[df_bfuse < 0.5].index
    active_vars = [tmodel._var_dict[k] for k in idx]
    active_rxn_id = [i.reaction.id for i in active_vars]
    min_size = len(active_rxn_id)
    if getAlternatives:
        #add integer cut constraints
        alternative=1
        while((tmodel.solver.status != INFEASIBLE) and (alternative<=max_alternative)):

            #todo look for primal from lUMPGEM
            df_sol = sol.raw
            df_flux=sol.fluxes
            df_bfuse = df_sol.loc[BFUSE_var_names]
            #todo check integer feasibility <1 but close to 1 or no !
            idx = df_bfuse[df_bfuse < 0.5].index
            active_vars=[tmodel._var_dict[k] for k in idx]
            df = pd.DataFrame(columns=['metabolites', 'flux', 'FU', 'BU'])
            active_rxn_id = [i.reaction.id for i in active_vars]

            if stop_cond:
                if (len(active_rxn_id)-min_size>criteria):
                    print('Stopping criteria is reached at alternative {} min+{} alternatives are exhausted for model'.format(alternative, criteria,tmodel.id))
                    break
            fu_vars = [tmodel.forward_use_variable.get_by_id(k).name for k in active_rxn_id]
            bu_vars = [tmodel.backward_use_variable.get_by_id(k).name for k in active_rxn_id]
            df['metabolites'] = active_rxn_id
            df['BFUSE']=df_bfuse[df_bfuse < 0.5].values
            'error was here this var is only for forward take net flux'
            df['flux'] = df_flux[active_rxn_id].values
            df['FU'] = df_sol[fu_vars].values
            df['BU'] = df_sol[bu_vars].values
            df['alternative'] = np.ones(len(active_rxn_id)) * alternative
            #add growth also as a column
            df['Growth']=np.ones(len(active_rxn_id)) * sol.fluxes[biomass]
            df['status']=len(active_rxn_id) * [sol.status]

            df_all['Alternative_{}'.format(alternative)]=df_sol.values

            #expr = [k.reverse_variable.name for k in tmodel.exchanges]
            #df['uptake_yield']=np.ones(len(active_rxn_id)) * df_all.loc[expr].sum()
            df['uptake_yield']=np.ones(len(active_rxn_id))*abs(df[df.BU>0.5].flux.sum())

            print('Generating alternative {} of size {} for model {}'.format(alternative,len(active_rxn_id),tmodel.id))
            list_.append(df)

            #todo fix NOW!
            #upper_bound = sol.objective_value
            # upper_bound = df_bfuse.sum()
            #
            expression = symbol_sum(BFUSE_VARS)
            # eps=1e-5
            # 'FIRST CONSTRAINT WE CAN ADD OOUTSIDE OF THE LOOP'
            # 'forbidden profile'
            # cons0 = tmodel.add_constraint(ModelConstraint,
            #                                      tmodel,
            #                                      expression,
            #                                      id_='imm_constraint',
            #                               #not to fix it there
            #                                     lb=upper_bound-size-eps,
            #                                     ub=upper_bound-size)
            #                                      #lb=upper_bound * 0.999999, ub=upper_bound)



            expression2 = sum(active_vars)
            cons1 = tmodel.add_constraint(ModelConstraint,
                                                 tmodel,
                                                 expression2,
                                                 id_='dime_constraint_integer_cut_{}'.format(alternative),
                                                #maybe later set a real upperbound
                                                 lb= 0.5, ub=len(active_rxn_id))

            #check if its so slow or not
            tmodel.repair()

            try:
               # tmodel.objective= sympy.S.Zero
                sol=tmodel.optimize()
                print("Solver status is {}".format(sol.status))
            except SolverError:
                print('Number of alternatives generated is ',alternative)
                break
            except OptimizationError:
                #todo in the case that thre s optim error the same solution will be saved twice
                 print("optimization error detected")
                 tmodel.objective_direction = 'min'
                 tmodel.objective_direction = 'max'
            alternative += 1

            #to store the values  in wgat type of frame work hdf5


    return list_,df_all



#for now stop_cond is built in and set for min+2 size
def imeanalysis_once_size(tmodel,reactions,max_alternative=1,biomass='Growth'):
    "tmodel:tfa model reactions: reaction objects to consider"
    #tmodel=_tmodel.copy()
    'for imm minimize active uptakes by maximizing non active ones'
    for rxn in reactions:
        FU = tmodel.forward_use_variable.get_by_id(rxn.id)
        BU = tmodel.backward_use_variable.get_by_id(rxn.id)


        F_ = tmodel.reactions.get_by_id(rxn.id).forward_variable
        B_ = tmodel.reactions.get_by_id(rxn.id).reverse_variable
        #its a reaction variable needs to be rxn type the name
        BFUSE = tmodel.add_variable(ForwardBackwardUseVariable,rxn)


        #expression based on use variables
        expression = FU+BU+50*BFUSE
        #whynot FU+BU+BFUSE=1
        #expression2= BFUSE+B_

        'if bfuse=1 BU=0 why we have this C --> ? from redGEMx '
        'not required' #todo remove later
        cons0 = tmodel.add_constraint(ActiveReaction,
                                     rxn,
                                     expression,
                                     lb = 0,ub = 50)




        #to ensure that when BFUSE=0 R_rxn carries a flux R_rxn!=0
        # eps=1e-7
        # cons1 = tmodel.add_constraint(ActiveReaction2,
        #                              rxn,
        #                              expression2,
        #                              lb = eps,ub = 50)

    tmodel.repair()
    #mytfa_imm = imeAnalysis(mytfa, mytfa.exchanges)
    BFUSE_VARS = tmodel._var_kinds[ForwardBackwardUseVariable.__name__]
    BFUSE_var_names = [x.name for x in BFUSE_VARS]

    expression = sum(BFUSE_VARS)

    tmodel.objective = expression
    tmodel.objective_direction = 'max'
    # todo define objective direction
    sol = tmodel.optimize()
    list_=[]
    df_all=pd.DataFrame()
    df_all['Alternative_1']=np.zeros(len(tmodel.variables))
    df_all.index=sol.raw.index

    #to get min_Size
    df_sol = sol.raw
    df_flux = sol.fluxes
    df_bfuse = df_sol.loc[BFUSE_var_names]
    # todo check integer feasibility <1 but close to 1 or no !
    idx = df_bfuse[df_bfuse < 0.5].index
    active_vars = [tmodel._var_dict[k] for k in idx]
    active_rxn_id = [i.reaction.id for i in active_vars]
    min_size = len(active_rxn_id)
    obj=sol.objective_value

    return min_size,obj



''' this was a trial in the end forcing does not work ! but as a worrkflow of if and if elses '''


def imeanalysis_all_size_w_min_mw(tmodel,reactions,getAlternatives=True,max_alternative=1000,biomass='Growth'):
    "tmodel:tfa model reactions: reaction objects to consider"
    'in addition to ime here consider satisfying the size and having exchanges of minimal mw'
    #tmodel=_tmodel.copy()
    'for imm minimize active uptakes by maximizing non active ones'

    "enforce basal flux there should be at least a flux if FU OR BU is positive"
    basal_flux=1e-5
    for rxn in reactions:
        FU = tmodel.forward_use_variable.get_by_id(rxn.id)
        BU = tmodel.backward_use_variable.get_by_id(rxn.id)


        F_ = tmodel.reactions.get_by_id(rxn.id).forward_variable
        B_ = tmodel.reactions.get_by_id(rxn.id).reverse_variable
        BFUSE = tmodel.add_variable(ForwardBackwardUseVariable,rxn)

        expression = FU+BU+BFUSE

        expressionBackward= B_-basal_flux*BU
        expressionForward= F_-basal_flux*FU


        'if bfuse=1 BU=0 why we have this C --> ? from redGEMx '
        'not required' #todo remove later
        # cons0 = tmodel.add_constraint(ActiveReaction,
        #                              rxn,
        #                              expression,
        #                              lb = 0,ub = 50)

        cons0 = tmodel.add_constraint(ActiveReaction,
                                     rxn,
                                     expression,
                                     lb = 1.0,ub = 1.0)

        #to ensure that when BFUSE=0 R_rxn carries a flux R_rxn!=0
        # eps=1e-7
        'to ensure that rxns carry flux if FU or BU is active'
        consBackward = tmodel.add_constraint(ActiveReactionBackward,
                                     rxn,
                                     expressionBackward,
                                     lb = 0)

        consForward = tmodel.add_constraint(ActiveReactionForward,
                                     rxn,
                                     expressionForward,
                                     lb = 0)

    tmodel.repair()
    BFUSE_VARS = tmodel._var_kinds[ForwardBackwardUseVariable.__name__]
    BFUSE_var_names = [x.name for x in BFUSE_VARS]

    expression = sum(BFUSE_VARS)

    tmodel.objective = expression
    tmodel.objective_direction = 'max'
    # todo define objective direction
    sol = tmodel.optimize()
    list_=[]
    df_all=pd.DataFrame()
    df_all['Alternative_1']=np.zeros(len(tmodel.variables))
    df_all.index=sol.raw.index
    addIntegerCuts=False
    if getAlternatives:
        #add integer cut constraints
        alternative=1
        while((tmodel.solver.status != INFEASIBLE) and (alternative<=max_alternative)):
            #todo do not store the first and add integer
            #todo look for primal from lUMPGEM
            df_sol = sol.raw
            df_bfuse = df_sol.loc[BFUSE_var_names]

            upper_bound = df_bfuse.sum()
            expression = symbol_sum(BFUSE_VARS)
            eps=1e-5
            'FIRST CONSTRAINT WE CAN ADD OOUTSIDE OF THE LOOP'
            'forbidden profile'
            size=0
            check_constraint_size=tmodel.constraints.get('MODC_ime_size_constraint')
            if check_constraint_size is None:
                cons1 = tmodel.add_constraint(ModelConstraint,
                                                     tmodel,
                                                     expression,
                                                     id_='ime_size_constraint',
                                              #not to fix it there
                                                    #lb=0.0,
                                                     lb=upper_bound - size - eps,
                                                    ub=upper_bound-size)


                tmodel.repair()

            else:
                tmodel.constraints.get('MODC_ime_size_constraint').lb=upper_bound-eps
                #tmodel.constraints.get('MODC_ime_size_constraint').lb=0.0

                tmodel.constraints.get('MODC_ime_size_constraint').ub=upper_bound
                tmodel.repair()

            check_constraint_mw=tmodel.constraints.get('MODC_mw_cons_size')
            if check_constraint_mw is None:
            #CALL CARBON SOURCES
                all_sources=[rxn.id for rxn in reactions]
                c_sources = listcarbonsources(tmodel)

                expr_mw=[(1.0-tmodel.variables.get('BFUSE_'+k))*(tmodel.metabolites.get_by_id(k.replace('EX_','')).formula_weight) for k in c_sources]
                expression_mw=symbol_sum(expr_mw)

                tmodel.objective=expression_mw
                tmodel.objective_direction='min'
                sol_min_mw=tmodel.optimize()
                tmodel.objective_direction='max'
                sol_max_mw=tmodel.optimize()
                print('minimum mw uptake is {} maximum uptake is {}'.format(sol_min_mw.objective_value,sol_max_mw.objective_value))

                tol_mw=1.5

                cons_min_mw=tmodel.add_constraint(ModelConstraint,
                                                 tmodel,
                                                 expression_mw,
                                                 id_='mw_cons_size',#.format(len(active_vars)),
                                                 lb=(sol_min_mw.objective_value),
                                                  ub=(sol_min_mw.objective_value)*tol_mw)

                print('Constraint for uptake is minimum  is {} maximum is {}'.format(sol_min_mw.objective_value,
                                                                        sol_min_mw.objective_value*tol_mw))

                #check variability
                from pytfa.analysis.variability import variability_analysis, _variability_analysis_element
                var = _variability_analysis_element(tmodel, tmodel.reactions.EX_cit_e, sense='min')
                print('min for citrate is {}'.format(var))


            "for integer cuts"
            if addIntegerCuts:
                expression2 = sum(active_vars)
                cons2 = tmodel.add_constraint(ModelConstraint,
                                                     tmodel,
                                                     expression2,
                                                     id_='dime_constraint_integer_cut_{}'.format(alternative-1),
                                                    #maybe later set a real upperbound
                                                     lb= 0.5, ub=len(active_rxn_id))

            tmodel.repair()
            tmodel.objective=expression
            tmodel.objective_direction='max'
            try:
                'here store values'
                sol=tmodel.optimize()
                df_sol = sol.raw
                df_flux = sol.fluxes
                df_bfuse = df_sol.loc[BFUSE_var_names]
                idx = df_bfuse[df_bfuse < 0.5].index
                active_vars = [tmodel._var_dict[k] for k in idx]
                df = pd.DataFrame(columns=['metabolites', 'flux', 'FU', 'BU'])
                active_rxn_id = [i.reaction.id for i in active_vars]
                fu_vars = [tmodel.forward_use_variable.get_by_id(k).name for k in active_rxn_id]
                bu_vars = [tmodel.backward_use_variable.get_by_id(k).name for k in active_rxn_id]
                df['metabolites'] = active_rxn_id
                df['BFUSE'] = df_bfuse[df_bfuse < 0.5].values
                'error was here this var is only for forward take net flux'
                df['flux'] = df_flux[active_rxn_id].values
                df['FU'] = df_sol[fu_vars].values
                df['BU'] = df_sol[bu_vars].values
                df['alternative'] = np.ones(len(active_rxn_id)) * alternative
                # add growth also as a column
                df['Growth'] = np.ones(len(active_rxn_id)) * sol.fluxes[biomass]
                df_all['Alternative_{}'.format(alternative)] = df_sol.values
                print('Generating alternative {} of size {}'.format(alternative, len(active_rxn_id)))
                list_.append(df)
                addIntegerCuts = True
                alternative += 1

            except SolverError:
                #first remove constraints
                #tmodel.remove_constraint(tmodel.constraints.MODC_ime_size_constraint)
                tmodel.remove_constraint(tmodel.constraints.MODC_mw_cons_size)
                tmodel.constraints.get('MODC_ime_size_constraint').lb=0.0
                tmodel.constraints.get('MODC_ime_size_constraint').ub=upper_bound-1.0
                tmodel.repair()
                #tmodel.constraints.
                #try to solve again you have ime cuts but u dont have mw constraint or the size
                #possibilities1) it can go to smallest size w mw violating options this shouldnt happen
                #possibilities  you want it to go to a size at least higher than the previous one
                try:
                    sol=tmodel.optimize()
                    print('i m here and the solution is {}'.format(sol.objective_value))
                    addIntegerCuts = False
                except SolverError:
                    print('Number of alternatives generated is ',alternative)
                    break


            #to store the values  in wgat type of frame work hdf5


    return list_,df_all