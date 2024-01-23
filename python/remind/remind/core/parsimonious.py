from pytfa.optim.utils import symbol_sum
import re
import pandas as pd


'A function to minimize the uptake of a specific element'
def minimize_element_upt(model, list_resources):
    rxn_vars = [model.reactions.get_by_id(k.id).reverse_variable for k in list_resources]
    expression = symbol_sum(rxn_vars)
    model.objective = expression
    model.objective_direction = 'min'
    sol = model.optimize()
    return expression , sol


from pytfa.optim.utils import symbol_sum
import re


' if you want to constrain carbon sources'
def constrainCarbonSources(model) :
    listPotCarbonSources=[]
    for r in model.reactions:
        if 'EX_' in r.id:
            for m in r.metabolites:
                if 'C' in m.formula:
                    listPotCarbonSources.append(r.id)
                    model.reactions.get_by_id(r.id).lower_bound=-20
    return model



#TODO here it looks for exchanges through EX_ annotation might be misleading in some cases!
'get a list of potential carbon sources from exchange reactions'
def listcarbonsources(model) :
    listPotCarbonSources=[]
    # you dont want Cl or cOBALT but c followed by a capital letter or number
    pattern = r"[C]\d|[C][A-Z]"
    for r in model.reactions:
        if 'EX_' in r.id:
            for m in r.metabolites:
                result=re.findall(pattern,m.formula)
                if any(result):
                    listPotCarbonSources.append(r.id)
                   # model.reactions.get_by_id(r.id).lower_bound=-20
    return listPotCarbonSources



#TODO here it looks for exchanges through EX_ annotation might be misleading in some cases!
'get a list of potential carbon sources from exchange reactions'
def listcarbonsources_all_rxns(model) :
    listPotCarbonSources=[]
    # you dont want Cl or cOBALT but c followed by a capital letter or number
    pattern = r"[C]\d|[C][A-Z]"
    for r in model.reactions:
        for m in r.metabolites:
            result=re.findall(pattern,m.formula)
            if any(result):
                listPotCarbonSources.append(r.id)
                   # model.reactions.get_by_id(r.id).lower_bound=-20
    return listPotCarbonSources


#get a list of carbon sources from a given dataframe which has metabolite\
# and a formula column
def listcarbonsources_from_df(df,met_col='mets',formula_col='formula') :
    listPotCarbonSources=[]
    # you dont want Cl or cOBALT but c followed by a capital letter or number
    pattern = r"[C]\d|[C][A-Z]"
    for met in df[met_col].unique():
        formula=df[df[met_col]==met][formula_col].iloc[0]
        result=re.findall(pattern,formula)
        if any(result):
            listPotCarbonSources.append(met)
                   # model.reactions.get_by_id(r.id).lower_bound=-20

    return listPotCarbonSources
def listnitrogensources(model) :
    listPotNitrogenSources=[]
    # you dont want Cl or cOBALT but c followed by a capital letter or number
    pattern = r"[N]\d|[N][A-Z]"
    for r in model.reactions:
        if 'EX_' in r.id:
            for m in r.metabolites:
                result=re.findall(pattern,m.formula)
                if any(result):
                    listPotNitrogenSources.append(r.id)
                   # model.reactions.get_by_id(r.id).lower_bound=-20

    return listPotNitrogenSources



#todo to be fixed to a more general notation now its EXC kind of
def listcarbonsources_community(model) :
    listPotCarbonSources=[]
    # you dont want Cl or cOBALT but c followed by a capital letter or number
    pattern = r"[C]\d|[C][A-Z]"
    for r in model.reactions:
        if 'EXC_' in r.id:
            for m in r.metabolites:
                result=re.findall(pattern,m.formula)
                if any(result):
                    listPotCarbonSources.append(r.id)
                   # model.reactions.get_by_id(r.id).lower_bound=-20

    return listPotCarbonSources



#under development
#here idea is to use ids as they ll be fe2_e fe3_e
def list_ions(model):
    listions=[]
    # you dont want Cl or cOBALT but c followed by a capital letter or number
    pattern = r"^[a-z]{2}\d$|^[a-z]{2}$"
    pattern = r"^[A-Z][a-z][+]|^[A-Z][a-z]\d[+]"
    for r in model.reactions:
        if 'EX_' in r.id:
            for m in r.metabolites:
                result=re.findall(pattern,m.name)
                if any(result):
                    listions.append(r.id)
    return listions


'minimize uptake of carbon sources'
'you need to minimize reverse variable for uptakes'
def minimizecarbonuptakes(model,listcarbonsources):
    rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]
    expression=symbol_sum(rxn_vars)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol


def minimizealluptakes(model):
    rxn_vars=[k.reverse_variable for k in model.exchanges]
    expression=symbol_sum(rxn_vars)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol


#this function is to do max yield constraint  on considering only moles of C
def minimizecarbonuptakes_cmoles_only(model,list_carbon_sources):
   # rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]

    #here we use a different patern because we want to track all digits
    pattern = r"[C]\d{1,3}|[C][A-Z]"
    rxns=[model.reactions.get_by_id(k) for k in list_carbon_sources]
    expr=[]
    for this_reaction in rxns:
        for m in this_reaction.metabolites:
            result=re.findall(pattern,m.formula)
            #returns a list access the string by [0]
            #if the next  number is digit remove C and take the float numver
            if result[0][1].isdigit():
               # print('yes')
                C_number=result[0].replace('C','')
                expr.append(float(C_number)*this_reaction.reverse_variable)
              # 'this means it has H or other molecule after so C number is 1'
            else:
                expr.append(1*this_reaction.reverse_variable)

    expression=symbol_sum(expr)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol


#this function is to do max yield constraint  on considering only moles of N
def minimizecarbonuptakes_n_moles_only(model,list_nitrogen_sources):
   # rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]

    #here we use a different patern because we want to track all digits
    pattern = r"[N]\d{1,3}|[N][A-Z]"
    rxns=[model.reactions.get_by_id(k) for k in list_nitrogen_sources]
    expr=[]
    for this_reaction in rxns:
        for m in this_reaction.metabolites:
            result=re.findall(pattern,m.formula)
            #returns a list access the string by [0]
            #if the next  number is digit remove C and take the float numver
            if result[0][1].isdigit():
               # print('yes')
                C_number=result[0].replace('N','')
                expr.append(float(C_number)*this_reaction.reverse_variable)
              # 'this means it has H or other molecule after so C number is 1'
            else:
                expr.append(1*this_reaction.reverse_variable)

    expression=symbol_sum(expr)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol


def minimizecarbonuptakes_cmoles_only_high_C(model,listcarbonsources,molecular_weight='formula_weight'):
   # rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]

    #here we use a different patern because we want to track all digits
    pattern = r"[C]\d{1,3}|[C][A-Z]"
    #pattern2 = r"[C]\d{1,3}|[C][A-Z]
    rxns=[model.reactions.get_by_id(k) for k in listcarbonsources]
    expr=[]
    list_sources=[]
    for this_reaction in rxns:
        for m in this_reaction.metabolites:
            result=re.findall(pattern,m.formula)
            #returns a list access the string by [0]
            #if the next  number is digit remove C and take the float numver
            if result[0][1].isdigit():
               # print('yes')
                C_number=result[0].replace('C','')
                if float(C_number)>6:
                    #expr.append(float(C_number)*this_reaction.reverse_variable)
                    expr.append(getattr(m,molecular_weight)*this_reaction.reverse_variable)
                    list_sources.append(this_reaction.id)
              # 'this means it has H or other molecule after so C number is 1'

    expression=symbol_sum(expr)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol,list_sources


def minimizecarbonandnitrogen_sources(model,listcarbonsources,molecular_weight='formula_weight'):
     # rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]
    #here we use a different patern because we want to track all digits
    pattern1 = r"[C]\d{1,3}|[C][A-Z]"
    pattern2 = r"[N]\d{1,3}|[N][A-Z]"
    rxns=[k for k in model.exchanges]
    expr=[]
    check=[]
    name=[]
    for this_reaction in rxns:
        for m in this_reaction.metabolites:
            result1=re.findall(pattern1,m.formula)
            result2=re.findall(pattern2,m.formula)
            #returns a list access the string by [0]
            #if the next  number is digit remove C and take the float numver
            if any(result1) and any(result2):
               # print('yes')

                expr.append(getattr(m,molecular_weight)*this_reaction.reverse_variable)
              # 'this means it has H or other molecule after so C number is 1'
                check.append(m.formula)
                name.append(m.name)
    expression=symbol_sum(expr)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol



def listcarbonandnitrogensources(model):
    # rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]
    # here we use a different patern because we want to track all digits
    pattern1 = r"[C]\d{1,3}|[C][A-Z]"
    pattern2 = r"[N]\d{1,3}|[N][A-Z]"
    rxns = [k for k in model.exchanges]
    expr = []
    check = []
    name = []
    list_sources=[]
    for this_reaction in rxns:
        for m in this_reaction.metabolites:
            result1 = re.findall(pattern1, m.formula)
            result2 = re.findall(pattern2, m.formula)
            # returns a list access the string by [0]
            if any(result1) and any(result2):
                # print('yes')
                list_sources.append(this_reaction.id)

    return list_sources


def minimize_uptake_list(model,list_sources,molecular_weight='formula_weight'):
    rxns=[model.reactions.get_by_id(k) for k in list_sources]
    expr=[]
    for this_reaction in rxns:
        for m in this_reaction.metabolites:
            expr.append(getattr(m, molecular_weight) * this_reaction.reverse_variable)

    expression=symbol_sum(expr)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol

'for carveme models it is formula_weight'
#todo put these as options to all uptakes and c uptakes function
#this function is to do max yield constraint  on considering only moles of C
def minimizecarbonuptakes_by_weight(model,listcarbonsources,molecular_weight='formula_weight'):
   # rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]

    #here we use a different patern because we want to track all digits
    rxns=[model.reactions.get_by_id(k) for k in listcarbonsources]
    expr=[]
    for this_reaction in rxns:
        for m in this_reaction.metabolites:
             expr.append(getattr(m,molecular_weight)*this_reaction.reverse_variable)

    expression=symbol_sum(expr)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol



#take all uptakes and penalize each substrate by its molecular weight
def minimizealluptakes_by_weight(model, reaction_list=None, molecular_weight='formula_weight'):
   # rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]
    #here we use a different patern because we want to track all digits
   # rxns=[model.reactions.get_by_id(k) for k in listcarbonsources]
    expr=[]
    if reaction_list is None:
        reaction_list = model.exchanges
   #todo beware of excanges
    for this_reaction in reaction_list:
        for m in this_reaction.metabolites:
            expr.append(getattr(m,molecular_weight)*this_reaction.reverse_variable)

    expression=symbol_sum(expr)
    model.objective=expression
    model.objective_direction='min'
    sol=model.optimize()
    return expression,sol


def get_C_number (frame,met_column='mets',formula_column='formula'):
   # rxn_vars = [model.reactions.get_by_id(k).reverse_variable for k in listcarbonsources]
    #here we use a different patern because we want to track all digits
    pattern = r"[C]\d{1,3}|[C][A-Z]"
    mets=[k for k in frame[met_column].unique()]
    expr=[]
    for this_met in mets:
        formula=frame[frame[met_column]==this_met][formula_column].iloc[0]
        result=re.findall(pattern,formula)
        if any(result):
            if result[0][1].isdigit():

                C_number = result[0].replace('C', '')
                expr.append(float(C_number))

            else:
                expr.append(1)
        else:
            expr.append(0)


    fr=pd.DataFrame()
    fr['mets']=mets
    fr['c_num'] = expr

    return fr




def get_C_number_from_an_exch_rxn(this_reaction):
    pattern = r"[C]\d{1,3}|[C][A-Z]"
    for m in this_reaction.metabolites:
        result=re.findall(pattern,m.formula)
        #returns a list access the string by [0]
        #if the next  number is digit remove C and take the float numver
        if any(result):
            if result[0][1].isdigit():
               # print('yes')
                C_number=float(result[0].replace('C',''))

            else:
                C_number=1.0


    return C_number






