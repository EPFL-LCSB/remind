
import pandas as pd
import numpy as np

from ..optim.variables import YieldUse, UptakeActivation, SecretionActivation


def get_subs_prods(df):
    #todo can do better with np.close but this shoudl also work
    #after you can screen for subs and prods at the same time

    s=df[df.BU>0.5].metabolites.unique()
    p=df[df.FU>0.5].metabolites.unique()
    substrates=[k for k in s]
    products=[k for k in p]
    substrates.sort()
    products.sort()
    return substrates, products

def get_unique_mets(df,col='metabolites'):
    unique_mets=df.metabolites.unique().tolist()
    unique_mets.sort()
    return unique_mets


def get_substrates(df,col='metabolites'):
    substrates = df[df.BU > 0.5].metabolites.unique().tolist()
    # unique_mets=df[col].unique().tolist()
    substrates.sort()
    return substrates

def get_products(df,col='metabolites'):
    products = df[df.FU > 0.5].metabolites.unique().tolist()
    # unique_mets=df[col].unique().tolist()
    products.sort()
    return products


def get_subs_prods_dIME(df,org_name='Ecoli'):
    #todo can do better with np.close but this shoudl also work
    #after you can screen for subs and prods at the same time

    s_o = df[df['FU or BU'].str.startswith('UA')]['FU or BU'].unique()
    p_o = df[df['FU or BU'].str.startswith('SA')]['FU or BU'].unique()
    s = [k.replace('UA_{}_'.format(org_name),'') for k in s_o]
    p = [k.replace('SA_{}_'.format(org_name),'') for k in p_o]
    substrates=s
    products=p
    substrates.sort()
    products.sort()
    return substrates, products


#also add the total ids
def get_subs_prods_per_group(df,gr=['alternative'],name='subs_prods',option='iME'):
    if option=='iME':
        frame = df.groupby(gr).apply(lambda gr: get_subs_prods(gr)).reset_index(name=name)
    if option=='diME':
        frame = df.groupby(gr).apply(lambda gr: get_subs_prods_dIME(gr)).reset_index(name=name)

    #check this
    frame['id'] = df.groupby(gr).apply(lambda gr: tuple(get_unique_mets(gr))).reset_index(
        name='id')['id']
    frame['subs_list'] = [tuple(frame.iloc[k][name][0]) for k in range(frame.shape[0])]
    frame['prods_list'] = [tuple(frame.iloc[k][name][1]) for k in range(frame.shape[0])]
    frame['alt_size'] = [len(frame.subs_list.iloc[k]) + len(frame.prods_list.iloc[k]) for k in
                               range(frame.shape[0])]

    #alter_ids_o_unique = alter_ids_o.groupby(['subs_list', 'prods_list']).first().reset_index()
    #alter_ids_o_unique.alt_size.value_counts()

    return frame

def get_subs_uptake(df,model):
    '''get the total uptake from flux *mws if
    it was not calculated in the ime calculation'''
    uptakes=[]
    for k in df.alternative.unique():
        d=df[df.alternative==k]
        d_s=d[d.BU>0.5]
        uptake_g=abs(sum([d_s.iloc[k].flux * model.metabolites.get_by_id(d_s.iloc[k].metabolites.replace('EX_', '')).formula_weight for \
           k in range(d_s.shape[0])]))
        uptakes.append(uptake_g*np.ones(d.shape[0]))
    upts=np.concatenate(uptakes)
    df['uptake_mw']=upts
    return df



def get_subs_prods_with_id(df):
    #todo can do better with np.close but this shoudl also work
    #after you can screen for subs and prods at the same time

    substrates=df[df.BU>0.5].metabolites.unique()
    products=df[df.FU>0.5].metabolites.unique()
    substrates=[k.replace('EX_','') for k in substrates]
    products=[k.replace('EX_','') for k in products]

    return substrates, products

def get_subs_prods_with_name(df):
    #todo can do better with np.close but this shoudl also work
    #after you can screen for subs and prods at the same time

    substrates=df[df.BU>0.5].metabolites.unique()
    products=df[df.FU>0.5].metabolites.unique()
    substrates=[salvi.metabolites.get_by_id(k.replace('EX_','')).name for k in substrates]
    products=[salvi.metabolites.get_by_id(k.replace('EX_','')).name for k in products]

    return substrates, products

def get_subs_prods_tva(df,exchanges):
    dfexc=df.loc[exchanges]
    'get all substrates and products'
    #after you can screen for subs and prods at the same time
    subs=dfexc[dfexc.minimum<0].index
    prods = dfexc[dfexc.maximum > 0].index
    return subs, prods

'always needs to be secereted or uptaken '
def get_essential_subs_prods_tva(df,exchanges):
    dfexc=df.loc[exchanges]
    'get all substrates and products'
    subs=dfexc[(dfexc.minimum<0) & (dfexc.maximum<0)].index
    prods = dfexc[(dfexc.maximum > 0) &((dfexc.minimum > 0))].index
    return subs, prods

'adds size of the alternative to the dataframe' \
'aim: easily filter out sizewise what you want'
'group can be two fold if run for different growth rates'
def add_size(df,group=[''],name='label'):
    df[name] = df.groupby(by=group)[group[0]].transform(np.size)



'groupby and for each alternative get substrates products and their numbers'
def get_subs_prods_stats(df):
    #todo can do better with np.close but this shoudl also work
    #after you can screen for subs and prods at the same time
    df_stats=pd.DataFrame()
    substrates=df[df.BU>0.5].metabolites.unique()
    products=df[df.FU>0.5].metabolites.unique()
    df_stats['subs']=[substrates]
    df_stats['prods']=[products]
    df_stats['no_subs']=df[df.BU>0.5].shape[0]
    df_stats['no_prods']=df[df.FU>0.5].shape[0]
    #df_stats['subs']=[substrates]
    #df_stats['prods']=[products]
    return df_stats

'for each group get it'
def get_df_for_subs_prods(df,group=['Growth','alt_size','alternative']):
    df_stats=df.groupby(group).apply(lambda gr: get_subs_prods_stats(gr)).reset_index()
    return df_stats

def get_dime_coupling(data,yield_col='yield'):
    '''
    
    coupling based on grouping
    1)yield
    2)element (if applicable) THIS NEEDS TO BE FIXED!
    3)alternative
    4)substrates and products
    Parameters
    ----------
    data : pd.DataFarme
        The ouptput of iMe function for the alternative DIMEs.
        The columns must contain ['metabolites', 'alternative']

    Returns
    -------
    coupling_dict: dict
        A dict to couple each DIME with the exchange reactions.

    '''

    coupling_dict = dict()
    yield_dict = dict(tuple(data.groupby(yield_col))) # a dict where data grouped for each DIME
    
    for yield_, yield_data in yield_dict.items():
        #todo give here options if element is not present dont do
        
        dime_dict = dict(tuple(yield_data.groupby('alternative')))
        for dime_id, dime in dime_dict.items():
            ret = {}
            ret['reactions'] = dime['metabolites'].tolist()
            ret['substrates'] = dime[np.isclose(dime['BU'],1)]['metabolites'].tolist()
            ret['products'] = dime[np.isclose(dime['FU'],1)]['metabolites'].tolist()
            ret['yield'] = yield_
            coupling_dict['{}_{}'.format(yield_,int(dime_id))] = ret

    return coupling_dict

def get_indispensability_score(data):
    '''
    A fucntion taking an DIME DataFrame and returns a dict of scores proportional to the frequency
    that the metabolites are used in DIMEs. The score is abs(num_secret. - num_upt.), where
    num_secret. and num_upt. are the number of DIMEs in which the metabolite was secreted
    and uptaken, respectively.

    Parameters
    ----------
    data : pd.DataFrame 
        The ouptput of iMe function for the alternative DIMEs.
        The columns must contain ['metabolites', ''FU, 'BU', 'alternative'].

    Returns
    -------
    score_dict : dict
        A dict where the keys are exchange rxn IDs and the values are the scores.

    '''
    score_dict = dict()
    
    rxn_dict = dict(tuple(data.groupby('metabolites'))) # a dict where data grouped for each reaction (metabolite)
    
    for rxn, rxn_data in rxn_dict.items():
        num_sec = len(rxn_data[np.isclose(rxn_data['FU'],1)])
        num_upt = len(rxn_data[np.isclose(rxn_data['BU'],1)])
        indispensability = abs(num_upt - num_sec) # the idea is if a met is sometdimes secreted and sometdimes uptaken, it's not indispensable.
        score_dict[rxn] = {'uptake'         : num_upt,
                           'secretion'     : num_sec,
                           'indispensability': indispensability}
    
    return score_dict

def take_intersect(data_1, data_2):
    '''
    returns the DIMEs that are identical between two datasets
    '''
    
    dime_dict_1 = dict(tuple(data_1.groupby('alternative'))) # a dict where data grouped for each DIME
    dime_dict_2 = dict(tuple(data_2.groupby('alternative'))) # a dict where data grouped for each DIME
    
    ret_1 = []
    ret_2 = []
    
    for dime_id, dime in dime_dict_1.items():
        
        ret_1.append(dime['metabolites'].tolist())
        
    for dime_id, dime in dime_dict_2.items():
        
        ret_2.append(dime['metabolites'].tolist())
        
    intersection_set = [x for x in ret_1 if x in ret_2]
    
    return intersection_set

def get_active_rxn_yield(model, species_id):
    '''
    this functions checks which yield, substrates, and products are active for a specific organism in the community

    Parameters
    ----------
    model : InteractionModel
        DESCRIPTION.
    species_id : str
        DESCRIPTION.

    Returns
    -------
    yield_ : float
        DESCRIPTION.
    substrates : list
        DESCRIPTION.
    products : list
        DESCRIPTION.

    '''
    # check that the model is solved and feasible
    
    if model.problem.__name__ == 'optlang.cplex_interface':
        try:
            assert model.solver.problem.solution.get_status() == 101 #solved and feasible
        except AssertionError:
            raise ValueError('The model must be solved and feasible.')
    if model.problem.__name__ == 'optlang.gurobi_interface':
        try:
            assert model.solver.problem.Status == 2 #solved and feasible
        except AssertionError:
            raise ValueError('The model must be solved and feasible.')
    
    yield_use = model.get_variables_of_type(YieldUse)
    sec_act = model.get_variables_of_type(SecretionActivation)
    upt_act = model.get_variables_of_type(UptakeActivation)
    
    delimiter = species_id + '_'
    # take species-specific variables and check if theye are active
    active_yield = [var for var in yield_use if var.species == species_id \
                    and np.isclose(var.variable.primal,1)]
    substrates = [var.name for var in upt_act if var.species == species_id \
                    and np.isclose(var.variable.primal,1)]
    products = [var.name for var in sec_act if var.species == species_id \
                    and np.isclose(var.variable.primal,1)]
    
    assert len(active_yield) == 1
    yield_ = float(active_yield[0].name.split(delimiter)[-1])
    products = [p.split(delimiter)[-1] for p in products]
    substrates = [s.split(delimiter)[-1] for s in substrates]
    
    return yield_, substrates, products

def remove_redundant_sols(df,groups=['yield','alternative']):
    '''
    If one solution appears at different yield cuts, we can just keep the one that
    appear at the highest yeild cut.

    Parameters
    ----------
    df : pd.DataFrame
        DESCRIPTION.

    Returns
    -------
    data : pd.DataFrame
        DESCRIPTION.

    '''
    
    yield_dict = dict(tuple(df.groupby(groups)))
    
    subs_prod_pair = {}
    last_key2subs_prod = {}
    removable_keys = []
    
    for key,frame in yield_dict.items():
        yield_ = key[0]
        
        # subs_list = [x for x in frame[frame['BU']==1]['metabolites'].unique()]
        # prod_list = [x for x in frame[frame['FU']==1]['metabolites'].unique()]
        # s_p_key = (subs_list, prod_list)

        #should always be sorted
        unique_mets=[x for x in frame['metabolites'].unique()]
        unique_mets.sort()
        s_p_key = tuple(unique_mets)
        
        if s_p_key in subs_prod_pair:
            seen_yield = subs_prod_pair[s_p_key]
            if yield_ > seen_yield:
                subs_prod_pair[s_p_key] = yield_
                removable_keys += [last_key2subs_prod[s_p_key]]
                last_key2subs_prod[s_p_key] = key
            else:
                 removable_keys += [key]   
        else:
            subs_prod_pair[s_p_key] = yield_
            last_key2subs_prod[s_p_key] = key
    
    for k in removable_keys:
        yield_dict.pop(k)
        
    data = pd.concat(yield_dict.values())
    
    return data


#get data for 7 species or more general


def get_species_data(d,models_dict={},len_models=7):

    vars_PI = [k for k in d.variables.unique() if k.startswith('PI_')]
    vars_NI = [k for k in d.variables.unique() if k.startswith('NI_')]

    vars_UA = [k for k in d.variables.unique() if k.startswith('UA_')]
    vars_AV = [k for k in d.variables.unique() if k.startswith('AV_')]

    vars_SA = [k for k in d.variables.unique() if k.startswith('SA_')]
    vars_yield = [k for k in d.variables.unique() if k.startswith('YU_')]

    vars_dimes = [k for k in d.variables.unique() if 'MXV_' in k]

    pos_int = d[d.variables.isin(vars_PI)]
    neg_int = d[d.variables.isin(vars_NI)]
    uptakes_dict={}
    dimes_dict={}
    for key,value in models_dict.items():
        uptakes_dict[value]=[k for k in vars_UA if value in k]
        dimes_dict[value] = [k for k in vars_dimes if 'MXV_{}'.format(value) in k]

    pos_int = d[d.variables.isin(vars_PI)]
    neg_int = d[d.variables.isin(vars_NI)]
    uptakes = d[d.variables.isin(vars_UA)]
    uptakes['upt'] = [k.split('EX_')[1] for k in uptakes.variables]
    secretions = d[d.variables.isin(vars_SA)]
    yield_df = d[d.variables.isin(vars_yield)]
    abiotic = d[d.variables.isin(vars_AV)]

    dime_no = d[d.variables.isin(vars_dimes)]

    grouping = ['alternative']

    if "modelofint" in d.columns:
        grouping = ['alternative', 'modelofint']

    frame_int = pos_int.groupby(grouping).apply(
        lambda gr: tuple([k.replace('PI__', '') for k in gr.variables.unique()])).reset_index(name='pos_int')
    frame_int['neg_int'] = neg_int.groupby(grouping).apply(
        lambda gr: tuple(np.sort([k.replace('NI__', '') for k in gr.variables.unique()]))).reset_index(name='neg_int')[
        'neg_int']

    for key,value in models_dict.items():
        frame_int['{}_uptake'.format(value)] = uptakes.groupby(grouping).apply(lambda gr: tuple(
            np.sort([k.replace('UA_{}_EX_'.format(value), '') for k in gr.variables.unique() if
                     value in k]))).reset_index(
            name='{}_uptake'.format(value))['{}_uptake'.format(value)]
        frame_int['{}_secretions'.format(value)] = secretions.groupby(grouping).apply(lambda gr: tuple(
            np.sort([k.replace('SA_{}_EX_'.format(value), '') for k in gr.variables.unique() if
                     value in k]))).reset_index(
            name='{}_secretion'.format(value))['{}_secretion'.format(value)]

        frame_int['{}_yield'.format(value)] = yield_df.groupby(grouping).apply(lambda gr: float(
            [k.replace('YU_{}_'.format(value), '') for k in gr.variables.unique() if
             value in k][0])).reset_index(
            name='{}_yield'.format(value))['{}_yield'.format(value)]

        frame_int['{}_alt'.format(value)] = dime_no.groupby(grouping).apply(
            lambda gr: float(
                [k.split('_')[-1] for k in gr.variables.unique() if value in k][0])).reset_index(
            name='{}_alt'.format(value))['{}_alt'.format(value)]

        frame_int['{}_dimes'.format(value)] = frame_int['{}_uptake'.format(value)] + frame_int['{}_secretions'.format(value)]



    #general
    frame_int['upt_counts'] = \
        uptakes.groupby(grouping).apply(lambda gr: dict(gr.upt.value_counts())).reset_index(name='uptake_num')[
                'uptake_num']

    if len(vars_AV):
        frame_int['abiotic'] = abiotic.groupby(grouping).apply(
            lambda gr: tuple(np.sort([k.replace('AV__', '') for k in gr.variables.unique()]))).reset_index(
            name='abiotic')[
            'abiotic']


    return frame_int



#trial to get the unique dimes
#after grouping wrt model and yield
# grouped=gr_frame_nu.groupby(['model','yield_perc'])
#
# name_list=[]
# group_nums=[]
# model_name=[]
# yield_perc=[]
# for name, gr,  in grouped:
#     name_list+=[name]
#     model_name+=[gr.model.unique()]
#     yield_perc+=[gr.yield_perc.unique()]
#     group_nums +=[gr.groupby(['subs_list','prods_list']).ngroups]
#     # print(gr.groupby(['subs_list','prods_list']).ngroups)
#
# dd=pd.DataFrame()
# dd['model']=model_name
# dd['yield_perc']=yield_perc
# dd['group_nums']=group_nums