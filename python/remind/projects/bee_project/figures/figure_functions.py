import os
import glob
from sys import argv
import pandas as pd
import numpy as np
import networkx as nx
from remind.utils.postprocessing import *
import matplotlib.pyplot as plt
import math

alpha_chosen= 0.9

yield_cmap = "rainbow"
fig, ax = plt.subplots(1, 2,figsize=(12, 12),gridspec_kw={'width_ratios': [1, 10]})

"this one is to create a video from interactions considering all the dimes for 2 members"
def get_interaction_map(frame, positions, intr,background_black=False):
    """returns the figure with nodes
    input: frame where there is a column for the model including combined dimes for all members in the community of \
    interest"""

    no_mets = frame.metabolites.nunique()
    unique_mets = frame.metabolites.unique()
    unique_mets = [k.replace('EX_', '') for k in unique_mets]
    # create the dictionary after create it for all unique mets
    # pos=get_coordinates_in_circle(no_mets)
    # positions=dict(zip(unique_mets,pos))

    G = nx.DiGraph()
    model_map = {'Snodgrassella_alvi_wkB2GF': 'S.alvi',
                 'Gilliamella_apicola_wkB1GF': 'G.apicola'}
    models = [model_map[model] for model in frame.model.unique()]
    # here to make it automatic
    # for k in range(len(models)):
    #     positions[models[k]]=
    positions[models[0]] = (0.5, 0)
    positions[models[1]] = (-0.5, 0)
    no_models = len(models)
    'different types of nodes'
    G.add_nodes_from(models, size=1000)
    G.add_nodes_from(unique_mets, size=500)
    all_nodes = np.append(models, unique_mets)
    # create a function to make this Digraph
    subs_dict = {}
    prods_dict = {}
    for model in frame.model.unique():
        subs, prods = get_subs_prods(frame.groupby('model').get_group(model))
        subs = [k.replace('EX_', '') for k in subs]
        prods = [k.replace('EX_', '') for k in prods]
        # subs = [k.replace('_e', '') for k in subs]
        # prods = [k.replace('_e', '') for k in prods]
        model = model_map[model]
        subs_dict[model] = subs
        prods_dict[model] = prods
        G.add_edges_from(zip([model] * len(prods), prods), color='darkred')
        G.add_edges_from(zip(subs, [model] * len(subs)), color='seagreen')

    # fig = plt.figure(figsize=(25, 25))
    # if not video:
    fig, ax = plt.subplots(1, 2,figsize=(12, 12),gridspec_kw={'width_ratios': [1, 10]})
    # pos = nx.spring_layout(G, k=2.5)  # positions for all nodes
    # fixed_positions = positions  # dict with two of the positions set
    # get the subdict
    nodes_of_int = [k for k in G.nodes]
    fixed_positions = {k: v for k, v in positions.items() if k in nodes_of_int}

    fixed_nodes = fixed_positions.keys()
    pos = nx.spring_layout(G, pos=fixed_positions, fixed=fixed_nodes)

    nodesize = np.ones(len(G.nodes)) * 2500
    nodesize[0:no_models] = 3000
    colormap = ['lightgray'] * (len(G.nodes))
    colormap[0:no_models] = ['peachpuff'] * no_models

    colors = {'abiotic': 'khaki',
              'pos_int': 'lightgreen',
              'secretions': 'skyblue',
              'model': 'lightgray'
              }

    labels = []
    for k in G.nodes:
        if k in intr.abiotic:
            label = 'abiotic'

        elif k in intr.pos_int:
            label = 'pos_int'
        elif k in (intr.Salvi_secretions + intr.Gapicola_secretions):
            if k not in intr.pos_int:
                label = 'secretions'

        else:
            label = 'model'

        labels.append(label)

    label = dict(zip(G.nodes, labels))

    colormap = [colors[label[k]] for k in G.nodes]
    colors = nx.get_edge_attributes(G, 'color').values()
    nodes_size = nx.get_node_attributes(G, 'size').values()
    nx.draw_networkx_nodes(G, pos, node_size=nodesize, node_color=colormap)
    nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos, edge_color=colors, node_size=nodesize, arrows=True, arrowsize=15,
                           connectionstyle='arc3, rad = 0.15',
                           width=np.ones(len(G.edges)) * 3, alpha=0.3)  #

    plt.title('Alternative {}, S.alvi yield: {} ,G.apicola yield: {}'.format(intr.label_abiotic, intr.Salvi_yield,
                                                                             intr.Gapicola_yield),color="black")

    plt.box(False)
    plt.sca(ax[0])

    df_pivoted = pd.DataFrame(intr[['Gapicola_yield', 'Salvi_yield']], dtype=float)

    # ax1=ax[1]

    plt.imshow(df_pivoted, cmap=yield_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)
    plt.sca(ax[1])
    ax[0].set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                                 df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
    ax[0].set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
                                 df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])

    ax[0].set_yticklabels([str(k.split('_')[0]) for k in df_pivoted.index.to_list()], rotation=90, fontsize=12,
                          style='italic')
    ax[0].set_xticklabels([" " for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12, style='italic')

    for y in range(df_pivoted.shape[0]):
        for x in range(df_pivoted.shape[1]):
            text = (df_pivoted.iloc[y, x])
            if not math.isnan(text):
                text = str(np.round(float(text), 1))
                ax[0].text(x, y, text,
                           horizontalalignment='center',
                           verticalalignment='center',
                           fontsize=15)

    # return fig,ax
    return


# normal version


def get_interaction_map_video(frame, positions, intr,alt_no=0,background_black=False):
    """returns the figure with nodes
    input: frame where there is a column for the model including combined dimes for all members in the community of \
    interest"""

    if background_black:
        fontcolor = 'white'

    else:
        fontcolor = 'black'
    no_mets = frame.metabolites.nunique()
    unique_mets = frame.metabolites.unique()
    unique_mets = [k.replace('EX_', '') for k in unique_mets]
    # create the dictionary after create it for all unique mets
    # pos=get_coordinates_in_circle(no_mets)
    # positions=dict(zip(unique_mets,pos))

    G = nx.DiGraph()
    model_map = {'Snodgrassella_alvi_wkB2GF': 'S.alvi',
                 'Gilliamella_apicola_wkB1GF': 'G.apicola'}
    models = [model_map[model] for model in frame.model.unique()]
    # here to make it automatic
    # for k in range(len(models)):
    #     positions[models[k]]=
    positions[models[0]] = (0.5, 0)
    positions[models[1]] = (-0.5, 0)
    no_models = len(models)
    'different types of nodes'
    G.add_nodes_from(models, size=1000)
    G.add_nodes_from(unique_mets, size=500)
    all_nodes = np.append(models, unique_mets)
    # create a function to make this Digraph
    subs_dict = {}
    prods_dict = {}
    for model in frame.model.unique():
        subs, prods = get_subs_prods(frame.groupby('model').get_group(model))
        subs = [k.replace('EX_', '') for k in subs]
        prods = [k.replace('EX_', '') for k in prods]

        # here it can be screening
        # subs = [k.replace('_e', '') for k in subs]
        # prods = [k.replace('_e', '') for k in prods]
        model = model_map[model]
        subs_dict[model] = subs
        prods_dict[model] = prods
        G.add_edges_from(zip([model] * len(prods), prods), color='darkred')
        G.add_edges_from(zip(subs, [model] * len(subs)), color='seagreen')

    # fig = plt.figure(figsize=(25, 25))
    # if not video:
    #     fig, ax = plt.subplots(1, 2,figsize=(12, 12),gridspec_kw={'width_ratios': [1, 10]})
    # pos = nx.spring_layout(G, k=2.5)  # positions for all nodes
    # fixed_positions = positions  # dict with two of the positions set
    # get the subdict
    nodes_of_int = [k for k in G.nodes]
    fixed_positions = {k: v for k, v in positions.items() if k in nodes_of_int}

    fixed_nodes = fixed_positions.keys()
    pos = nx.spring_layout(G, pos=fixed_positions, fixed=fixed_nodes)

    nodesize = np.ones(len(G.nodes)) * 3000
    nodesize[0:no_models] = 6000
    colormap = ['lightgray'] * (len(G.nodes))
    colormap[0:no_models] = ['peachpuff'] * no_models

    colors = {'abiotic': 'khaki',
              'pos_int': 'lightgreen',
              'secretions': 'skyblue',
              'model': 'lightgray'
              }

    labels = []
    for k in G.nodes:
        if k in intr.abiotic:
            label = 'abiotic'

        elif k in intr.pos_int:
            label = 'pos_int'
        elif k in (intr.Salvi_secretions + intr.Gapic_secretions):
            if k not in intr.pos_int:
                label = 'secretions'

        else:
            label = 'model'

        labels.append(label)

    label = dict(zip(G.nodes, labels))

    colormap = [colors[label[k]] for k in G.nodes]
    colors = nx.get_edge_attributes(G, 'color').values()
    nodes_size = nx.get_node_attributes(G, 'size').values()
    nx.draw_networkx_nodes(G, pos, node_size=nodesize, node_color=colormap,alpha=0.8)
    nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos, edge_color=colors, node_size=nodesize, arrows=True, arrowsize=15,
                           connectionstyle='arc3, rad = 0.15',
                           width=np.ones(len(G.edges)) * 3, alpha=1.0)  #
    #
    plt.title('Alternative {}'.format(alt_no),color=fontcolor)

    plt.box(False)
    plt.sca(ax[0])

    df_pivoted = pd.DataFrame(intr[['Gapicola_yield', 'Salvi_yield']], dtype=float)

    # ax1=ax[1]

    plt.imshow(df_pivoted, cmap=yield_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)
    plt.sca(ax[1])
    ax[0].set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                                 df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
    ax[0].set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
                                 df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])

    ax[0].set_yticklabels([str(k.split('_')[0]) for k in df_pivoted.index.to_list()], rotation=0, fontsize=15,
                          style='italic',color=fontcolor)
    ax[0].set_xticklabels(["Yields " for k in df_pivoted.columns.to_list()], rotation=0, fontsize=15,color=fontcolor)

    for y in range(df_pivoted.shape[0]):
        for x in range(df_pivoted.shape[1]):
            text = (df_pivoted.iloc[y, x])
            if not math.isnan(text):
                text = str(np.round(float(text), 1))
                ax[0].text(x, y, text,
                           horizontalalignment='center',
                           verticalalignment='center',
                           fontsize=15)

    if background_black:
        fig.set_facecolor("#00000F")

    return fig, ax
    # return



"this one is to create a video from interactions considering all the dimes for 2 members"
"to be modified to be colorede differently for pos interactions and others"
def get_interaction_map_positive_int(frame, positions, intr):
    """returns the figure with nodes
    input: frame where there is a column for the model including combined dimes for all members in the community of \
    interest"""
    no_mets = frame.metabolites.nunique()
    unique_mets = frame.metabolites.unique()
    unique_mets = [k.replace('EX_', '') for k in unique_mets]
    # create the dictionary after create it for all unique mets
    # pos=get_coordinates_in_circle(no_mets)
    # positions=dict(zip(unique_mets,pos))

    G = nx.DiGraph()
    model_map = {'Snodgrassella_alvi_wkB2GF': 'S.alvi',
                 'Gilliamella_apicola_wkB1GF': 'G.apicola'}
    models = [model_map[model] for model in frame.model.unique()]
    # here to make it automatic
    # for k in range(len(models)):
    #     positions[models[k]]=
    positions[models[0]] = (0.5, 0)
    positions[models[1]] = (-0.5, 0)
    no_models = len(models)
    'different types of nodes'
    G.add_nodes_from(models, size=1000)
    G.add_nodes_from(unique_mets, size=500)
    all_nodes = np.append(models, unique_mets)
    # create a function to make this Digraph
    subs_dict = {}
    prods_dict = {}
    for model in frame.model.unique():
        subs, prods = get_subs_prods(frame.groupby('model').get_group(model))
        subs = [k.replace('EX_', '') for k in subs]
        prods = [k.replace('EX_', '') for k in prods]
        # subs = [k.replace('_e', '') for k in subs]
        # prods = [k.replace('_e', '') for k in prods]
        model = model_map[model]
        subs_dict[model] = subs
        prods_dict[model] = prods
        G.add_edges_from(zip([model] * len(prods), prods), color='darkred')
        G.add_edges_from(zip(subs, [model] * len(subs)), color='seagreen')

    # fig = plt.figure(figsize=(25, 25))
    # if not video:
    # fig, ax = plt.subplots(1, 2,figsize=(20, 20),gridspec_kw={'width_ratios': [1, 10]})
    # pos = nx.spring_layout(G, k=2.5)  # positions for all nodes
    # fixed_positions = positions  # dict with two of the positions set
    # get the subdict
    nodes_of_int = [k for k in G.nodes]
    fixed_positions = {k: v for k, v in positions.items() if k in nodes_of_int}

    fixed_nodes = fixed_positions.keys()
    pos = nx.spring_layout(G, pos=fixed_positions, fixed=fixed_nodes)

    nodesize = np.ones(len(G.nodes)) * 2500
    nodesize[0:no_models] = 3000
    colormap = ['lightgray'] * (len(G.nodes))
    colormap[0:no_models] = ['peachpuff'] * no_models

    colors = {'abiotic': 'khaki',
              'pos_int': 'lightgreen',
              'secretions': 'skyblue',
              'model': 'lightgray',
              'other':'salmon'
              }

    alphas_ = {'pos_int': 1.0,
              'model': 1.0,
              'other': 0.1,
              }

    labels = []
    for k in G.nodes:
        # if k in intr.abiotic:
        #     label = 'abiotic'
        if k in intr.pos_int:
            label = 'pos_int'

        elif k.count('_')<1:
            label = 'model'

        else:
            label = 'other'

        labels.append(label)

    label = dict(zip(G.nodes, labels))

    alphas=[]
    #here do the same for edges

    for k in G.edges:
        if set(k).intersection(intr.pos_int):
            alpha = 'pos_int'
        else:
            alpha = 'other'

        alphas.append(alpha)

    alpha = dict(zip(G.edges, alphas))


    colormap = [colors[label[k]] for k in G.nodes]
    alpha_for_edges=[alphas_[alpha[k]] for k in G.edges]
    alpha_for_nodes = [alphas_[label[k]] for k in G.nodes]
    colors = nx.get_edge_attributes(G, 'color').values()
    nodes_size = nx.get_node_attributes(G, 'size').values()
    nx.draw_networkx_nodes(G, pos, node_size=nodesize, node_color=colormap,alpha=alpha_for_nodes)
    nx.draw_networkx_labels(G, pos)
    arcs=nx.draw_networkx_edges(G, pos, edge_color=colors, node_size=nodesize, arrows=True, arrowsize=15,
                           connectionstyle='arc3, rad = 0.15',
                           width=np.ones(len(G.edges)) * 3)#,alpha=0.3)  #

    for i, arc in enumerate(arcs):  # change alpha values of arcs
        arc.set_alpha(alpha_for_edges[i])

    plt.title('Alternative {}, S.alvi yield: {} ,G.apicola yield: {}'.format(intr.label_pos_int, intr.Salvi_yield,
                                                                             intr.Gapicola_yield))

    plt.box(False)
    plt.sca(ax[0])

    df_pivoted = pd.DataFrame(intr[['Gapicola_yield', 'Salvi_yield']], dtype=float)

    # ax1=ax[1]

    plt.imshow(df_pivoted, cmap=yield_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)
    plt.sca(ax[1])
    ax[0].set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                                 df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
    ax[0].set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
                                 df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])

    ax[0].set_yticklabels([str(k.split('_')[0]) for k in df_pivoted.index.to_list()], rotation=90, fontsize=12,
                          style='italic')
    ax[0].set_xticklabels([" " for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12, style='italic')

    for y in range(df_pivoted.shape[0]):
        for x in range(df_pivoted.shape[1]):
            text = (df_pivoted.iloc[y, x])
            if not math.isnan(text):
                text = str(np.round(float(text), 1))
                ax[0].text(x, y, text,
                           horizontalalignment='center',
                           verticalalignment='center',
                           fontsize=20)

    # return fig,ax
    return



def get_3_species_data(d):
    # d=d[d.objective==3]
    model_dict={'Gapicola':'Gapic',
                      'Salvi':'Salvi',
                      'Bifido':'Bifido'}
    vars_PI = [k for k in d.variables.unique() if k.startswith('PI_')]
    vars_NI = [k for k in d.variables.unique() if k.startswith('NI_')]

    vars_UA = [k for k in d.variables.unique() if k.startswith('UA_')]
    vars_AV = [k for k in d.variables.unique() if k.startswith('AV_')]

    vars_SA = [k for k in d.variables.unique() if k.startswith('SA_')]
    vars_yield = [k for k in d.variables.unique() if k.startswith('YU_')]

    Gapic_uptake = [k for k in vars_UA if model_dict['Gapicola'] in k]
    Salvi_uptake = [k for k in vars_UA if model_dict['Salvi'] in k]

    Bifido_uptake = [k for k in vars_UA if model_dict['Bifido'] in k]


    vars_dimes = [k for k in d.variables.unique() if 'MXV_' in k]

    vars_dimes_Salvi = [k for k in d.variables.unique() if 'MXV_{}'.format(model_dict['Salvi']) in k]
    vars_dimes_Gapic = [k for k in d.variables.unique() if 'MXV_{}'.format(model_dict['Gapicola']) in k]
    vars_dimes_Bifid = [k for k in d.variables.unique() if 'MXV_{}'.format(model_dict['Bifido']) in k]

    k = d[d.variables.isin(vars_PI + vars_NI)]

    #
    # obj_max=d[d.objective==2]
    # obj_min=d[d.objective==1]
    pos_int = d[d.variables.isin(vars_PI)]
    neg_int = d[d.variables.isin(vars_NI)]
    uptakes = d[d.variables.isin(vars_UA)]
    uptakes['upt'] = [k.split('EX_')[1] for k in uptakes.variables]
    secretions = d[d.variables.isin(vars_SA)]
    yield_df = d[d.variables.isin(vars_yield)]
    abiotic = d[d.variables.isin(vars_AV)]

    dime_no = d[d.variables.isin(vars_dimes)]
    k.groupby('alternative').variables.unique()

    # normally
    grouping = ['alternative']

    if "modelofint" in d.columns:
        grouping = ['alternative', 'modelofint']

    # withn+1

    # extract also yield_use
    frame_int = pos_int.groupby(grouping).apply(
        lambda gr: tuple([k.replace('PI__', '') for k in gr.variables.unique()])).reset_index(name='pos_int')
    frame_int['neg_int'] = neg_int.groupby(grouping).apply(
        lambda gr: tuple(np.sort([k.replace('NI__', '') for k in gr.variables.unique()]))).reset_index(name='neg_int')[
        'neg_int']
    frame_int['Gapic_uptake'] = uptakes.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('UA_{}_EX_'.format(model_dict['Gapicola']), '') for k in gr.variables.unique() if model_dict['Gapicola'] in k]))).reset_index(
        name='Gapicola_uptake')['Gapicola_uptake']
    frame_int['Salvi_uptake'] = uptakes.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('UA_{}_EX_'.format(model_dict['Salvi']), '') for k in gr.variables.unique() if model_dict['Salvi'] in k]))).reset_index(
        name='Salvi_uptake')['Salvi_uptake']

    frame_int['Bifido_uptake'] = uptakes.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('UA_{}_EX_'.format(model_dict['Bifido']), '') for k in gr.variables.unique() if model_dict['Bifido'] in k]))).reset_index(
        name='Bifido_uptake')['Bifido_uptake']

    frame_int['Gapic_secretions'] = secretions.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('SA_{}_EX_'.format(model_dict['Gapicola']), '') for k in gr.variables.unique() if model_dict['Gapicola'] in k]))).reset_index(
        name='Gapicola_secretion')['Gapicola_secretion']
    frame_int['Salvi_secretions'] = secretions.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('SA_{}_EX_'.format(model_dict['Salvi']), '') for k in gr.variables.unique() if model_dict['Salvi'] in k]))).reset_index(
        name='Salvi_secretion')['Salvi_secretion']
    frame_int['Bifido_secretions'] = secretions.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('SA_{}_EX_'.format(model_dict['Bifido']), '') for k in gr.variables.unique() if model_dict['Bifido'] in k]))).reset_index(
        name='Bifido_secretion')['Bifido_secretion']




    frame_int['Salvi_yield'] = yield_df.groupby(grouping).apply(
        lambda gr: float([k.replace('YU_{}_'.format(model_dict['Salvi']), '') for k in gr.variables.unique() if model_dict['Salvi'] in k][0])).reset_index(
        name='Salvi_yield')['Salvi_yield']
    frame_int['Gapicola_yield'] = yield_df.groupby(grouping).apply(lambda gr: float(
        [k.replace('YU_{}_'.format(model_dict['Gapicola']), '') for k in gr.variables.unique() if model_dict['Gapicola'] in k][0])).reset_index(
        name='Gapicola_yield')['Gapicola_yield']

    frame_int['Bifido_yield'] = yield_df.groupby(grouping).apply(lambda gr: float(
        [k.replace('YU_{}_'.format(model_dict['Bifido']), '') for k in gr.variables.unique() if model_dict['Bifido'] in k][0])).reset_index(
        name='Bifido_yield')['Bifido_yield']

    frame_int['upt_counts'] = \
        uptakes.groupby(grouping).apply(lambda gr: dict(gr.upt.value_counts())).reset_index(name='uptake_num')[
            'uptake_num']
    frame_int['Gapicola_alt'] = dime_no.groupby(grouping).apply(
        lambda gr: float([k.split('_')[-1] for k in gr.variables.unique() if model_dict['Gapicola'] in k][0])).reset_index(
        name='Gapicola_alt')['Gapicola_alt']
    frame_int['Salvi_alt'] = dime_no.groupby(grouping).apply(
        lambda gr: float([k.split('_')[-1] for k in gr.variables.unique() if model_dict['Salvi'] in k][0])).reset_index(
        name='Salvi_alt')['Salvi_alt']

    frame_int['Bifido_alt'] = dime_no.groupby(grouping).apply(
        lambda gr: float([k.split('_')[-1] for k in gr.variables.unique() if model_dict['Bifido'] in k][0])).reset_index(
        name='Bifido_alt')['Bifido_alt']

    frame_int['Salvi_dimes'] = frame_int.Salvi_uptake + frame_int.Salvi_secretions
    frame_int['Gapic_dimes'] = frame_int.Gapic_uptake + frame_int.Gapic_secretions
    frame_int['Bifido_dimes'] = frame_int.Bifido_uptake + frame_int.Bifido_secretions

    if len(vars_AV):
        frame_int['abiotic'] = abiotic.groupby(grouping).apply(
            lambda gr: tuple(np.sort([k.replace('AV__', '') for k in gr.variables.unique()]))).reset_index(
            name='abiotic')[
            'abiotic']

    # frame_int['Gapic_dimes'] = frame_int.Gapic_uptake + frame_int.Gapic_secretions
    # frame_int['Salvi_dimes'] = frame_int.Salvi_uptake + frame_int.Salvi_secretions
    # col_of_int = ['Salvi_uptake', 'Gapic_uptake']
    # col_of_int=['Salvi_uptake']

    return frame_int




"get interaction map 3 species"


"this one is to create a video from interactions considering all the dimes for 2 members"
def get_interaction_map_3_species(frame, positions, intr):
    """returns the figure with nodes
    input: frame where there is a column for the model including combined dimes for all members in the community of \
    interest"""
    no_mets = frame.metabolites.nunique()
    unique_mets = frame.metabolites.unique()
    unique_mets = [k.replace('EX_', '') for k in unique_mets]
    # create the dictionary after create it for all unique mets
    # pos=get_coordinates_in_circle(no_mets)
    # positions=dict(zip(unique_mets,pos))

    G = nx.DiGraph()
    model_map = {'Snodgrassella_alvi_wkB2GF': 'S.alvi',
                 'Gilliamella_apicola_wkB1GF': 'G.apicola',
                 'Bifidobacterium_asteroides_PRL2011GF': 'Bifido'}
    models = [model_map[model] for model in frame.model.unique()]
    # here to make it automatic
    # for k in range(len(models)):
    #     positions[models[k]]=
    positions[models[0]] = (0.5, 0)
    positions[models[1]] = (-0.5, 0)
    positions[models[2]] = (0.0, 0)

    no_models = len(models)
    'different types of nodes'
    G.add_nodes_from(models, size=1000)
    G.add_nodes_from(unique_mets, size=500)
    all_nodes = np.append(models, unique_mets)
    # create a function to make this Digraph
    subs_dict = {}
    prods_dict = {}
    for model in frame.model.unique():
        subs, prods = get_subs_prods(frame.groupby('model').get_group(model))
        subs = [k.replace('EX_', '') for k in subs]
        prods = [k.replace('EX_', '') for k in prods]
        # subs = [k.replace('_e', '') for k in subs]
        # prods = [k.replace('_e', '') for k in prods]
        model = model_map[model]
        subs_dict[model] = subs
        prods_dict[model] = prods
        G.add_edges_from(zip([model] * len(prods), prods), color='darkred')
        G.add_edges_from(zip(subs, [model] * len(subs)), color='seagreen')

    # fig = plt.figure(figsize=(25, 25))
    # if not video:
    # fig, ax = plt.subplots(1, 2,figsize=(20, 20),gridspec_kw={'width_ratios': [1, 10]})
    # pos = nx.spring_layout(G, k=2.5)  # positions for all nodes
    # fixed_positions = positions  # dict with two of the positions set
    # get the subdict
    nodes_of_int = [k for k in G.nodes]
    fixed_positions = {k: v for k, v in positions.items() if k in nodes_of_int}

    fixed_nodes = fixed_positions.keys()
    pos = nx.spring_layout(G, pos=fixed_positions, fixed=fixed_nodes)

    nodesize = np.ones(len(G.nodes)) * 2500
    nodesize[0:no_models] = 3000
    colormap = ['lightgray'] * (len(G.nodes))
    colormap[0:no_models] = ['peachpuff'] * no_models

    colors = {'abiotic': 'khaki',
              'pos_int': 'lightgreen',
              'secretions': 'skyblue',
              'model': 'lightgray'
              }

    labels = []
    for k in G.nodes:
        if k in intr.abiotic:
            label = 'abiotic'

        elif k in intr.pos_int:
            label = 'pos_int'
        elif k in (intr.Salvi_secretions + intr.Gapicola_secretions+intr.Bifido_secretions):
            if k not in intr.pos_int:
                label = 'secretions'

        else:
            label = 'model'

        labels.append(label)

    label = dict(zip(G.nodes, labels))

    colormap = [colors[label[k]] for k in G.nodes]
    colors = nx.get_edge_attributes(G, 'color').values()
    nodes_size = nx.get_node_attributes(G, 'size').values()
    nx.draw_networkx_nodes(G, pos, node_size=nodesize, node_color=colormap)
    nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos, edge_color=colors, node_size=nodesize, arrows=True, arrowsize=15,
                           connectionstyle='arc3, rad = 0.15',
                           width=np.ones(len(G.edges)) * 3, alpha=0.3)  #

    plt.title('Alternative {}, S.alvi yield: {} ,G.apicola yield: {}'.format(intr.label_abiotic, intr.Salvi_yield,
                                                                             intr.Gapicola_yield))
    #
    # plt.box(False)
    # plt.sca(ax[0])

    df_pivoted = pd.DataFrame(intr[['Gapicola_yield', 'Salvi_yield','Bifido_yield']], dtype=float)

    # ax1=ax[1]

    plt.imshow(df_pivoted, cmap=yield_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)
    # plt.sca(ax[1])
    # ax[0].set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
    #                              df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
    # ax[0].set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
    #                              df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])
    #
    # ax[0].set_yticklabels([str(k.split('_')[0]) for k in df_pivoted.index.to_list()], rotation=90, fontsize=12,
    #                       style='italic')
    # ax[0].set_xticklabels([" " for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12, style='italic')
    #
    # for y in range(df_pivoted.shape[0]):
    #     for x in range(df_pivoted.shape[1]):
    #         text = (df_pivoted.iloc[y, x])
    #         if not math.isnan(text):
    #             text = str(np.round(float(text), 1))
    #             ax[0].text(x, y, text,
    #                        horizontalalignment='center',
    #                        verticalalignment='center',
    #                        fontsize=15)

    # return fig,ax
    return



def get_3_species_data(d):
    # d=d[d.objective==3]
    model_dict={'Gapicola':'Gapic',
                      'Salvi':'Salvi',
                      'Bifido':'Bifido'}
    vars_PI = [k for k in d.variables.unique() if k.startswith('PI_')]
    vars_NI = [k for k in d.variables.unique() if k.startswith('NI_')]

    vars_UA = [k for k in d.variables.unique() if k.startswith('UA_')]
    vars_AV = [k for k in d.variables.unique() if k.startswith('AV_')]

    vars_SA = [k for k in d.variables.unique() if k.startswith('SA_')]
    vars_yield = [k for k in d.variables.unique() if k.startswith('YU_')]

    Gapic_uptake = [k for k in vars_UA if model_dict['Gapicola'] in k]
    Salvi_uptake = [k for k in vars_UA if model_dict['Salvi'] in k]

    Bifido_uptake = [k for k in vars_UA if model_dict['Bifido'] in k]


    vars_dimes = [k for k in d.variables.unique() if 'MXV_' in k]

    vars_dimes_Salvi = [k for k in d.variables.unique() if 'MXV_{}'.format(model_dict['Salvi']) in k]
    vars_dimes_Gapic = [k for k in d.variables.unique() if 'MXV_{}'.format(model_dict['Gapicola']) in k]
    vars_dimes_Bifid = [k for k in d.variables.unique() if 'MXV_{}'.format(model_dict['Bifido']) in k]

    k = d[d.variables.isin(vars_PI + vars_NI)]

    #
    # obj_max=d[d.objective==2]
    # obj_min=d[d.objective==1]
    pos_int = d[d.variables.isin(vars_PI)]
    neg_int = d[d.variables.isin(vars_NI)]
    uptakes = d[d.variables.isin(vars_UA)]
    uptakes['upt'] = [k.split('EX_')[1] for k in uptakes.variables]
    secretions = d[d.variables.isin(vars_SA)]
    yield_df = d[d.variables.isin(vars_yield)]
    abiotic = d[d.variables.isin(vars_AV)]

    dime_no = d[d.variables.isin(vars_dimes)]
    k.groupby('alternative').variables.unique()

    # normally
    grouping = ['alternative']

    if "modelofint" in d.columns:
        grouping = ['alternative', 'modelofint']

    # withn+1

    # extract also yield_use
    frame_int = pos_int.groupby(grouping).apply(
        lambda gr: tuple([k.replace('PI__', '') for k in gr.variables.unique()])).reset_index(name='pos_int')
    frame_int['neg_int'] = neg_int.groupby(grouping).apply(
        lambda gr: tuple(np.sort([k.replace('NI__', '') for k in gr.variables.unique()]))).reset_index(name='neg_int')[
        'neg_int']
    frame_int['Gapic_uptake'] = uptakes.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('UA_{}_EX_'.format(model_dict['Gapicola']), '') for k in gr.variables.unique() if model_dict['Gapicola'] in k]))).reset_index(
        name='Gapicola_uptake')['Gapicola_uptake']
    frame_int['Salvi_uptake'] = uptakes.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('UA_{}_EX_'.format(model_dict['Salvi']), '') for k in gr.variables.unique() if model_dict['Salvi'] in k]))).reset_index(
        name='Salvi_uptake')['Salvi_uptake']

    frame_int['Bifido_uptake'] = uptakes.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('UA_{}_EX_'.format(model_dict['Bifido']), '') for k in gr.variables.unique() if model_dict['Bifido'] in k]))).reset_index(
        name='Bifido_uptake')['Bifido_uptake']

    frame_int['Gapic_secretions'] = secretions.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('SA_{}_EX_'.format(model_dict['Gapicola']), '') for k in gr.variables.unique() if model_dict['Gapicola'] in k]))).reset_index(
        name='Gapicola_secretion')['Gapicola_secretion']
    frame_int['Salvi_secretions'] = secretions.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('SA_{}_EX_'.format(model_dict['Salvi']), '') for k in gr.variables.unique() if model_dict['Salvi'] in k]))).reset_index(
        name='Salvi_secretion')['Salvi_secretion']
    frame_int['Bifido_secretions'] = secretions.groupby(grouping).apply(lambda gr: tuple(
        np.sort([k.replace('SA_{}_EX_'.format(model_dict['Bifido']), '') for k in gr.variables.unique() if model_dict['Bifido'] in k]))).reset_index(
        name='Bifido_secretion')['Bifido_secretion']




    frame_int['Salvi_yield'] = yield_df.groupby(grouping).apply(
        lambda gr: float([k.replace('YU_{}_'.format(model_dict['Salvi']), '') for k in gr.variables.unique() if model_dict['Salvi'] in k][0])).reset_index(
        name='Salvi_yield')['Salvi_yield']
    frame_int['Gapicola_yield'] = yield_df.groupby(grouping).apply(lambda gr: float(
        [k.replace('YU_{}_'.format(model_dict['Gapicola']), '') for k in gr.variables.unique() if model_dict['Gapicola'] in k][0])).reset_index(
        name='Gapicola_yield')['Gapicola_yield']

    frame_int['Bifido_yield'] = yield_df.groupby(grouping).apply(lambda gr: float(
        [k.replace('YU_{}_'.format(model_dict['Bifido']), '') for k in gr.variables.unique() if model_dict['Bifido'] in k][0])).reset_index(
        name='Bifido_yield')['Bifido_yield']

    frame_int['upt_counts'] = \
        uptakes.groupby(grouping).apply(lambda gr: dict(gr.upt.value_counts())).reset_index(name='uptake_num')[
            'uptake_num']
    frame_int['Gapicola_alt'] = dime_no.groupby(grouping).apply(
        lambda gr: float([k.split('_')[-1] for k in gr.variables.unique() if model_dict['Gapicola'] in k][0])).reset_index(
        name='Gapicola_alt')['Gapicola_alt']
    frame_int['Salvi_alt'] = dime_no.groupby(grouping).apply(
        lambda gr: float([k.split('_')[-1] for k in gr.variables.unique() if model_dict['Salvi'] in k][0])).reset_index(
        name='Salvi_alt')['Salvi_alt']

    frame_int['Bifido_alt'] = dime_no.groupby(grouping).apply(
        lambda gr: float([k.split('_')[-1] for k in gr.variables.unique() if model_dict['Bifido'] in k][0])).reset_index(
        name='Bifido_alt')['Bifido_alt']

    frame_int['Salvi_dimes'] = frame_int.Salvi_uptake + frame_int.Salvi_secretions
    frame_int['Gapic_dimes'] = frame_int.Gapic_uptake + frame_int.Gapic_secretions
    frame_int['Bifido_dimes'] = frame_int.Bifido_uptake + frame_int.Bifido_secretions

    if len(vars_AV):
        frame_int['abiotic'] = abiotic.groupby(grouping).apply(
            lambda gr: tuple(np.sort([k.replace('AV__', '') for k in gr.variables.unique()]))).reset_index(
            name='abiotic')[
            'abiotic']

    # frame_int['Gapic_dimes'] = frame_int.Gapic_uptake + frame_int.Gapic_secretions
    # frame_int['Salvi_dimes'] = frame_int.Salvi_uptake + frame_int.Salvi_secretions
    # col_of_int = ['Salvi_uptake', 'Gapic_uptake']
    # col_of_int=['Salvi_uptake']

    return frame_int




"get interaction map 3 species"


"this one is to create a video from interactions considering all the dimes for 2 members"
def get_interaction_map_3_species_video(frame, positions, intr):
    """returns the figure with nodes
    input: frame where there is a column for the model including combined dimes for all members in the community of \
    interest"""
    no_mets = frame.metabolites.nunique()
    unique_mets = frame.metabolites.unique()
    unique_mets = [k.replace('EX_', '') for k in unique_mets]
    # create the dictionary after create it for all unique mets
    # pos=get_coordinates_in_circle(no_mets)
    # positions=dict(zip(unique_mets,pos))

    G = nx.DiGraph()
    model_map = {'Snodgrassella_alvi_wkB2GF': 'S.alvi',
                 'Gilliamella_apicola_wkB1GF': 'G.apicola',
                 'Bifidobacterium_asteroides_PRL2011GF': 'Bifido'}
    models = [model_map[model] for model in frame.model.unique()]
    # here to make it automatic
    # for k in range(len(models)):
    #     positions[models[k]]=
    positions[models[0]] = (0.5, 0)
    positions[models[1]] = (-0.5, 0)
    positions[models[2]] = (0.0, 0)

    no_models = len(models)
    'different types of nodes'
    G.add_nodes_from(models, size=1000)
    G.add_nodes_from(unique_mets, size=500)
    all_nodes = np.append(models, unique_mets)
    # create a function to make this Digraph
    subs_dict = {}
    prods_dict = {}
    for model in frame.model.unique():
        subs, prods = get_subs_prods(frame.groupby('model').get_group(model))
        subs = [k.replace('EX_', '') for k in subs]
        prods = [k.replace('EX_', '') for k in prods]
        # subs = [k.replace('_e', '') for k in subs]
        # prods = [k.replace('_e', '') for k in prods]
        model = model_map[model]
        subs_dict[model] = subs
        prods_dict[model] = prods
        G.add_edges_from(zip([model] * len(prods), prods), color='darkred')
        G.add_edges_from(zip(subs, [model] * len(subs)), color='seagreen')

    # fig = plt.figure(figsize=(25, 25))
    # if not video:
    # fig, ax = plt.subplots(1, 2,figsize=(20, 20),gridspec_kw={'width_ratios': [1, 10]})
    # pos = nx.spring_layout(G, k=2.5)  # positions for all nodes
    # fixed_positions = positions  # dict with two of the positions set
    # get the subdict
    nodes_of_int = [k for k in G.nodes]
    fixed_positions = {k: v for k, v in positions.items() if k in nodes_of_int}

    fixed_nodes = fixed_positions.keys()
    pos = nx.spring_layout(G, pos=fixed_positions, fixed=fixed_nodes)

    nodesize = np.ones(len(G.nodes)) * 2500
    nodesize[0:no_models] = 3000
    colormap = ['lightgray'] * (len(G.nodes))
    colormap[0:no_models] = ['peachpuff'] * no_models

    colors = {'abiotic': 'khaki',
              'pos_int': 'lightgreen',
              'secretions': 'skyblue',
              'model': 'lightgray'
              }

    labels = []
    for k in G.nodes:
        if k in intr.abiotic:
            label = 'abiotic'

        elif k in intr.pos_int:
            label = 'pos_int'
        elif k in (intr.Salvi_secretions + intr.Gapic_secretions+intr.Bifido_secretions):
            if k not in intr.pos_int:
                label = 'secretions'

        else:
            label = 'model'

        labels.append(label)

    label = dict(zip(G.nodes, labels))

    colormap = [colors[label[k]] for k in G.nodes]
    colors = nx.get_edge_attributes(G, 'color').values()
    nodes_size = nx.get_node_attributes(G, 'size').values()
    nx.draw_networkx_nodes(G, pos, node_size=nodesize, node_color=colormap)
    nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos, edge_color=colors, node_size=nodesize, arrows=True, arrowsize=15,
                           connectionstyle='arc3, rad = 0.15',
                           width=np.ones(len(G.edges)) * 3, alpha=0.3)  #

    plt.title('Alternative {}, S.alvi yield: {} ,G.apicola yield: {}'.format(intr.label_abiotic, intr.Salvi_yield,
                                                                             intr.Gapicola_yield))

    plt.box(False)
    plt.sca(ax[0])

    df_pivoted = pd.DataFrame(intr[['Gapicola_yield', 'Salvi_yield','Bifido_yield']], dtype=float)

    # ax1=ax[1]

    plt.imshow(df_pivoted, cmap=yield_cmap, vmin=0.0, vmax=1.0, alpha=alpha_chosen)
    plt.sca(ax[1])
    ax[0].set_xticks(np.linspace(0, df_pivoted.columns.nunique() - 1,
                                 df_pivoted.columns.nunique()))  # , [str(k) for k in df_pivoted.index.to_list()])
    ax[0].set_yticks(np.linspace(0, df_pivoted.index.nunique() - 1,
                                 df_pivoted.index.nunique()))  # ,[str(k) for k in df_pivoted.columns.to_list()])

    ax[0].set_yticklabels([str(k.split('_')[0]) for k in df_pivoted.index.to_list()], rotation=90, fontsize=12,
                          style='italic')
    ax[0].set_xticklabels([" " for k in df_pivoted.columns.to_list()], rotation=90, fontsize=12, style='italic')

    for y in range(df_pivoted.shape[0]):
        for x in range(df_pivoted.shape[1]):
            text = (df_pivoted.iloc[y, x])
            if not math.isnan(text):
                text = str(np.round(float(text), 1))
                ax[0].text(x, y, text,
                           horizontalalignment='center',
                           verticalalignment='center',
                           fontsize=15)

    # return fig,ax
    return


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
    secretions['secr'] = [k.split('EX_')[1] for k in secretions.variables]

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

    frame_int['secretion_counts'] = \
        secretions.groupby(grouping).apply(lambda gr: dict(gr.secr.value_counts())).reset_index(name='secretion_num')[
                'secretion_num']
    if len(vars_AV):
        frame_int['abiotic'] = abiotic.groupby(grouping).apply(
            lambda gr: tuple(np.sort([k.replace('AV__', '') for k in gr.variables.unique()]))).reset_index(
            name='abiotic')[
            'abiotic']


    return frame_int

