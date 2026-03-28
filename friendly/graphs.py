import pandas as pd
import numpy as np
import networkx as nx
from itertools import combinations
from tqdm.auto import tqdm
import sys

from astropy.table import vstack

sys.path.append('../')
import friendly.ellipses as ellipses
import friendly.gaussians as gaussians
import friendly.entropy as entropy
import friendly.matching as matching
from friendly.utils import Group, FCatalog


def draw_nodes_edges(group: Group, param: dict, infos, link_type: str='ellipse'):
    """
    Initializes a NetworkX graph.
    Draws the nodes from galaxy and object row indices.
    Draws the edges using the ellipticity overlap test if link_type == 'ellipse', otherwise the overlap between two 2D Gaussians i.e. Gaussian overlap.

    Args:
        group (list): Friends-of-Friends group
        param (dict): A, B, C, D, E, F ellipse parameters of galaxies and objects of the group.
        infos (dict): x, y, a, b, theta information of the ellipses
        link_type (str): Overlap definition. 'ellipse' for the overlap between two ellipses, otherwise the overlap between two 2D Gaussian distributions.

    Returns:
        NetworkX graph: NetworkX graph of the corresponding Friends-of-Friends group.

    """
    idx1 = group.idx1
    idx2 = group.idx2    

    G = nx.Graph()
    
    #Nodes
    G.add_nodes_from(idx1, galaxy=False)
    G.add_nodes_from(idx2, galaxy=True)

    if link_type == 'ellipse':
    #Edges
        e_params = ['A', 'B', 'C', 'D', 'E', 'F']
        for i,j in combinations(G, 2):
            if ellipses.is_overlapping(param.loc[i][e_params], param.loc[j][e_params]):
                G.add_edges_from([(j, i)])

    else:
        g_params = ['RA', 'DEC', 'A', 'B', 'THETA']
        for i,j in combinations(G, 2):
            gauss_overlap = gaussians.gaussian_overlap(infos.loc[i][g_params], infos.loc[j][g_params])
            if gauss_overlap > 1e-2:
                G.add_edges_from([(j, i)])

    return G


def add_magnitude(G, cat1: FCatalog, cat2: FCatalog, params: dict ):
    """
    Add galaxy and object magnitudes in the nodes of the NetworkX graph G.

    Args:
        G (NetworkX graph): NetworkX graph corresponding to one Friends-of-Friends group
        cat1 (FCatalog): Ground Catalog
        cat2 (FCatalog): Space Catalog

    Returns:
        None
    """
    gal_nodes = [i for i in G if G.nodes[i]['galaxy']==True]
    obj_nodes = [i for i in G if G.nodes[i]['galaxy']==False]

    if len(gal_nodes) == 0:
        for m in obj_nodes:
            G.nodes[m]['magnitude'] = cat1.get_quantity(params['MAG1'], m, ndx=True)

    elif len(obj_nodes) == 0:
        for n in gal_nodes:
            G.nodes[n]['magnitude'] = cat2.get_quantity(params['MAG2'], n, ndx=True)

    else: 
        for n in gal_nodes:
            G.nodes[n]['magnitude'] = cat2.get_quantity(params['MAG2'], n, ndx=True)
        for m in obj_nodes:
            G.nodes[m]['magnitude'] = cat1.get_quantity(params['MAG1'], m, ndx=True)
    
    return None


def add_blendedness(G, cat1: FCatalog, params: dict):
    """
    Add object blendedness in the corresponding nodes of the NetworkX graph G.

    Args:
        G (NetworkX graph): NetworkX graph corresponding to one Friends-of-Friends group
        obj_cat (dict): Object catalog

    Returns:
        None
    """
    
    obj_nodes = [i for i in G if G.nodes[i]['galaxy']==False]

    if len(obj_nodes) > 0:
        for m in obj_nodes:
            G.nodes[m]['blendedness'] = cat1.get_quantity(params['BLENDEDNESS1'], m, ndx=True)

    else:
        return None
    
    return None


def add_purity(G, infos: dict):
    """
    Calculates the purity of the nodes in a NetworkX graph G based on wether they represent a galaxy or an 
    object. The purity is computed using the overlap between ellipses represented by 2D Gaussian profiles.
    The purity is then added in each node of the graph.

    Args:
        G (NetworkX graph): NetworkX graph corresponding to one Friends-of-Friends group
        infos (dict): x, y, a, b, theta information of the ellipses

    Returns:
        None.
    """

    gal_nodes = [i for i in G if G.nodes[i]['galaxy']==True]
    obj_nodes = [i for i in G if G.nodes[i]['galaxy']==False]
    info_names = ['RA', 'DEC', 'A', 'B', 'THETA']

    if len(gal_nodes) == 0:
        if len(obj_nodes) == 1:
            G.nodes[obj_nodes[0]]['purity'] = 1

        else:
            for n in obj_nodes:
                infos1 = infos.loc[n][info_names]
                int_ = [gaussians.gaussian_square_int(infos1)]
                
                for c in list(G.adj[n]):
                    if G.nodes[c]['galaxy']==False:
                        infos2 = infos.loc[c][info_names]
                        int_.append(gaussians.gaussian_overlap(infos1, infos2))
                
                purity = int_[0]/np.sum(int_)
                G.nodes[n]['purity'] = purity

    elif len(obj_nodes) == 0:
        if len(gal_nodes) == 1:
            G.nodes[gal_nodes[0]]['purity'] = 1

        else:
            for n in gal_nodes:
                infos1 = infos.loc[n][info_names]
                int_ = [gaussians.gaussian_square_int(infos1)]
                
                for c in list(G.adj[n]):
                    if G.nodes[c]['galaxy']==True:
                        infos2 = infos.loc[c][info_names]
                        int_.append(gaussians.gaussian_overlap(infos1, infos2))
                
                purity = int_[0]/np.sum(int_)
                G.nodes[n]['purity'] = purity


    else:
        for flag in [True, False]:
            for n in [i for i in G if G.nodes[i]['galaxy']==flag]:
                infos1 = infos.loc[n][info_names]
                int_ = [gaussians.gaussian_square_int(infos1)]
                
                for c in list(G.adj[n]):
                    if G.nodes[c]['galaxy']==flag:
                        infos2 = infos.loc[c][info_names]
                        int_.append(gaussians.gaussian_overlap(infos1, infos2))
                
                purity = int_[0]/np.sum(int_)
                G.nodes[n]['purity'] = purity
    
    return None


def add_overlap_fraction(G, e_conic, e_params, link_type='ellipse', eps=1):
    """
    Computes and adds the "overlap fraction" attribute to the edges of the NetworkX graph G.
    The overlap fraction is calculated based on the overlap surface of ellipses corresponding to the nodes of G.
    
    Args:
        G (NetworkX graph): NetworkX graph corresponding to one Friends-of-Friends group
        param (dict): A, B, C, D, E, F ellipse parameters
        infos (dict): x, y, a, b, theta information of the ellipses
        link_type (str): Overlap definition. 'ellipse' for the overlap between two ellipses, 
                                              otherwise overlap between two 2D Gaussian distributions.
        Returns:
        None
    """
    if link_type == 'ellipse':
        conic_names = ['A', 'B', 'C', 'D', 'E', 'F']
        if G.number_of_edges() > 0:
            for i, j, data in G.edges(data=True):
    
                #overlap fraction on the corresponded edges
                x_center = np.mean([e_params.loc[i]['RA'], e_params.loc[j]['RA']]) #x mean position of the two ellipses
                y_center = np.mean([e_params.loc[i]['DEC'], e_params.loc[j]['DEC']]) #y mean position of the two ellipses

                # x_center = np.mean([infos[i][0], infos[j][0]]) #x mean position of the two ellipses
                # y_center = np.mean([infos[i][1], infos[j][1]]) #y mean position of the two ellipses
    
                A_total, A_overlap = ellipses.overlap_area_MC(e_conic.loc[i][conic_names], e_conic.loc[j][conic_names],
                                                              [x_center-eps, x_center+eps], [y_center-eps, y_center+eps])
                fraction = A_overlap/A_total
                G[i][j]['overlap_fraction'] = fraction

    else:
        if G.number_of_edges() > 0:
            for i, j, data in G.edges(data=True):
                gauss_overlap = gaussians.gaussian_overlap(infos[i], infos[j])
                G[i][j]['overlap_fraction'] = gauss_overlap/np.sqrt(gaussians.gaussian_square_int(infos[i])*gaussians.gaussian_square_int(infos[j]))
        
    
    return None

def NetworkX_graph(group: list, truth_cat, obj_cat, infos: dict, param: dict, naming_params: dict, link_type: str='ellipse', blend=True):
    """
    Constructs a NetworkX graph representing the relationships between detected objects and
    true galaxies in a given group.

    Parameters
    ----------
    group : list
        Friends-of-Friends group
    truth_cat : pandas.DataFrame
        Truth (galaxy) catalog
    obj_cat : pandas.DataFrame
        Object catalog
    infos : dict
        x, y, a, b, theta information of the ellipses
    param : dict
        A, B, C, D, E, F ellipse parameters
    naming_params: dict
        Column names in FCatalogs for subprocesses
    link_type : str, optional
        Overlap definition. 'ellipse' for the overlap between two ellipses, otherwise the overlap between two 2D Gaussian distributions.
        Default is `'ellipse'`.


    Returns
    -------
    G : networkx.Graph
        A graph where:
        - Nodes represent either galaxies (from `truth_cat`) or detected objects (from `obj_cat`).
        - Edges represent overlaps between detected objects and their closest true galaxies.
        - Node attributes include:
            - `'magnitude'`: The magnitude of the galaxy or detected object.
            - `'blendedness'`: A measure of how blended an object is.
            - `'purity'`: A measure of how isolated an object/galaxy is.
        - Edge attributes include:
            - `'overlap_fraction'`: The fraction of the detected object that overlaps with the
              true galaxy. Computed either from ellipse or 2D Gaussians definitions.
    """

    #Initialization of the graph
    G = draw_nodes_edges(group, param, infos, link_type=link_type)
    
    #Add magnitude and blendedness
    # add_magnitude(G, truth_cat, obj_cat, naming_params)
    add_magnitude(G, obj_cat, truth_cat, naming_params)
    if blend:
        add_blendedness(G, obj_cat, naming_params)
    
    #Add purity
    add_purity(G, infos)
    
    #Add overlap fraction
    add_overlap_fraction(G, param, infos, link_type=link_type)
    
    return G


def refine_groups(G):
    """
    Splits a given graph into its connected components and returns a list of subgraphs.

    Parameters
    ----------
    G : networkx.Graph
        A NetworkX graph, which may contain multiple disconnected components.

    Returns
    -------
    refined_graphs : list of networkx.Graph
        A list of subgraphs, where each subgraph corresponds to a connected component
        of the original graph.
    """
    
    refined_graphs = []
    
    for S in [G.subgraph(c) for c in nx.connected_components(G)]:
        refined_graphs.append(S)
    
    return refined_graphs

def f_group2graph(group, truth_data, object_data, link_type='ellipse'):
    ellipseinfo_params = {'RA1': 'coord_ra', 'DEC1': 'coord_dec', 'A1': 'A_moment', 'B1': 'B_moment', 'THETA1':'THETA_moment',
                      'RA2': 'RA', 'DEC2': 'DEC', 'A2': 'A_ARCSEC', 'B2': 'B_ARCSEC', 'THETA2':'THETA_IMAGE'}

    obj_infos = ellipses.ellipse_infos(group, object_data, truth_data, ellipseinfo_params)
    obj_params = [ellipses.ellipse_parameters(oi) for oi in obj_infos]

    all_infos = vstack(obj_infos)
    all_infos.add_index('ID')
    all_params = vstack(obj_params)
    all_params.add_index('ID')

    graph_params = {'MAG1': 'i_mag', 'MAG2': 'F814_MAG', 'BLENDEDNESS1': 'i_blendedness'}
    g_t = NetworkX_graph(group, truth_data, object_data, all_infos, all_params, naming_params=graph_params, link_type='ellipse')

    S = refine_groups(g_t)
    for s in S:
        entropy.blending_entropy(s)
    return S

def group2graph(group, truth_data, object_data, link_type='ellipse'):
    """
    Constructs NetworkX graph(s) corresponding to an initial FoF group. Refines the group if necessary. Adds the probability of maching and
    the blending entropy.

    Parameters
    ----------
    group : tuple of lists
        Friends-of-Friends group
    truth_data : pandas.DataFrame
        Truth (galaxy) catalog.
    object_data : pandas.DataFrame
        Object catalog.
    link_type : str, optional
        Overlap definition. 'ellipse' for the overlap between two ellipses, 
                             otherwise overlap between two 2D Gaussian distributions.
    """
        
        # Check for doublons in galaxy and object ids. 
        # If it is the case: remove the galaxy
    for m in group[1]:
        if m in group[0]:
            group[0].remove(m)

    if len(group[0]) == 0: #0-m systems
        obj_infos = ellipses.ellipse_infos(group, truth_data, object_data, dc2_type='object')
        obj_param = ellipses.ellipse_parameters(obj_infos)
        
        G = NetworkX_graph(group, truth_data, object_data, obj_infos.T, obj_param.T, link_type=link_type)
        
        entropy.blending_entropy(G)
        return G
        
    elif len(group[1]) == 0: #n-0 systems
        gal_infos = ellipses.ellipse_infos(group, truth_data, object_data, dc2_type='galaxy')
        gal_param = ellipses.ellipse_parameters(gal_infos)
        
        G = NetworkX_graph(group, truth_data, object_data, gal_infos.T, gal_param.T, link_type=link_type)
        
        entropy.blending_entropy(G)
        return G
        
    else:
        obj_infos = ellipses.ellipse_infos(group, truth_data, object_data, dc2_type='object')
        obj_param = ellipses.ellipse_parameters(obj_infos)

        gal_infos = ellipses.ellipse_infos(group, truth_data, object_data, dc2_type='galaxy')
        gal_param = ellipses.ellipse_parameters(gal_infos)

        infos = pd.merge(gal_infos.T, obj_infos.T, how='outer', left_index=True, right_index=True)
        param = pd.merge(gal_param.T, obj_param.T, how='outer', left_index=True, right_index=True)

        G = NetworkX_graph(group, truth_data, object_data, infos, param, link_type=link_type)

        S = refine_groups(G)
        for s in S:
            entropy.blending_entropy(s)
        return S


def friendly_graphs(truth_data, object_data, link_type='ellipse'):
    """
    Constructs a list of graphs of (blended) galaxy-object systems using 
    Friends-of-Friends (FoF) algorithm, shape information and blending entropy calculations.

    Parameters
    ----------
    truth_data : pandas.DataFrame
        Truth (galaxy) catalog.
    
    object_data : pandas.DataFrame
        Object catalog.
    
    link_type : str, optional, default='ellipse'
        Overlap definition. 'ellipse' for the overlap between two ellipses, 
                             otherwise overlap between two 2D Gaussian distributions.

    Returns
    -------
    graphs : list
        A list of NetworkX graphs, where each graph represents a (blended) group of nearby galaxies 
        and detected objects. Each graph includes computed blending entropy (for detected objects only), overlap fraction, magnitudes in i-band,
        purity and blendedness.
    """
    
    graphs = []
    
    #FoF groups
    print('FoF is starting')
    FoF_groups = matching.FoF_matching(truth_data, object_data, 2.)
    print('FoF is finished')
    
    for group in tqdm(FoF_groups):
        graphs.append(group2graph(group, truth_data, object_data, link_type=link_type))
        
    return graphs
