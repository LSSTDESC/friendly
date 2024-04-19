import pandas as pd
import numpy as np
import networkx as nx
from itertools import combinations
from tqdm.auto import tqdm
import sys

sys.path.append('../')
import friendly.ellipses as ellipses
import friendly.gaussians as gaussians
import friendly.entropy as entropy
import friendly.matching as matching

def draw_nodes_edges_ell(group: list, param: dict):
    """
    Initializes a NetworkX graph.
    Draws the nodes from galaxy and object row indices.
    Draws the edges using the ellipticity overlap test.

    Args:
        group (list): Friends-of-Friends group
        param (dict): A, B, C, D, E, F ellipse parameters of galaxies and objects of the group.

    Returns:
        NetworkX graph: NetworkX graph of the corresponding Friends-of-Friends group.

    """
    idx1 = group[0]
    idx2 = group[1]
    
    G = nx.Graph()
    
    #Nodes
    G.add_nodes_from(idx1, galaxy=True)
    G.add_nodes_from(idx2, galaxy=False)
    
    #Edges
    for i,j in combinations(G, 2):
        if ellipses.is_overlapping(param[i], param[j]):
            G.add_edges_from([(j, i)])

    return G


def add_magnitude(G, truth_cat: dict, obj_cat: dict):
    """
    Add galaxy and object magnitudes in the nodes of the NetworkX graph G.

    Args:
        G (NetworkX graph): NetworkX graph corresponding to one Friends-of-Friends group
        truth_cat (dict): Truth (galaxy) catalog
        object_cat (dict): Object catalog

    Returns:
        None
    """
    gal_nodes = [i for i in G if G.nodes[i]['galaxy']==True]
    obj_nodes = [i for i in G if G.nodes[i]['galaxy']==False]

    if len(gal_nodes) == 0:
        for m in obj_nodes:
            G.nodes[m]['magnitude'] = obj_cat['mag_i_cModel'][m]

    elif len(obj_nodes) == 0:
        for n in gal_nodes:
            G.nodes[n]['magnitude'] = truth_cat['mag_i'][n]

    else: 
        for n in gal_nodes:
            G.nodes[n]['magnitude'] = truth_cat['mag_i'][n]
        for m in obj_nodes:
            G.nodes[m]['magnitude'] = obj_cat['mag_i_cModel'][m]
    
    return None


def add_blendedness(G, obj_cat: dict):
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
            G.nodes[m]['blendedness'] = obj_cat['blendedness'][m]

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

    if len(gal_nodes) == 0:
        if len(obj_nodes) == 1:
            G.nodes[obj_nodes[0]]['purity'] = 1

        else:
            for n in obj_nodes:
                infos1 = infos[n]
                int_ = [gaussians.gaussian_square_int(infos1)]
                
                for c in list(G.adj[n]):
                    if G.nodes[c]['galaxy']==False:
                        infos2 = infos[c]
                        int_.append(gaussians.gaussian_overlap(infos1, infos2))
                
                purity = int_[0]/np.sum(int_)
                G.nodes[n]['purity'] = purity

    elif len(obj_nodes) == 0:
        if len(gal_nodes) == 1:
            G.nodes[gal_nodes[0]]['purity'] = 1

        else:
            for n in gal_nodes:
                infos1 = infos[n]
                int_ = [gaussians.gaussian_square_int(infos1)]
                
                for c in list(G.adj[n]):
                    if G.nodes[c]['galaxy']==True:
                        infos2 = infos[c]
                        int_.append(gaussians.gaussian_overlap(infos1, infos2))
                
                purity = int_[0]/np.sum(int_)
                G.nodes[n]['purity'] = purity


    else:
        for flag in [True, False]:
            for n in [i for i in G if G.nodes[i]['galaxy']==flag]:
                infos1 = infos[n]
                int_ = [gaussians.gaussian_square_int(infos1)]
                
                for c in list(G.adj[n]):
                    if G.nodes[c]['galaxy']==flag:
                        infos2 = infos[c]
                        int_.append(gaussians.gaussian_overlap(infos1, infos2))
                
                purity = int_[0]/np.sum(int_)
                G.nodes[n]['purity'] = purity
    
    return None



def add_overlap_fraction(G, param, infos):
    """
    Computes and adds the "overlap fraction" attribute to the edges of the NetworkX graph G.
    The overlap fraction is calculated based on the overlap surface of ellipses corresponding to the nodes of G.
    
    Args:
        G (NetworkX graph): NetworkX graph corresponding to one Friends-of-Friends group
        param (dict): A, B, C, D, E, F ellipse parameters
        infos (dict): x, y, a, b, theta information of the ellipses

    Returns:
        None
    """
    
    if G.number_of_edges() > 0:
        for i, j, data in G.edges(data=True):

            #overlap fraction on the corresponded edges
            eps = 100
            x_center = np.mean([infos[i][0], infos[j][0]]) #x mean position of the two ellipses
            y_center = np.mean([infos[i][1], infos[j][1]]) #y mean position of the two ellipses

            A_total, A_overlap = gaussians.overlap_area_MC(param[i], param[j], [x_center-eps, x_center+eps], [y_center-eps, y_center+eps])
            fraction = A_overlap/A_total
            G[i][j]['overlap_fraction'] = fraction
    
    return None


def NetworkX_graph(group, truth_cat, obj_cat, infos, param):

    #Initialization of the graph
    G = draw_nodes_edges_ell(group, param)
    
    #Add magnitude and blendedness
    add_magnitude(G, truth_cat, obj_cat)
    add_blendedness(G, obj_cat)
    
    #Add purity
    add_purity(G, infos)
    
    #Add overlap fraction
    add_overlap_fraction(G, param, infos)
    
    return G


def refine_groups(G):
    
    refined_graphs = []
    
    for S in [G.subgraph(c) for c in nx.connected_components(G)]:
        refined_graphs.append(S)
    
    return refined_graphs


def group2graph(group, truth_data, object_data):
        
        # Check for doublons in galaxy and object ids. 
        # If it is the case: remove the galaxy
        for m in group[1]:
            if m in group[0]:
                group[0].remove(m)
    
        if len(group[0]) == 0: #0-m systems
            obj_infos = ellipses.ellipse_infos(group, truth_data, object_data, dc2_type='object')
            obj_param = ellipses.ellipse_parameters(obj_infos)
            
            G = NetworkX_graph(group, truth_data, object_data, obj_infos.T, obj_param.T)
            
            entropy.blending_entropy(G)
            return G
            
        elif len(group[1]) == 0: #n-0 systems
            gal_infos = ellipses.ellipse_infos(group, truth_data, object_data, dc2_type='galaxy')
            gal_param = ellipses.ellipse_parameters(gal_infos)
            
            G = NetworkX_graph(group, truth_data, object_data, gal_infos.T, gal_param.T)
            
            entropy.blending_entropy(G)
            return G
            
        else:
            obj_infos = ellipses.ellipse_infos(group, truth_data, object_data, dc2_type='object')
            obj_param = ellipses.ellipse_parameters(obj_infos)

            gal_infos = ellipses.ellipse_infos(group, truth_data, object_data, dc2_type='galaxy')
            gal_param = ellipses.ellipse_parameters(gal_infos)

            infos = pd.merge(gal_infos.T, obj_infos.T, how='outer', left_index=True, right_index=True)
            param = pd.merge(gal_param.T, obj_param.T, how='outer', left_index=True, right_index=True)

            G = NetworkX_graph(group, truth_data, object_data, infos, param)

            S = refine_groups(G)
            for s in S:
                entropy.blending_entropy(s)
            return S
        
        
def friendly_graphs(truth_data, object_data):
    
    graphs = []
    
    #FoF groups
    print('FoF is starting')
    FoF_groups = matching.FoF_matching(truth_data, object_data, 2.)
    print('FoF is finished')
    
    for group in tqdm(FoF_groups):
        graphs.append(group2graph(group, truth_data, object_data))
        
    return graphs