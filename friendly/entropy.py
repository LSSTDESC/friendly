import numpy as np
from scipy.special import xlogy

def entropy_simplex(x: list):
    """
    Computes the entropy of a probability simplex.

    Parameters
    ----------
    x : array-like
        A vector of non-negative values.

    Returns
    -------
    entropy : float
        The Shannon entropy of the normalized probability distribution derived from `x`.
    """

    p = x/np.sum(x) #matching probability
    entropy = -np.sum(xlogy(p,p))
    if entropy == -0.0: entropy = -entropy
    return entropy


def blending_entropy(G):
    """
    Computes the blending entropy for each object node in a NetworkX graph.

    Parameters
    ----------
    G : networkx.Graph
        A graph where nodes represent galaxies (truth sources) and detected objects.
        Each edge contains an `'overlap_fraction'` attribute indicating the
        overlap surface between the object and the galaxy.

    Returns
    -------
    None
        The function modifies the input graph `G` in place by adding the following attributes
        to object nodes:
        
        - `'proba'`: A dictionary mapping adjacent galaxies to their computed probability of matching the detected object.
        - `'blending_entropy'`: The entropy of the probability distribution over galaxies.
    """
    
    gal_nodes = [i for i in G if G.nodes[i]['galaxy']==True]
    obj_nodes = [i for i in G if G.nodes[i]['galaxy']==False]
    
    if len(gal_nodes) == 0: #0-m systems
        for m in obj_nodes:
            G.nodes[m]['proba'] = 0
            G.nodes[m]['blending_entropy'] = -1 
        
    elif len(obj_nodes) == 0: #n-0 systems
        return None
        
    else:
        for m in obj_nodes:

            adj = list(G.adj[m])
            if all(G.nodes[node]['galaxy'] == False for node in adj): #all the adjacent nodes of the object are objects - cannot compute the entropy
                    G.nodes[m]['proba'] = 0
                    G.nodes[m]['blending_entropy'] = -1
                    
            else:
                x = []
                adj_gal = []
                for c in adj:
                    if (c in gal_nodes) & (G[m][c]['overlap_fraction'] != 0):
                        adj_gal.append(str(c))
                        x.append(G[m][c]['overlap_fraction']*np.exp(-np.abs(G.nodes[m]['magnitude'] - G.nodes[c]['magnitude'])))

                G.nodes[m]['proba'] = {gal:proba for gal, proba in zip(adj_gal, x/np.sum(x))}
                G.nodes[m]['blending_entropy'] = entropy_simplex(x)            
    return None
