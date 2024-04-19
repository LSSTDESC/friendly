import numpy as np
from scipy.special import xlogy

def entropy_simplex(x):
    p = x/np.sum(x) #probability of overlap
    entropy = -np.sum(xlogy(p,p)) #cost for one object
    if entropy == -0.0: entropy = -entropy
    return entropy


def blending_entropy(G):
    
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

