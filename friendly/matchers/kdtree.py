from ..matcher import Matcher
from ..utils import Group

import numpy as np
from scipy.spatial import KDTree

class FKDTree(Matcher):  
    # Where does it make sense to define required_params?
    # Should setup() and init() be two different functions?

    def __init__(self, tune_params: dict):
        self.required_params = ['RA1', 'DEC1', 'RA2', 'DEC2']
        self.radius = tune_params.get('search_rad', 1/(60.**2))


    def __call__(self,cat1, cat2,  params, *args) -> tuple[list[Group], list[float]]:
        """
        Matches cat1 objects to cat2.
        cat1
            Ground/observed based catalog
        cat2
            Space/truth based catalog 
        params
            Dictionary mapping the required params to relevant column names
        """

        coords1 = np.vstack((cat1.get_quantity(params['RA1']), cat1.get_quantity(params['DEC1']))).T
        coords2 = np.vstack((cat2.get_quantity(params['RA2']), cat2.get_quantity(params['DEC2']))).T
        tree1 = KDTree(coords1)
        tree2 = KDTree(coords2)

        space_results = tree1.query_ball_tree(tree2, r=self.radius)
        ground_results = tree1.query_ball_tree(tree1, r=self.radius) 

        groups = []
        weights = []
        for gr, sr in zip(ground_results, space_results):
            ground_ndx = cat1.get_quantity(cat1.ndx_name, gr).tolist()
            space_ndx = cat2.get_quantity(cat2.ndx_name, sr).tolist()
            groups.append(Group(ground_ndx, space_ndx))
            weights.append([1 for i in range(len(sr))]) # Could I just use np.ones_like(sr)?
            # TODO: Reweight such that sum is 1 

        return groups, weights

        

