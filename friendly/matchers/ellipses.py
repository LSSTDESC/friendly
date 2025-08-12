from ..matcher import Matcher
from ..utils import Group
from friendly.ellipses import ab2AB_np, is_overlapping

import pandas as pd
import numpy as np
from scipy.spatial import KDTree

class FEllipse(Matcher):
    arcsec = 1/60.**2
    def __init__(self, tune_params: dict):
        self.required_params = ['RA1', 'DEC1', 'RA2', 'DEC2', 'A1', 'B1', 'THETA1', 'A2', 'B2', 'THETA2']
        self.search_radius = tune_params.get('search_rad', 5/(60.**2))

    def __call__(self, cat1, cat2, params, scale=1) -> tuple[list[Group], list[float]]:
        """
        Matches cat1 objects to cat2.
        cat1
            Ground/observed based catalog
        cat2
            Space/truth based catalog 
        params
            Dictionary mapping the required params to relevant column names
        scale=1
            Scaling factor applied to A1,B1,A2,B2
        """

        coords1 = np.vstack((cat1.get_quantity(params['RA1']), cat1.get_quantity(params['DEC1']))).T
        coords2 = np.vstack((cat2.get_quantity(params['RA2']), cat2.get_quantity(params['DEC2']))).T
        tree1 = KDTree(coords1)
        tree2 = KDTree(coords2)

        space_results = tree1.query_ball_tree(tree2, r=self.search_radius)
        ground_results = tree1.query_ball_tree(tree1, r=self.search_radius) 


        ellipse1 = ab2AB_np(cat1.get_quantity(params['RA1']), cat1.get_quantity(params['DEC1']),
                         cat1.get_quantity(params['A1']) * scale, cat1.get_quantity(params['B1']) * scale,
                         cat1.get_quantity(params['THETA1']), sky=True)

        ellipse2 = ab2AB_np(cat2.get_quantity(params['RA2']), cat2.get_quantity(params['DEC2']),
                         cat2.get_quantity(params['A2']) * scale, cat2.get_quantity(params['B2']) * scale,
                         cat2.get_quantity(params['THETA2']), sky=True)

        groups = []
        weights = []
        gndx = np.arange(len(ground_results))
        for gr, sr, gn in zip(ground_results, space_results, gndx):
            g_start = gn
            ground_group = []
            space_group = []

            e1_start = ellipse1[:, g_start]

            if np.isnan(e1_start).any():
                groups.append(Group(ground_group, space_group))
                weights.append([1 for i in range(len(space_group))]) 
                continue

            if len(gr) >= 1:
                for i, g_ndx in enumerate(gr):
                    e2 = ellipse1[:, g_ndx]
                    if np.isnan(e2).any():
                        continue
                    if is_overlapping(e1_start, e2):
                        ids = cat1.get_quantity(cat1.ndx_name, g_ndx)
                        clean_ids = FEllipse.cast_int(ids)
                        ground_group.append(clean_ids)


            if len(sr) >= 1:
                for i, s_ndx in enumerate(sr):
                    e2 = ellipse2[:, s_ndx]
                    if np.isnan(e2).any():
                        continue
                    if is_overlapping(e1_start, e2):
                        ids = cat2.get_quantity(cat2.ndx_name, s_ndx)
                        clean_ids = FEllipse.cast_int(ids)
                        space_group.append(clean_ids)

#             if len(gr) >= 1 and len(sr) >= 1:
#                 g_start = gn
#                 ground_group = []
#                 space_group = []
#                 # Very unhappy with how this is written to guarantee using python int
#                 # instead of np.int
#                 for i, g_ndx in enumerate(gr):
#                     if is_overlapping(ellipse1[:,g_start], ellipse1[:,g_ndx]):
#                         ids = cat1.get_quantity(cat1.ndx_name, g_ndx)
#                         clean_ids = FEllipse.cast_int(ids)
#                         ground_group.append(clean_ids)
# 
#                 for i, s_ndx in enumerate(sr):
#                     if is_overlapping(ellipse1[:, g_start], ellipse2[:, s_ndx]):
#                         ids = cat2.get_quantity(cat2.ndx_name, s_ndx)
#                         clean_ids = FEllipse.cast_int(ids)
#                         space_group.append(clean_ids)

            groups.append(Group(ground_group, space_group))
            weights.append([1 for i in range(len(space_group))]) 

        return groups, weights


    @staticmethod
    def cast_int(ids):
        if isinstance(ids, int):
            return ids
        if isinstance(ids, np.integer):
            return int(ids)
        if isinstance(ids, np.floating):
            return int(ids)
        if isinstance(ids, list):
            return [int(i) for i in ids]
        print(ids, type(ids))
