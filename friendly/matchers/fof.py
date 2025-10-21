from ..matcher import Matcher

import FoFCatalogMatching
from ..utils import Group
import numpy as np
import pandas as pd

class FoF(Matcher):     
    def __init__(self, tune_params: dict):
        self.linking_length = tune_params.get('linking_length', 1.)
        
    def __call__(self, cat1, cat2, params, *args) -> list[Group]:
        results = FoFCatalogMatching.match(catalog_dict={'cat1':cat1.cat, 'cat2':cat2.cat},
            linking_lengths=self.linking_length,
            ra_label=params['RA'], dec_label=params['DEC']).to_pandas()

        results['is_cat1'] = (results['catalog_key']=='cat1')
        res = results.drop('catalog_key', axis=1)

        groups = []
        for group_id, rows in res.groupby('group_id'):
            idx1 = list(rows['row_index'][rows['is_cat1']])
            idx2= list(rows['row_index'][~rows['is_cat1']])
            # Convert the index of the match to the actual ID 
            ground_ndx = cat1.get_quantity(cat1.ndx_name, idx1).tolist()
            space_ndx = cat2.get_quantity(cat2.ndx_name, idx2).tolist()
            groups.append(Group(ground_ndx, space_ndx))

        return groups

