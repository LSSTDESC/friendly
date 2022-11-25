from .. import Matcher

import FoFCatalogMatching
from ..utils import Group
import numpy as np

class FoF(Matcher):     
    def __init__(self, tune_params: dict):
        self.linking_length = tune_params.get('linking_length', 1.)
        
    def __call__(self, cat1, cat2, *args) -> List[Group]:
        results = FoFCatalogMatching.match(catalog_dict={'cat1':cat1, 'cat2':cat2}, linking_lengths=self.linking_length)
        results['is_object'] = (results['catalog_key']=='object')
        res = results.drop('catalog_key', axis=1)

        groups = []
        for group_id, rows in res.groupby('group_id'):
            id_gal = list(rows['row_index'][~rows['is_object']])
            id_obj = list(rows['row_index'][rows['is_object']])
            groups.append([id_gal, id_obj])

