from .. import Matcher

import FoFCatalogMatching
from ..utils import Group
import numpy as np

class FoF(Matcher):     
    def __init__(self, tune_params: dict):
        self.linking_length = tune_params.get('linking_length', 1.)
        
    def __call__(self, cat1, cat2, *args) -> List[Group]:
        results = FoFCatalogMatching.match(catalog_dict={'cat1':cat1, 'cat2':cat2}, linking_lengths=self.linking_length)

        # first we need to know which rows are from the truth catalog and which are from the object
        mask1 = results['catalog_key'] == 'cat1'
        mask2 = ~mask1
        idx1 = 

        # then np.bincount will give up the number of id occurrences (like historgram but with integer input)
        n_groups = results['group_id'].max() + 1

