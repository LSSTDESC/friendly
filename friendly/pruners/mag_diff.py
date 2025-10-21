from ..pruner import Pruner
from ..utils import Group

from collections.abc import Iterable
import numpy as np

class MagDiffPruner(Pruner):

    def __init__(self, tune_params: dict):
        self.required_params = ['ground_mag_name', 'space_mag_name']
        self.ground_mag_limit = tune_params.get('ground_mag_limit', 25)
        self.space_mag_limit = tune_params.get('space_mag_limit', 27)
        self.delta_mag_limit = tune_params.get('delta_mag_limit', 2)
    
    def __call__(self, cat1, cat2, params, groups)->list[Group]:
        """
        Prune based on the magnitude difference of the space objects
        after applying basic limiting cuts on both catalogs
        """
        pruned_groups = []

        for group in groups:
            # If there are no space match, skip this entry
            if len(group.idx2)==0:
                pruned_groups.append(group)
                continue

            ground_mags, space_mags = group2mags(group, cat1, cat2, params['ground_mag_name'], params['space_mag_name'])

            # ground_filt = ground_mags < self.ground_mag_limit # The objects we want to keep
            # new_ground_idx = apply_bool_list(group.idx1, ground_filt)
            new_ground_idx = group.idx1

            space_cut1 = space_mags < self.space_mag_limit # The objects we want to keep
            space_cut2 = (space_mags - space_mags.min()) < 2 # The objects we want to keep
            space_filt = space_cut1 & space_cut2
            if not isinstance(space_filt, Iterable):
                space_filt = [space_filt]
            new_space_idx = apply_bool_list(group.idx2, space_filt)

            # if type(space_filt) is np.bool:
            #     new_space_idx = group.idx2 if space_filt else []
            # else:
            #     new_space_idx = apply_bool_list(group.idx2, space_filt)
            
            pruned_groups.append(Group(new_ground_idx, new_space_idx))

        return pruned_groups


def group2mags(group, cat1, cat2, mag_name1, mag_name2):
    ground_mags = cat1.get_quantity(mag_name1, idx=group.idx1, ndx=True)
    space_mags = cat2.get_quantity(mag_name2, idx=group.idx2, ndx=True)
    return ground_mags, space_mags

def apply_bool_list(plist, filt):
    return [int(p) for p,f in zip(plist, filt) if f]



