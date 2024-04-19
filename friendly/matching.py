import FoFCatalogMatching
import pandas as pd

""
def FoF_matching(truth_cat: dict, obj_cat: dict, linking_length: float=2.) -> list:
    """
    Friends-of-Friends (FoF) matching between two catalogs.
    Returns the resulted FoF groups.

    Args:
        truth_cat (dict): Truth (galaxy) catalog
        obj_cat (dict): Object catalog
        linking_length (float, optional): Linking length in arcseconds. Defaults to 2.

    Returns:
        list: list of the FoF groups. Each group is composed of two sublists.
              The first sublist corresponds to the row indices of truth catalog components.
              The second sublist corresponds to the row indices of object catalog components.
    """

    FoF_groups = []
    
    FoF_res = FoFCatalogMatching.match(catalog_dict={'galaxy':pd.DataFrame(truth_cat), 'object':pd.DataFrame(obj_cat)}, 
                                  linking_lengths=linking_length).to_pandas()
    
    FoF_res['is_galaxy'] = (FoF_res['catalog_key']=='galaxy')
    res = FoF_res.drop('catalog_key', axis=1)
    
    for group_id, rows in res.groupby('group_id'):  #this takes time
        gal_id = list(rows['row_index'][rows['is_galaxy']])
        obj_id = list(rows['row_index'][~rows['is_galaxy']])
        FoF_groups.append([gal_id, obj_id])
    
    return FoF_groups
