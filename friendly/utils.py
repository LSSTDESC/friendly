from collections import namedtuple
from pandas import DataFrame
from astropy.table import Table
import numpy as np

Group = namedtuple('Group', ['idx1', 'idx2'])

class FCatalog:
    def __init__(self, cat, ndx_name, pixel_scale=0.2, columns=None):
        """To be replaced with table_io

        Parameters
        ----------
        cat
            _description_
        ndx_name
            Column name with index 
        pixel_scale
            Scale of pixels in arcseconds (LSST is 0.2'' per pixel)
        columns
            Name+Description of columns in catalog
        """
        self.cat = cat
        self.ndx_name = ndx_name
        self.pixel_scale = pixel_scale
        self.columns = columns
    
    def get_quantity(self, key, idx=None):
        if isinstance(self.cat, DataFrame):
            if idx is None: #There should be a cleaner way to do this...
                return self.cat[key]
            else:
                return self.cat[key][idx]
        elif isinstance(self.cat, Table):
            if idx is None: #There should be a cleaner way to do this...
                return self.cat[key]
            else:
                return self.cat[key][idx]
        else:
            raise TypeError

    def rename_col(self, old_colname, new_colname):
        """
        Rename column

        Currently just adds a new column instead of overwriting the previous one.
            Should discuss what to do here
        """
        if isinstance(self.cat, Table):
            self.cat[new_colname] = self.cat[old_colname]
            return None
        elif isinstance(self.cat, DataFrame):
            self.cat = self.cat.assign(**{new_colname: self.cat[old_colname].values})
            return None
        else:
            raise TypeError

    def __len__(self):
        # return np.shape(self.cat)[0]
        return len(self.cat)

    # def get_quantity(self, key, idx):
    #     if isinstance(self.cat, GCRCatalogs):
    #         pass
    #     elif isinstance(self.cat, DataFrame):
    #         self.cat[key][idx]
    #     elif isinstance(self.cat, Table):
    #         self.cat[key][idx]
    #     else:
    #         raise TypeError
                    
