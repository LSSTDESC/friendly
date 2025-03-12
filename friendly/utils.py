from collections import namedtuple
from pandas import DataFrame
from astropy.table import Table

Group = namedtuple('Group', ['idx1', 'idx2'])

class FCatalog:
    def __init__(self, cat, pixel_scale=0.2, columns=None):
        """To be replaced with table_io

        Parameters
        ----------
        cat
            _description_
        pixel_scale
            Scale of pixels in arcseconds (LSST is 0.2'' per pixel)
        columns
            Name+Description of columns in catalog
        """
        self.cat = cat
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
                    
    # def get_quantity(self, key, idx):
    #     if isinstance(self.cat, GCRCatalogs):
    #         pass
    #     elif isinstance(self.cat, DataFrame):
    #         self.cat[key][idx]
    #     elif isinstance(self.cat, Table):
    #         self.cat[key][idx]
    #     else:
    #         raise TypeError
                    
