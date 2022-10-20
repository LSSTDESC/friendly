from collections import namedtuple
import GCRCatalogs
from pandas import DataFrame

Group = namedtuple('Group', ['idx1', 'idx2'])

class Catalog:
    def __init__(self, cat):
        """To be replaced with table_io

        Parameters
        ----------
        cat
            _description_
        """
        self.cat = cat
    
    def get_quantity(self, key, idx):
        if isinstance(self.cat, GCRCatalogs):
            pass
        elif isinstance(self.cat, DataFrame):
            self.cat[key][idx]
        else:
            raise TypeError
                    