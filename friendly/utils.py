from collections import namedtuple
import GCRCatalogs

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
        else:
            raise TypeError
                    