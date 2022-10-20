from .utils import Group
from typing import List

class Matcher:
    def __init__(self, tune_params: dict):
        """Matcher superclass

        Parameters
        ----------
        tune_params
            dict
        """
        
    def __call__(self, cat1, cat2, *args) -> List[Group]:
        """match catalogs given the matching method and tuning parameter(s)

        Parameters
        ----------
        cat1
            _description_
        cat2
            _description_
        """
        raise(NotImplementedError)
    
    
class FoF(Matcher):        
    pass
    
class KDTree(Matcher):  
    pass 
        
        
    

    
    