from .utils import Group
from typing import List

class Matcher:
    arcsec = 1/(60.**2)
    def __init__(self, required_params: List, tune_params: dict):
        """Matcher superclass

        Parameters
        ----------
        required_params
            List of required parameters for a matching script
            
        tune_params
            dict
        """
        self.requried_params = required_params
        
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
    

        
        
    

    
    
