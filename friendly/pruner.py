from .utils import Group
from typing import List
    
class Pruner:
    def __init__(self, tune_params:dict) -> None:
        self.tune_params = tune_params
    
    def __call__(self, cat1, cat2, groups)->List[Group]:
        """_summary_

        Parameters
        ----------
        cat1
            _description_
        cat2
            _description_
        groups
            initial matches as list of Group objects
        """
        raise NotImplementedError
    

class EllipseOverlap(Pruner):
    def __init__(self):
        pass