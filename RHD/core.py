# -*- coding: utf-8 -*-
from . import helpers

class Hydro(object):
    """
    Hydro simulation. Holds everything from grid to results
    arguments:
    - grid: Initialised grid object (geometry and values).
    - solver: solver object that we need to evolve the grid.
    """
    def __init__(self, setup, solver):
        super(Hydro, self).__init__()
        self.setup  = setup
        self.solver = solver


    def run(cls):
        """
        Runs the simulation and dumps the data in a results directory.
        Needs an initialised Hydro object to run.
        """
        return True

        
class Setup(object):
    """
    Grid setting and initialisation tools.
    """
    def __init__(self, ncells=100, xL=0., xR=1.):
        super(Setup, self).__init__()
        self.ncells = ncells
        self.xL = xL
        self.xR = xR
            

class Setup_ST(Setup):
    """Shock Tube setup"""
    def __init__(self, rhoL, vL, pL, rhoR, pR, vR, xL=0., xR=1., xmid=0.5, ncells=100):
        Setup.__init__(self, ncells, xL, xR)
        super(Setup, self).__init__()
        self.rhoL = rhoL
        self.vL = vL
        self.pL = pL
        self.rhoR = rhoR
        self.pR = pR
        self.vR = vR
        self.xL = xL
        self.xR = xR
        self.xmid = xmid
        


def create_hydro():
    return Hydro(None, "HLL")