# -*- coding: utf-8 -*-
from . import helpers

import numpy as np

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


class State(object):
    """docstring for State"""
    def __init__(self, rho, v, p, D, m, E, lfac, h, cs):
        super(State, self).__init__()
        self.rho = rho
        self.v = v
        self.p = p
        self.D = D
        self.m = m
        self.E = E
        self.lfac = lfac
        self.h = h
        self.cs = cs
        


class Grid(object):
    """docstring for Grid"""
    def __init__(self, xL=0., xR=1., ncells=100.):
        super(Grid, self).__init__()
        self.ncells = ncells
        self.xL = float(xL)
        self.xR = float(xR)

        self.x   = xL + (np.arange(ncells) + 0.5) * (xR - xL) / float(ncells)
        self.xI  = xL + (np.arange(ncells+1)) * (xR - xL) / float(ncells)

        self.rho= np.zeros(ncells)
        self.v  = np.zeros(ncells)
        self.p  = np.zeros(ncells)

        self.D  = np.zeros(ncells)
        self.m  = np.zeros(ncells)
        self.E  = np.zeros(ncells)

        self.lfac= np.zeros(ncells)
        self.h   = np.zeros(ncells)
        self.cs  = np.zeros(ncells)


    def getState(self, ix):
        return(State(self.rho[ix], 
            self.v[ix], 
            self.p[ix], 
            self.D[ix], 
            self.m[ix], 
            self.E[ix], 
            self.lfac[ix], 
            self.h[ix], 
            self.cs[ix]))


class Setup(object):
    """
    Grid setting and initialisation tools.
    """
    def __init__(self, ncells=100, xL=0., xR=1.):
        super(Setup, self).__init__()
        self.ncells = ncells
        self.xL = xL
        self.xR = xR

        self.grid = Grid(xL, xR, ncells)
  

    def fillGrid():
        pass


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

        self.fillGrid()
        

    def fillGrid(self):
        self.grid.rho[self.grid.x < self.xmid] = self.rhoL
        self.grid.v[self.grid.x < self.xmid] = self.vL
        self.grid.p[self.grid.x < self.xmid] = self.pL

        self.grid.rho[self.grid.x >= self.xmid] = self.rhoR
        self.grid.v[self.grid.x >= self.xmid] = self.vR
        self.grid.p[self.grid.x >= self.xmid] = self.pR

        






