# -*- coding: utf-8 -*-
from . import helpers

import numpy as np
from scipy import optimize

GAMMA_ = 4./3.

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


class Solver(object):
    """docstring for Solver"""
    def __init__(self, clf = 0.2):
        super(Solver, self).__init__()
        
    def initialize(grid):
        pass


class Solver_HLL(Solver):
    """docstring for Solver_HLL"""
    def __init__(self, cfl = 0.2):
        Solver.__init__(self)
        super(Solver_HLL, self).__init__()
        self.cfl = cfl

    def initialize(grid):
        lL = np.zeros(grid.ncells)
        lR = np.zeros(grid.ncells)




class State(object):
    """docstring for State"""
    def __init__(self, rho=0, v=0, p=1., D=0, m=0, E=0, lfac=1., u=0, h=0, cs=0):
        super(State, self).__init__()
        self.rho = rho
        self.v = v
        self.p = p
        self.D = D
        self.m = m
        self.E = E
        self.lfac = lfac
        self.u = u
        self.h = h
        self.cs = cs

    @classmethod
    def fromPrim(cls, rho, v, p):
        S = cls(rho=rho, v=v, p=p)
        S.prim2cons()
        return S

    @classmethod
    def fromCons(cls, D, m, E):
        S = cls(D=D, m=m, E=E)
        S.cons2prim()
        return S

    def parse(self):
        return(self.rho, 
            self.v, 
            self.p, 
            self.D, 
            self.m, 
            self.E, 
            self.lfac, 
            self.u, 
            self.h, 
            self.cs)

    def prim2cons(self):
        self.lfac = np.sqrt(1. / (1. - self.v*self.v))
        self.u = self.lfac * self.v
        self.h = 1. + self.p * GAMMA_ / (GAMMA_ - 1.) / self.rho
        self.cs = np.sqrt(GAMMA_ * self.p / (self.rho*self.h))

        self.D = self.lfac * self.rho
        self.m = self.D * self.h * self.lfac * self.v
        self.E = self.D * self.h * self.lfac - self.p

    def prim2aux(self):
        self.lfac = np.sqrt(1. / (1. - self.v*self.v))
        self.u = self.lfac * self.v
        self.h = 1. + self.p * GAMMA_ / (GAMMA_ - 1.) / self.rho
        self.cs = np.sqrt(GAMMA_ * self.p / (self.rho*self.h))
    

    def f(self, p, D, m, E):
        lfac = 1. / np.sqrt(1. - (m * m) / ((E + p) * (E + p)))
        return (E + p - D * lfac - (GAMMA_ * p * lfac * lfac) / (GAMMA_ - 1.))


    def cons2prim(self):

        p = self.p
        params = self.D, self.m, self.E

        f_init = self.f(p, *params);

        if abs(f_init) > 1.e-15:
            p_lo = max(0., (1. + 1e-13) * abs(self.m) - self.E);
            if f_init < 0:
                p_hi = 1. * p
            else:
                i = 0                
                p_hi = p
                while self.f(p_hi, *params) > 0:
                    i += 1
                    p_hi *= 10
                    if i == 14:
                        exit(20)

            self.p = optimize.brentq(self.f, p_lo, p_hi, args=(params))

        self.lfac = 1. \
            / np.sqrt(1 - (self.m * self.m) / (self.E + self.p)**2)
        self.rho = self.D / self.lfac
        if self.m == 0:
            self.v == 0
        if self.m > 0:
            self.v = np.sqrt(1. - 1. / self.lfac**2);
        if self.m < 0:
            self.v = - np.sqrt(1. - 1. / self.lfac**2);


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


    def writeState(self, state, ix):
        self.rho[ix] = state.rho
        self.v[ix] = state.v
        self.p[ix] = state.p
        self.D[ix] = state.D
        self.m[ix] = state.m
        self.E[ix] = state.E
        self.lfac[ix] = state.lfac
        self.h[ix] = state.h
        self.cs[ix] = state.cs


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

        


if __name__ == '__main__':
    setup = Setup_ST(rhoL=1.,vL=.1,pL=1.,rhoR=0.1,vR=0.,pR=0.1)
    solver = Solver('HLL', cfl=0.2)
    hydro = Hydro(setup, solver)




