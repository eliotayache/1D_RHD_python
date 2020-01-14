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
    def __init__(self, setup, solver, maxiter = 1000, dumpstep = 10):
        super(Hydro, self).__init__()
        self.setup  = setup
        self.solver = solver
        self.maxiter = maxiter
        self.dumpstep = dumpstep

        self.grid = self.setup.grid

        self.t_arr = []
        self.dt_arr = []
        self.t = 0
        self.dt = 0

        self.it = 0
        self.keep_running = True

    def run(self):
        """
        Runs the simulation and dumps the data in a results directory.
        Needs an initialised Hydro object to run.
        """

        while self.keep_running:
            self.info()
            self.recordTime(*self.solver.timescale(self.grid, self.t))
            self.grid = self.solver.evolve(self.grid)
            self.dump()
            self.checkForEnd()

        return True

    def recordTime(self, t, dt):
        self.t = t
        self.dt = dt
        self.t_arr = np.append(self.t_arr, t)
        self.dt_arr = np.append(self.dt_arr, dt)

    def dump(self):
        pass

    def info(self):
        if self.it%1 == 0.:
            print self.it

    def checkForEnd(self):
        self.it += 1
        if self.it >= self.maxiter:
            self.keep_running = False




class Solver(object):
    """docstring for Solver"""
    def __init__(self, clf = 0.2):
        super(Solver, self).__init__()
        
    def initialize(self, grid):
        pass

    def evolve(self, grid):
        self.computeFluxes(grid)
        grid.D[1:-1] += (self.Fd[:-1] - self.Fd[1:]) * self.dt
        grid.m[1:-1] += (self.Fm[:-1] - self.Fm[1:]) * self.dt
        grid.E[1:-1] += (self.Fe[:-1] - self.Fe[1:]) * self.dt
        grid.cons2prim()
        grid.prim2cons() # for consistency
        return grid


class Solver_HLL(Solver):
    """docstring for Solver_HLL"""
    def __init__(self, cfl = 0.2):
        Solver.__init__(self)
        super(Solver_HLL, self).__init__()
        self.cfl = cfl

    def initialize(self, grid):
        self.lL = np.zeros(grid.ncells)
        self.lR = np.zeros(grid.ncells)

        # left fluxes
        self.Fd = np.zeros(grid.ncells)
        self.Fm = np.zeros(grid.ncells)
        self.Fe = np.zeros(grid.ncells)

    def timescale(self, grid, t):
        #Piecewise constant reconstruction for now
        sigS = grid.cs**2 / (grid.lfac**2 * (1.- grid.cs**2))
        l1   = (grid.v - np.sqrt(sigS * (1. - grid.v**2 + sigS))) / (1.+sigS)
        self.lL = np.minimum(l1[:-1], l1[1:])

        l2   = (grid.v + np.sqrt(sigS * (1. - grid.v**2 + sigS))) / (1.+sigS)
        self.lR = np.minimum(l2[:-1], l2[1:])
        
        dtR = np.min(np.abs(grid.dx[1:] / self.lR))
        dtL = np.min(np.abs(grid.dx[:-1] / self.lL))
        dt  = self.cfl * min(dtR, dtL)
 
        self.dt = dt
        return t + dt, dt

    def computeFluxes(self, grid):
        # Grid Fluxes
        Fgd = grid.D * grid.v
        Fgm = grid.m * grid.v + grid.p
        Fge = grid.m

        # Star region fluxes
        Fsd = (self.lR*Fgd[:-1] - self.lL*Fgd[1:] \
            + self.lR*self.lL*(grid.D[1:]-grid.D[:-1])) / (self.lR - self.lL)
        Fsm = (self.lR*Fgm[:-1] - self.lL*Fgm[1:] \
            + self.lR*self.lL*(grid.m[1:]-grid.m[:-1])) / (self.lR - self.lL)
        Fse = (self.lR*Fge[:-1] - self.lL*Fge[1:] \
            + self.lR*self.lL*(grid.E[1:]-grid.E[:-1])) / (self.lR - self.lL)

        condlist = [0 < self.lL, np.logical_and(self.lL <= 0, 0 < self.lR), self.lR <= 0.]
        dchoice  = [Fgd[:-1], Fsd, Fgd[1:]]
        mchoice  = [Fgm[:-1], Fsm, Fgm[1:]]
        echoice  = [Fge[:-1], Fse, Fge[1:]]

        self.Fd = np.select(condlist, dchoice)
        self.Fm = np.select(condlist, mchoice)
        self.Fe = np.select(condlist, echoice)

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

            self.p = optimize.brentq(self.f, p_lo, p_hi, args=(params), rtol=1.e-14)

        self.lfac = 1. \
            / np.sqrt(1 - (self.m * self.m) / (self.E + self.p)**2)
        self.rho = self.D / self.lfac
        if self.m == 0:
            self.v == 0
        if self.m > 0:
            self.v = np.sqrt(1. - 1. / self.lfac**2);
        if self.m < 0:
            self.v = - np.sqrt(1. - 1. / self.lfac**2);


    # def state2flux(self):





class Grid(object):
    """docstring for Grid"""
    def __init__(self, xL=0., xR=1., ncells=100):
        super(Grid, self).__init__()
        self.ncells = ncells
        self.xL = float(xL)
        self.xR = float(xR)

        self.x   = xL + (np.arange(ncells) + 0.5) * (xR - xL) / float(ncells)
        self.xI  = xL + (np.arange(ncells+1)) * (xR - xL) / float(ncells)
        self.dx  = self.xI[1:] - self.xI[:-1]

        self.rho= np.zeros(ncells)
        self.v  = np.zeros(ncells)
        self.p  = np.zeros(ncells)

        self.D  = np.zeros(ncells)
        self.m  = np.zeros(ncells)
        self.E  = np.zeros(ncells)

        self.lfac= np.zeros(ncells)
        self.h   = np.zeros(ncells)
        self.cs  = np.zeros(ncells)


    def readState(self, ix):
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

    def prim2cons(self):
        for ix in range(self.ncells):
            S = self.readState(ix)
            S.prim2cons()
            self.writeState(S, ix)

    def cons2prim(self):
        for ix in range(self.ncells):
            S = self.readState(ix)
            S.cons2prim()
            self.writeState(S, ix)





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

        self.grid.prim2cons()

        


if __name__ == '__main__':
    setup = Setup_ST(rhoL=1.,vL=.1,pL=1.,rhoR=0.1,vR=0.,pR=0.1)
    solver = Solver('HLL', cfl=0.2)
    hydro = Hydro(setup, solver)




