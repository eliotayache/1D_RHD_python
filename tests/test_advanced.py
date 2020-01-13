# -*- coding: utf-8 -*-

from .context import RHD
import numpy as np

import unittest


class AdvancedTestSuite(unittest.TestCase):
    """Advanced test cases."""

    def test_grid(self):
        grid = RHD.Grid(0., 1., 100)
        assert np.all(grid.x == (np.arange(100) + 0.5) / float(100))
        assert np.all(grid.xI == (np.arange(101)) / float(100))

    def test_setup(self):
        setup = RHD.Setup()
        assert np.all(setup.grid.x == (np.arange(100) + 0.5) / float(100))

    def test_setup_ST(self):
        setup = RHD.Setup_ST(rhoL=1.,vL=.1,pL=1.,rhoR=0.1,vR=0.,pR=0.1)
        tar_v = np.zeros(100)
        tar_v[:50] = 0.1
        tar_v[50:] = 0.
        assert np.all(setup.grid.v == tar_v)

    def test_getState(self):
        setup = RHD.Setup_ST(rhoL=1.,vL=.1,pL=1.,rhoR=0.1,vR=0.,pR=0.1)
        state = setup.grid.getState(20)
        assert state.rho == 1.
        state2 = RHD.State.fromPrim(1.,0.,1.)
        assert state2.rho == 1.

    def test_writeState(self):
        setup = RHD.Setup_ST(rhoL=1.,vL=.1,pL=1.,rhoR=0.1,vR=0.,pR=0.1)
        state = setup.grid.getState(20)
        state.rho = 2.
        setup.grid.writeState(state,21)
        assert setup.grid.rho[21] == 2.

    def test_prim2cons(self):
        S = RHD.State.fromPrim(1.,0.,1.)
        S.prim2aux()
        S.prim2cons()
        print S.parse()
        

    def test_run(self):
        setup = np.arange(100)
        solver = 'HLL'
        hydro = RHD.Hydro(setup=setup, solver=solver)
        self.assertTrue(hydro.run())


if __name__ == '__main__':
    unittest.main()
