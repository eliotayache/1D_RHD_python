# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2020-01-14 08:01:05
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2020-01-14 08:46:29

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    # Adding parent to Pythonpath

from RHD import *




if __name__ == '__main__':

    setup = Setup_ST(rhoL=1.,vL=.0,pL=1.,rhoR=0.1,vR=0.,pR=0.1)
    solver = Solver_HLL(cfl=0.2)        
    hydro = Hydro(setup=setup, solver=solver, maxiter=5000)
    hydro.run()

