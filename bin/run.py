# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2020-01-14 08:01:05
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2020-02-27 19:10:46

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    # Adding parent to Pythonpath

import matplotlib.pyplot as plt

from RHD import *


if __name__ == '__main__':

    setup = Setup_2DST(rhoL=1.,vL=0.9,pL=1.,rhoR=1.,vR=-0.9,pR=1.)