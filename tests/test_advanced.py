# -*- coding: utf-8 -*-

from .context import RHD
import numpy as np

import unittest


class AdvancedTestSuite(unittest.TestCase):
    """Advanced test cases."""

    def test_setup(self):
        setup = RHD.Setup_ST(1,1,1,1,1,1)
        print setup.ncells
        assert setup.ncells == 100

    def test_run(self):
        setup = np.arange(100)
        solver = 'HLL'
        hydro = RHD.Hydro(setup=setup, solver=solver)
        self.assertTrue(hydro.run())


if __name__ == '__main__':
    unittest.main()
