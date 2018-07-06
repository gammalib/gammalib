# ==========================================================================
# This module performs unit tests for the GammaLib optimizer module.
#
# Copyright (C) 2012-2018 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import gammalib
import test_support


# ======================================== #
# Test class for GammaLib optimizer module #
# ======================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib optimizer module.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Setup GOptimizerPars container
    def _setup_pars(self):
        """
        Setup GOptimizerPars container

        Returns
        -------
        pars : `~gammalib.GOptimizerPars`
            GOptimizerPars container
        """
        # Setup GOptimizerPars container
        pars = gammalib.GOptimizerPars()
        par  = gammalib.GOptimizerPar()
        for i in range(10):
            par.name('%s' % i)
            pars.append(par)

        # Return GOptimizerPars container
        return pars

    # Test GOptimizerPars class access operators
    def _test_pars_access(self):
        """
        Test GOptimizerPars class parameter access
        """
        # Setup GOptimizerPars container and parameter
        pars = self._setup_pars()
        par  = gammalib.GOptimizerPar()

        # Perform GOptimizerPars access tests
        test_support._container_access_index(self, pars)

        # Check parameter setting by index from start
        par.name('98')
        pars[3] = par
        self.test_value(pars[3].name(), '98')

        # Check parameter setting by index from end
        par.name('99')
        pars[-2] = par
        self.test_value(pars[-2].name(), '99')

        # Return
        return

    # Test GOptimizerPars class slicing
    def _test_pars_slicing(self):
        """
        Test GOptimizerPars class slicing
        """
        # Setup GOptimizerPars container
        pars = self._setup_pars()

        # Perform slicing tests
        test_support._container_slicing(self, pars)

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Perform pickeling tests of empty classes
        test_support._pickeling(self, gammalib.GOptimizerLM())
        test_support._pickeling(self, gammalib.GOptimizerPar())
        test_support._pickeling(self, gammalib.GOptimizerPars())

        # Setup tests
        pars = gammalib.GOptimizerPars()
        pars.append(gammalib.GOptimizerPar('Test1', 1.0))
        pars.append(gammalib.GOptimizerPar('Test2', 2.0))
        opt = gammalib.GOptimizerLM()
        opt.eps(0.1)

        # Perform pickeling tests of filled classes
        test_support._pickeling(self, gammalib.GOptimizerLM(opt))
        test_support._pickeling(self, gammalib.GOptimizerPar('Test', 2.1))
        test_support._pickeling(self, gammalib.GOptimizerPars(pars))

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('opt')

        # Append tests
        self.append(self._test_pars_access, 'Test GOptimizerPars parameter access')
        self.append(self._test_pars_slicing, 'Test GOptimizerPars slicing')
        self.append(self._test_pickeling, 'Test optimizer class pickeling')

        # Return
        return
