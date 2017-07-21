# ==========================================================================
# This module performs unit tests for the GammaLib model module
#
# Copyright (C) 2012-2017 Juergen Knoedlseder
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


# ==================================== #
# Test class for GammaLib model module #
# ==================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib model module
    """
    # Constructor
    def __init__(self):
        """
        Constructor
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Setup GModels container
    def _setup_models(self):
        """
        Setup GModels container

        Returns
        -------
        models : `~gammalib.GModels`
            Models container
        """
        # Setup model container
        models = gammalib.GModels()
        model  = gammalib.GModelSky()
        for i in range(10):
            model.name('%s' % i)
            models.append(model)

        # Return model container
        return models

    # Test GModels class access operators
    def _test_models_access(self):
        """
        Test GModels class model access
        """
        # Setup model container and sky model
        models = self._setup_models()
        model  = gammalib.GModelSky()

        # Perform model access tests
        test_support._container_access_index(self, models)

        # Check model getting by name
        self.test_value(models['7'].name(), '7')

        # Check model setting by name
        model.name('Test name')
        models['3'] = model
        self.test_value(models['Test name'].name(), 'Test name')

        # Check model setting by index from start
        model.name('Test index from start')
        models[1] = model
        self.test_value(models[1].name(), 'Test index from start')

        # Check model setting by index from end
        model.name('Test index from end')
        models[-2] = model
        self.test_value(models[-2].name(), 'Test index from end')

        # Return
        return

    # Test GModels class slicing
    def _test_models_slicing(self):
        """
        Test GModels class slicing
        """
        # Setup models container
        models = self._setup_models()

        # Perform slicing tests
        test_support._container_slicing(self, models)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('model')

        # Append tests
        self.append(self._test_models_access, 'Test GModels model access')
        self.append(self._test_models_slicing, 'Test GModels slicing')

        # Return
        return
