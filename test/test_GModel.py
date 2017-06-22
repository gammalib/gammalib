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
import os
import gammalib


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

        # Set test data directory
        self._datadir = os.environ['TEST_DATA']

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
        self.append(self.test_models_access, 'Test GModels class model access')
        self.append(self.test_models_slicing, 'Test GModels class slicing')

        # Return
        return

    # Test function
    def test_models_access(self):
        """
        Test GModels class model access
        """
        # Set model filename
        filename = self._datadir + '/model_point_plaw.xml'

        # Read model container and create a copy of its first model (the
        # copy is needed since the model container will go out of scope
        # later)
        models = gammalib.GModels(filename)
        model  = models[0].copy()
        name   = model.name()

        # Create an empty model container
        models = gammalib.GModels()

        # Append 10 identical models to model container
        for i in range(10):
            model.name(name+'_'+str(i))
            models.append(model)
        
        # Loop over all models using the container iterator and count the
        # number of iterations
        nmodels = 0
        for m in models:
            self.test_value(m.name(), name+'_'+str(nmodels))
            nmodels += 1

        # Check that looping was successful
        self.test_value(models.size(), 10, 'Check model container size')
        self.test_value(nmodels, 10, 'Check model iterator')

        # Check model access from end
        self.test_value(models[-1].name(), name+'_9')

        # Check model access by name
        self.test_value(models[name+'_7'].name(), name+'_7')

        # Check model setting by index
        model.name('Test index')
        models[0] = model
        self.test_value(models[0].name(), 'Test index')

        # Check model setting by name
        model.name('Test name')
        models[name+'_3'] = model
        self.test_value(models['Test name'].name(), 'Test name')

        # Return
        return


    # Test GModels class slicing
    def test_models_slicing(self):
        """
        Test GModels class slicing
        """
        # Setup model container
        models = gammalib.GModels()
        model  = gammalib.GModelSky()
        for i in range(10):
            model.name('%s' % i)
            models.append(model)

        # Test models[start:end]
        self.test_value(len(models[3:5]), 2)
        self.test_value(models[3:5][0].name(), '3')
        self.test_value(models[3:5][1].name(), '4')

        # Test models[start:]
        self.test_value(len(models[7:]), 3)
        self.test_value(models[7:][0].name(), '7')
        self.test_value(models[7:][1].name(), '8')
        self.test_value(models[7:][2].name(), '9')

        # Test models[:end]
        self.test_value(len(models[:2]), 2)
        self.test_value(models[:2][0].name(), '0')
        self.test_value(models[:2][1].name(), '1')

        # Test models[:]
        self.test_value(len(models[:]), 10)
        for i in range(10):
            self.test_value(models[:][i].name(), '%s' % i)

        # Test models[start:end:step]
        self.test_value(len(models[3:7:2]), 2)
        self.test_value(models[3:7:2][0].name(), '3')
        self.test_value(models[3:7:2][1].name(), '5')

        # Test models[start:end:step]
        self.test_value(len(models[6:3:-2]), 2)
        self.test_value(models[6:3:-2][0].name(), '6')
        self.test_value(models[6:3:-2][1].name(), '4')

        # Test models[-start:]
        self.test_value(len(models[-2:]), 2)
        self.test_value(models[-2:][0].name(), '8')
        self.test_value(models[-2:][1].name(), '9')

        # Test models[:-end]
        self.test_value(len(models[:-7]), 3)
        self.test_value(models[:-7][0].name(), '0')
        self.test_value(models[:-7][1].name(), '1')
        self.test_value(models[:-7][2].name(), '2')

        # Return
        return
