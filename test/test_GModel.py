# ==========================================================================
# This module performs unit tests for the GammaLib model module.
#
# Copyright (C) 2012-2015 Juergen Knoedlseder
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


# ==================================== #
# Test class for GammaLib model module #
# ==================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib model module.
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

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("model")

        # Append tests
        self.append(self.test_iterator, "Test Model iterator")

        # Return
        return

    # Test function
    def test_iterator(self):
        """
        Test function.
        """
        
        # Initialise model container
        models = gammalib.GModels()
        
        # Read model container
        models_test = gammalib.GModels('data/model_point_plaw.xml')
        model = models_test[0]
        name = model.name()
        
        for i in range(10):
            model.name(name+'_'+str(i))
            models.append(model)
        
        size = models.size()
        nmodels = 0
        for model in models:
            nmodels += 1
            
        self.test_assert(size == 10, 'Test model container size')
        self.test_assert(nmodels == 10, 'Test model iterator')
        
        # Return
        return
