# ==========================================================================
# This module performs unit tests for the GammaLib model module
#
# Copyright (C) 2012-2016 Juergen Knoedlseder
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

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('model')

        # Append tests
        self.append(self.test_models, 'Test GModels')

        # Return
        return

    # Test function
    def test_models(self):
        """
        Test GModels class
        """
        # Read model container and create a copy of its first model (the
        # copy is needed since the model container will go out of scope
        # later)
        models = gammalib.GModels('data/model_point_plaw.xml')
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
        for model in models:
            nmodels += 1

        # Check that looping was successful
        self.test_value(models.size(), 10, 'Check model container size')
        self.test_value(nmodels, 10, 'Check model iterator')
        
        # Return
        return
