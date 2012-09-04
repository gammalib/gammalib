#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the GammaLib Python bindings.
#
# Copyright (C) 2012 Juergen Knoedlseder
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
from gammalib import *
from math import *
import sys
import test_GApp
import test_GLinalg
import test_GModel
import test_GNumerics
import test_GObs
import test_GOpt
import test_GFits
import test_GSky
import test_GSupport
import test_GTest
import test_GXml


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
	"""
	Perform unit testing for Python interface.
	"""
	# Allocate test suites
	suites = GTestSuites("Python interface testing")

	# Allocate test suites and append them to the container
	suite1 = test_GApp.Test()
	suite2 = test_GFits.Test()
	suite3 = test_GLinalg.Test()
	suite4 = test_GModel.Test()
	suite5 = test_GNumerics.Test()
	suite6 = test_GObs.Test()
	suite7 = test_GOpt.Test()
	suite8 = test_GSky.Test()
	suite9 = test_GSupport.Test()
	suite10 = test_GTest.Test()
	suite11 = test_GXml.Test()

	# Setup unit tests
	suite1.set()
	suite2.set()
	suite3.set()
	suite4.set()
	suite5.set()
	suite6.set()
	suite7.set()
	suite8.set()
	suite9.set()
	suite10.set()
	suite11.set()

	# Append tests to container
	suites.append(suite1)
	suites.append(suite2)
	suites.append(suite3)
	suites.append(suite4)
	suites.append(suite5)
	suites.append(suite6)
	suites.append(suite7)
	suites.append(suite8)
	suites.append(suite9)
	suites.append(suite10)
	suites.append(suite11)

	# Run test suite
	success = suites.run()
	#print suites

	# Save test results
	suites.save("reports/GPython.xml")

	# Set return code
	if success:
		rc = 0
	else:
		rc = 1

	# Print a linefeed
	sys.stdout.write("\n")

	# Exit with return code
	sys.exit(rc)
	
