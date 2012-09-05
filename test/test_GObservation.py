# ==========================================================================
# This module performs unit tests for the GammaLib observation module.
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
import os


# ========================================== #
# Test class for GammaLib observation module #
# ========================================== #
class Test(GPythonTestSuite):
	"""
	Test class for GammaLib observation module.
	"""
	# Constructor
	def __init__(self):
		"""
		Constructor.
		"""
		# Call base class constructor
		GPythonTestSuite.__init__(self)

		# Return
		return


	# Set test functions
	def set(self):
		"""
		Set all test functions.
		"""
		# Set test name
		self.name("obs")

		# Append tests
		self.append(self.test, "Observation module dummy test")

		# Return
		return


	# Test function
	def test(self):
		"""
		Test function.
		"""
		# Return
		return
