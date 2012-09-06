# ==========================================================================
# This module performs unit tests for the GammaLib linalg module.
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


# ===================================== #
# Test class for GammaLib linalg module #
# ===================================== #
class Test(GPythonTestSuite):
	"""
	Test class for GammaLib linalg module.
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
		self.name("linalg")

		# Append tests
		self.append(self.test_matrix, "Test GMatrix")
		self.append(self.test_sym_matrix, "Test GSymMatrix")
		self.append(self.test_sparse_matrix, "Test GSparseMatrix")

		# Return
		return


	# Test GMatrix class
	def test_matrix(self):
		"""
		Test GMatrix class.
		"""
		# Set matrix size
		nrows = 3
		ncols = 5
		
		# Allocate matrix
		m = GMatrix(nrows,ncols)
		
		# Fill elements
		for i in range(nrows):
			for j in range(ncols):
				m[i,j] = i+j
		
		# Check elements
		for i in range(nrows):
			for j in range(ncols):
				ref = i+j
				self.test_value(m[i,j],ref,0.0,"Test matrix element access")
		
		# Return
		return


	# Test GSymMatrix class
	def test_sym_matrix(self):
		"""
		Test GSymMatrix class.
		"""
		# Set matrix size
		nrows = 3
		ncols = 3
		
		# Allocate matrix
		m = GSymMatrix(nrows,ncols)
		
		# Fill elements
		for i in range(nrows):
			for j in range(ncols):
				m[i,j] = i+j
		
		# Check elements
		for i in range(nrows):
			for j in range(ncols):
				ref = i+j
				self.test_value(m[i,j],ref,0.0,"Test matrix element access")
		
		# Return
		return


	# Test GSparseMatrix class
	def test_sparse_matrix(self):
		"""
		Test GSparseMatrix class.
		"""
		# Set matrix size
		nrows = 5
		ncols = 3
		
		# Allocate matrix
		m = GSparseMatrix(nrows,ncols)
		
		# Fill elements
		for i in range(nrows):
			for j in range(ncols):
				m[i,j] = i+j
		
		# Check elements
		for i in range(nrows):
			for j in range(ncols):
				ref = i+j
				self.test_value(m[i,j],ref,0.0,"Test matrix element access")
		
		# Return
		return
