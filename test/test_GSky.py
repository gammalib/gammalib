# ==========================================================================
# This module performs unit tests for the GammaLib sky module.
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


# =========== #
# Text parser #
# =========== #
def loadtxt(filename):
	"""
	Parse text file with columns of floats.

	We have files with multiple spaces as delimiter,
	which I think csv.reader can't handle.
	"""
	# Initialise result array
	data = []

	# Loop over all rows in file
	for row in [_.split() for _ in open(filename)]:
		row_data = []
		for item in row:
			row_data.append(float(item))
		data.append(row_data)

	# Return result
	return data


# ================================== #
# Test class for GammaLib sky module #
# ================================== #
class Test(GPythonTestSuite):
	"""
	Test class for GammaLib sky module.
	"""
	def __init__(self):
		"""
		Constructor.
		"""
		# Call base class constructor
		GPythonTestSuite.__init__(self)

		# Return
		return

	def set(self):
		"""
		Set all test functions.
		"""
		# Set test name
		self.name("sky")

		# Append tests
		self.append(self.test_skymap_healpix, "Test HEALPix map")
		self.append(self.test_skymap_car, "Test CAR projection map")
		self.append(self.test_skymap_tan, "Test TAN projection map")
		self.append(self.test_skymap_stg, "Test STG projection map")
		self.append(self.test_skymap_azp, "Test AZP projection map")
		self.append(self.test_fk5_to_galactic, "Test FK5 to Galactic coordinate conversion")
		# Return
		return

	# Generic skymap pixel transformation test
	def test_skymap_pixels(self, pixels, map):
		"""
		Control coordinate and pixel transformation methods.
		"""
		# Control pix2dir / dir2pix methods for all pixels
		for i in range(pixels.npix()):
			dir  = pixels.pix2dir(i)
			ipix = pixels.dir2pix(dir)
			msg  = map + " pix2dir/dir2pix check for pixel "+str(i)
			err  = map + " GSkymap trouble with pixel "+str(i)+" ("+str(ipix)+ \
			       "), RA="+str(dir.ra()*180/pi)+", Dec="+str(dir.dec()*180/pi)
			self.test_assert(i == ipix, msg, err)

		# Control SkyDir coordinate transformation for all pixels
		for i in range(pixels.npix()):
			dir     = pixels.pix2dir(i)
			dir_new = GSkyDir()
			dir_new.lb(dir.l(),dir.b())
			dra  = abs(dir.ra()  - dir_new.ra())
			if (dra >= 5.0):
				dra -= 2.0*pi
			ddec = abs(dir.dec() - dir_new.dec())
			msg  = map + " dir check for pixel "+str(i)
			err  = map + " GSkymap trouble with pixel "+str(i)+" ("+str(dra)+ \
			       ","+str(ddec)+")"
			self.test_assert(not (dra > 1.0e-9 or ddec > 1.0e-9), msg, err)

		# Return
		return

	# Generic skymap projection test
	def test_skymap_proj(self, proj):
		"""
		Test skymap projection.
		"""
		# Set filename
		file = "test_python_skymap_"+proj+".fits"

		# Remove test files
		try:
			os.remove(file)
		except:
			pass

		# Create skymap
		pixels = GSkymap(proj, "CEL", 83.6331, 22.0145, -3.7, 2.6, 5, 5, 20)
		for map in range(pixels.nmaps()):
			for i in range(pixels.npix()):
				pixels[i,map] = i+map*pixels.npix()
		pixels.save(file)

		# Test coordinate and pixel transformations
		self.test_skymap_pixels(pixels, proj)

		# Return
		return

	# Test HEALPix projection
	def test_skymap_healpix(self):
		"""
		Test HEALPix interface.
		"""
		# Set filenames
		file1 = "test_python_skymap_hpx_v1.fits"
		file2 = "test_python_skymap_hpx_v2.fits"

		# Remove test files
		try:
			os.remove(file1)
			os.remove(file2)
		except:
			pass

		# Create HEALPix skymap
		pixels = GSkymap("HPX", "GAL", 2, "RING", 2)
		for i in range(pixels.npix()):
			pixels[i]   = i+1.0
			pixels[i,1] = i+1.0 + 1000.0
		pixels.save(file1)

		# Load HEALPix skymap
		pixels = GSkymap(file1)

		# Control coordinate and pixel transformations
		self.test_skymap_pixels(pixels, "HEALPix")

		# Save HEALPix skymap twice. The second saving should fail.
		try:
			pixels.save(file2, True)
			pixels.save(file2)
		except RuntimeError:
			pass
		else:
			raise RuntimeError("*** TEST ERROR: FITS file overwritten!")
		pixels.save(file2, True)

		# Load again HEALPix skymap
		pixels = GSkymap(file1)

		# Return
		return

	# Test CAR projection
	def test_skymap_car(self):
		"""
		Test CAR projection.
		"""
		# Execute generic test
		self.test_skymap_proj("CAR")

		# Return
		return

	# Test TAN projection
	def test_skymap_tan(self):
		"""
		Test TAN projection.
		"""
		# Execute generic test
		self.test_skymap_proj("TAN")

		# Return
		return

	# Test STG projection
	def test_skymap_stg(self):
		"""
		Test STG projection.
		"""
		# Execute generic test
		self.test_skymap_proj("STG")

		# Return
		return

	# Test AZP projection
	def test_skymap_azp(self):
		"""
		Test AZP projection.
		"""
		# Execute generic test
		self.test_skymap_proj("AZP")

		# Return
		return

	# FK5 J2000 to Galactic coordinate tranformations
	def test_fk5_to_galactic(self):
		"""
		Test precision of FK5 J2000 to Galactic coordinate tranformations
		against pyast (http://dsberry.github.com/starlink/pyast.html)
		for 1000 random positions on the sky.

		Test txt files taken from:
		See https://github.com/astropy/coordinates-benchmark
		"""
		# Set parameters
		fk5_filename         = "data/initial_coords.txt"
		galactic_filename    = "data/fk5_to_galactic.txt"
		fk5_coordinates      = loadtxt(fk5_filename)
		galactic_coordinates = loadtxt(galactic_filename)
		coordinates          = zip(fk5_coordinates, galactic_coordinates)
		max_dist             = 0

		# Loop over coordinates
		for (ra, dec), (glon, glat) in coordinates:

			# Create point A in FK5 J2000 coordinates
			input_point = GSkyDir()
			input_point.radec_deg(ra, dec)

			# Convert point A to Galactic coordinates
			ll, bb = input_point.l_deg(), input_point.b_deg()
			actual_point = GSkyDir()
			actual_point.lb_deg(ll, bb)

			# Compute offset to reference result
			reference_point = GSkyDir()
			reference_point.lb_deg(glon, glat)

			# Compute distance and remember maximum
			this_dist = actual_point.dist_deg(reference_point)
			max_dist = max(max_dist, this_dist)

		# Convert max_dist from deg to milli-arcsec
		max_dist *= 1e3 * 3600

		# Results match within 10 milli-arcsec on my machine,
		# I put 20 arcsec here to avoid test failures on machines with
		# slightly different numerical results
		msg = ''
		err = ''
		self.test_assert(max_dist < 20, msg, err)

		# Return
		return
