# ==========================================================================
# This module performs unit tests for the GammaLib sky module.
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
import os
import math
import gammalib
import test_support


# =========== #
# Text parser #
# =========== #
def _loadtxt(filename):
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
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib sky module.
    """
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Setup GSkyRegions container
    def _setup_regions(self):
        """
        Setup GSkyRegions container

        Returns
        -------
        regions : `~gammalib.GSkyRegions`
            Sky regions container
        """
        # Setup regions container
        regions = gammalib.GSkyRegions()
        circle  = gammalib.GSkyRegionCircle()
        for i in range(10):
            circle.name('%s' % i)
            regions.append(circle)

        # Return regions container
        return regions

    # Test sky map methods
    def _test_methods(self):
        """
        Test the sky map methods.
        """
        # Setup sky maps
        map = gammalib.GSkyMap('CAR', 'CEL', 83.63, 22.01, -3.7, 2.6, 10, 8, 12)

        # Set sky map shape
        map.shape(3,4)

        # Test map dimensions and shape
        self.test_value(map.npix(), 80, "Check that sky map has 80 pixels")
        self.test_value(map.nx(), 10, "Check that sky map has X dimension of 10")
        self.test_value(map.ny(), 8, "Check that sky map has Y dimension of 8")
        self.test_value(map.nmaps(), 12, "Check that sky map contains 12 maps")
        self.test_value(map.ndim(), 2, "Check that sky map has two dimensions")
        self.test_value(len(map.shape()), 2, "Check that sky map has a shape size of 2")
        self.test_value(map.shape()[0], 3, "Check that sky map has 3 maps in first dimension")
        self.test_value(map.shape()[1], 4, "Check that sky map has 4 maps in second dimension")

        # Set new map shape
        map.shape(2,3,2);

        # Test map dimensions and shape
        self.test_value(map.nmaps(), 12, "Check that sky map contains 12 maps")
        self.test_value(map.ndim(), 3, "Check that sky map has three dimensions")
        self.test_value(len(map.shape()), 3, "Check that sky map has a shape size of 3")
        self.test_value(map.shape()[0], 2, "Check that sky map has 2 maps in first dimension")
        self.test_value(map.shape()[1], 3, "Check that sky map has 3 maps in second dimension")
        self.test_value(map.shape()[2], 2, "Check that sky map has 2 maps in third dimension")

        # Test array method
        array = map.array()
        self.test_value(len(array),     8, 'Check number of lists in list returned by array() method')
        self.test_value(len(array[0]), 10, 'Check number of elements in list returned by array() method')

        # Return
        return

    # Test sky map friend methods
    def _test_friend_methods(self):
        """
        Test the sky map friend methods
        """
        # Setup sky maps
        m = gammalib.GSkyMap('CAR', 'CEL', 83.63, 22.01, -3.7, 2.6, 10, 8, 12)

        # Fill all values of skymaps with numbers > 0
        values = range(1,m.npix()+1)
        for j in range(m.nmaps()):
            for i in range(m.npix()):
                m[i,j] = values[i]

        # Test map values after calling friend methods on one pixel
        ipix = 18
        self.test_value(m.sqrt()[ipix,0], values[ipix]**0.5, 'Check sqrt method')
        self.test_value(m.log()[ipix,0], math.log(values[ipix]), 'Check log method')
        self.test_value(m.log10()[ipix,0], math.log10(values[ipix]), 'Check log10 method')

        # Set the interesting value to a negative value
        m[ipix,0] = -values[ipix]
        self.test_value(m.abs()[ipix,0], values[ipix], 'Check abs method')
        self.test_value(m.sign()[ipix,0], -1, 'Check sign method')

        # Return
        return

    # Test sky map operators
    def _test_operators(self):
        """
        Test the sky map operators.
        """
        # Setup sky maps
        map    = gammalib.GSkyMap('CAR', 'CEL', 83.6331, 22.0145, -3.7, 2.6, 2, 2)
        map[0] = 1.0
        map[1] = 2.0
        map[2] = 3.0
        map[3] = 4.0
        map_b  = map.copy()
        
        # Addition operator
        map += map_b
        self.test_value(map[0], 2.0)
        self.test_value(map[1], 4.0)
        self.test_value(map[2], 6.0)
        self.test_value(map[3], 8.0)

        # Multiplication operator
        map *= map_b
        self.test_value(map[0],  2.0)
        self.test_value(map[1],  8.0)
        self.test_value(map[2], 18.0)
        self.test_value(map[3], 32.0)

        # Subtraction operator
        map -= map_b
        self.test_value(map[0],  1.0)
        self.test_value(map[1],  6.0)
        self.test_value(map[2], 15.0)
        self.test_value(map[3], 28.0)

        # Division operator
        map /= map_b
        self.test_value(map[0], 1.0)
        self.test_value(map[1], 3.0)
        self.test_value(map[2], 5.0)
        self.test_value(map[3], 7.0)

        # Scaling operator
        map *= 2.0
        self.test_value(map[0],  2.0)
        self.test_value(map[1],  6.0)
        self.test_value(map[2], 10.0)
        self.test_value(map[3], 14.0)

        # Division operator
        map /= 2.0
        self.test_value(map[0], 1.0)
        self.test_value(map[1], 3.0)
        self.test_value(map[2], 5.0)
        self.test_value(map[3], 7.0)

        # Access operator (tests also proper iteration)
        sum = 0.0
        for pix in map:
            sum += pix
        self.test_value(sum, 16.0)        

        # Addition operator
        map_res = map + map_b
        self.test_value(map_res[0],  2.0)
        self.test_value(map_res[1],  5.0)
        self.test_value(map_res[2],  8.0)
        self.test_value(map_res[3], 11.0)

        # Subtraction operator
        map_res = map - map_b
        self.test_value(map_res[0], 0.0)
        self.test_value(map_res[1], 1.0)
        self.test_value(map_res[2], 2.0)
        self.test_value(map_res[3], 3.0)

        # Multiplication operator
        map_res = map * map_b
        self.test_value(map_res[0],  1.0)
        self.test_value(map_res[1],  6.0)
        self.test_value(map_res[2], 15.0)
        self.test_value(map_res[3], 28.0)

        # Division operator
        map_res = map / map_b
        self.test_value(map_res[0], 1.0/1.0)
        self.test_value(map_res[1], 3.0/2.0)
        self.test_value(map_res[2], 5.0/3.0)
        self.test_value(map_res[3], 7.0/4.0)

        # Return
        return

    # Generic skymap pixel transformation test
    def _test_skymap_pixels(self, pixels, map):
        """
        Control coordinate and pixel transformation methods.
        """
        # Control inx2dir/dir2inx methods for all pixels
        for i in range(pixels.npix()):
            dir = pixels.inx2dir(i)
            inx = pixels.dir2inx(dir)
            msg = map + ' inx2dir/dir2inx check for pixel ' + str(i)
            err = map + ' GSkyMap trouble with pixel ' + str(i) + ' (' + str(inx) + \
                '), RA=' + str(dir.ra() * 180 / math.pi) + ', Dec=' + str(dir.dec() * 180 / math.pi)
            self.test_assert(i == inx, msg, err)

        # Control SkyDir coordinate transformation for all pixels
        for i in range(pixels.npix()):
            dir     = pixels.inx2dir(i)
            dir_new = gammalib.GSkyDir()
            dir_new.lb(dir.l(), dir.b())
            dra     = abs(dir.ra() - dir_new.ra())
            if (dra >= 5.0):
                dra -= 2.0 * math.pi
            ddec = abs(dir.dec() - dir_new.dec())
            msg = map + ' dir check for pixel ' + str(i)
            err = map + ' GSkyMap trouble with pixel ' + str(i) + ' (' + str(dra) + \
                ',' + str(ddec) + ')'
            self.test_assert(not (dra > 1.0e-9 or ddec > 1.0e-9), msg, err)

        # Return
        return

    # Generic skymap projection test
    def _test_skymap_proj(self, proj):
        """
        Test skymap projection.
        """
        # Set filename
        file = 'test_python_skymap_' + proj + '.fits'

        # Remove test files
        try:
            os.remove(file)
        except:
            pass

        # Create skymap
        pixels = gammalib.GSkyMap(proj, 'CEL', 83.6331, 22.0145, -3.7, 2.6, 5, 5, 20)
        for map in range(pixels.nmaps()):
            for i in range(pixels.npix()):
                pixels[i, map] = i + map * pixels.npix()
        pixels.save(file)

        # Test coordinate and pixel transformations
        self._test_skymap_pixels(pixels, proj)

        # Return
        return

    # Test HEALPix projection
    def _test_skymap_healpix(self):
        """
        Test HEALPix interface.
        """
        # Set filenames
        file1 = 'test_python_skymap_hpx_v1.fits'
        file2 = 'test_python_skymap_hpx_v2.fits'

        # Remove test files
        try:
            os.remove(file1)
            os.remove(file2)
        except:
            pass

        # Create HEALPix skymap
        pixels = gammalib.GSkyMap('GAL', 2, 'RING', 2)
        for i in range(pixels.npix()):
            pixels[i] = i + 1.0
            pixels[i, 1] = i + 1.0 + 1000.0
        pixels.save(file1)

        # Load HEALPix skymap
        pixels = gammalib.GSkyMap(file1)

        # Control coordinate and pixel transformations
        self._test_skymap_pixels(pixels, 'HEALPix')

        # Save HEALPix skymap twice. The second saving should fail.
        try:
            pixels.save(file2, True)
            pixels.save(file2)
        except RuntimeError:
            pass
        else:
            raise RuntimeError('*** TEST ERROR: FITS file overwritten!')
        pixels.save(file2, True)

        # Load again HEALPix skymap
        pixels = gammalib.GSkyMap(file1)

        # Return
        return

    # Test AIT projection
    def _test_skymap_ait(self):
        """
        Test AIT projection.
        """
        # Execute generic test
        self._test_skymap_proj('AIT')

        # Return
        return

    # Test ARC projection
    def _test_skymap_arc(self):
        """
        Test ARC projection.
        """
        # Execute generic test
        self._test_skymap_proj('ARC')

        # Return
        return

    # Test AZP projection
    def _test_skymap_azp(self):
        """
        Test AZP projection.
        """
        # Execute generic test
        self._test_skymap_proj('AZP')

        # Return
        return

    # Test CAR projection
    def _test_skymap_car(self):
        """
        Test CAR projection.
        """
        # Execute generic test
        self._test_skymap_proj('CAR')

        # Return
        return

    # Test GLS projection
    def _test_skymap_gls(self):
        """
        Test GLS projection.
        """
        # Execute generic test
        self._test_skymap_proj('GLS')

        # Return
        return

    # Test MER projection
    def _test_skymap_mer(self):
        """
        Test MER projection.
        """
        # Execute generic test
        self._test_skymap_proj('MER')

        # Return
        return

    # Test MOL projection
    def _test_skymap_mol(self):
        """
        Test MOL projection.
        """
        # Execute generic test
        self._test_skymap_proj('MOL')

        # Return
        return

    # Test SIN projection
    def _test_skymap_sin(self):
        """
        Test SIN projection.
        """
        # Execute generic test
        self._test_skymap_proj('SIN')

        # Return
        return

    # Test SFL projection
    def _test_skymap_sfl(self):
        """
        Test SFL projection.
        """
        # Execute generic test
        self._test_skymap_proj('SFL')

        # Return
        return

    # Test STG projection
    def _test_skymap_stg(self):
        """
        Test STG projection.
        """
        # Execute generic test
        self._test_skymap_proj('STG')

        # Return
        return

    # Test TAN projection
    def _test_skymap_tan(self):
        """
        Test TAN projection.
        """
        # Execute generic test
        self._test_skymap_proj('TAN')

        # Return
        return

    # FK5 J2000 to Galactic coordinate tranformations
    def _test_fk5_to_galactic(self):
        """
        Test precision of FK5 J2000 to Galactic coordinate tranformations
        against pyast (http://dsberry.github.com/starlink/pyast.html)
        for 1000 random positions on the sky.

        Test txt files taken from:
        See https://github.com/astropy/coordinates-benchmark
        """
        # Get test data directory
        datadir = os.environ['TEST_DATA'] + '/'

        # Set parameters
        fk5_filename         = datadir + 'initial_coords.txt'
        galactic_filename    = datadir + 'fk5_to_galactic.txt'
        fk5_coordinates      = _loadtxt(fk5_filename)
        galactic_coordinates = _loadtxt(galactic_filename)
        coordinates          = zip(fk5_coordinates, galactic_coordinates)
        max_dist             = 0

        # Loop over coordinates
        for (ra, dec), (glon, glat) in coordinates:

            # Create point A in FK5 J2000 coordinates
            input_point = gammalib.GSkyDir()
            input_point.radec_deg(ra, dec)

            # Convert point A to Galactic coordinates
            ll, bb = input_point.l_deg(), input_point.b_deg()
            actual_point = gammalib.GSkyDir()
            actual_point.lb_deg(ll, bb)

            # Compute offset to reference result
            reference_point = gammalib.GSkyDir()
            reference_point.lb_deg(glon, glat)

            # Compute distance and remember maximum
            this_dist = actual_point.dist_deg(reference_point)
            max_dist  = max(max_dist, this_dist)

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

    # Test GSkyRegions class access operators
    def _test_regions_access(self):
        """
        Test GSkyRegions class observation access
        """
        # Setup regions container
        regions = self._setup_regions()
        circle  = gammalib.GSkyRegionCircle()

        # Perform regions access tests
        test_support.container_access_index(self, regions)

        # Check regions setting by index from start
        circle.name('98')
        regions[3] = circle
        self.test_value(regions[3].name(), '98')

        # Check observation setting by index from end
        circle.name('99')
        regions[-2] = circle
        self.test_value(regions[-2].name(), '99')

        # Return
        return

    # Test GSkyRegions class slicing
    def _test_regions_slicing(self):
        """
        Test GSkyRegions class slicing
        """
        # Setup regions container
        regions = self._setup_regions()

        # Perform slicing tests
        test_support.container_slicing(self, regions)

        # Return
        return

    # Test class pickeling
    def _test_pickeling(self):
        """
        Test class pickeling
        """
        # Perform pickeling tests of empty classes
        test_support.pickeling(self, gammalib.GHealpix())
        test_support.pickeling(self, gammalib.GHorizDir())
        test_support.pickeling(self, gammalib.GSkyDir())
        test_support.pickeling(self, gammalib.GSkyMap())
        test_support.pickeling(self, gammalib.GSkyPixel())
        test_support.pickeling(self, gammalib.GSkyRegionCircle())
        test_support.pickeling(self, gammalib.GSkyRegionMap())
        test_support.pickeling(self, gammalib.GSkyRegions())
        test_support.pickeling(self, gammalib.GWcsAIT())
        test_support.pickeling(self, gammalib.GWcsARC())
        test_support.pickeling(self, gammalib.GWcsAZP())
        test_support.pickeling(self, gammalib.GWcsCAR())
        test_support.pickeling(self, gammalib.GWcsGLS())
        test_support.pickeling(self, gammalib.GWcsMER())
        test_support.pickeling(self, gammalib.GWcsMOL())
        test_support.pickeling(self, gammalib.GWcsSFL())
        test_support.pickeling(self, gammalib.GWcsSIN())
        test_support.pickeling(self, gammalib.GWcsSTG())
        test_support.pickeling(self, gammalib.GWcsTAN())

        # Setup test
        dir      = gammalib.GSkyDir()
        dir.lb_deg(1.0, 2.0)
        map      = gammalib.GFilename(os.environ['TEST_DATA']+'/test_cube.fits')
        regions  = gammalib.GSkyRegions()
        regions.append(gammalib.GSkyRegionCircle(dir,2.0))
        regions.append(gammalib.GSkyRegionMap(map))
        horizdir = gammalib.GHorizDir()
        horizdir.altaz_deg(2.0, 3.0)

        # Perform pickeling tests of filled classes
        test_support.pickeling(self, gammalib.GHealpix(4))
        test_support.pickeling(self, horizdir)
        test_support.pickeling(self, gammalib.GSkyDir(dir))
        test_support.pickeling(self, gammalib.GSkyMap(map))
        test_support.pickeling(self, gammalib.GSkyPixel(1.0))
        test_support.pickeling(self, gammalib.GSkyPixel(2.0,3.0))
        test_support.pickeling(self, gammalib.GSkyRegionCircle(dir,2.0))
        test_support.pickeling(self, gammalib.GSkyRegionMap(map))
        test_support.pickeling(self, gammalib.GSkyRegions(regions))
        test_support.pickeling(self, gammalib.GWcsAIT('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsARC('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsAZP('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsCAR('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsGLS('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsMER('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsMOL('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsSFL('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsSIN('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsSTG('GAL',1.0,2.0,1.0,1.0,0.1,0.1))
        test_support.pickeling(self, gammalib.GWcsTAN('GAL',1.0,2.0,1.0,1.0,0.1,0.1))

        # Return
        return

    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name('sky')

        # Append tests
        self.append(self._test_methods, 'Test sky map methods')
        self.append(self._test_friend_methods, 'Test sky map friend methods')
        self.append(self._test_operators, 'Test sky map operators')
        self.append(self._test_skymap_healpix, 'Test HEALPix map')
        self.append(self._test_skymap_ait, 'Test AIT projection map')
        self.append(self._test_skymap_arc, 'Test ARC projection map')
        self.append(self._test_skymap_azp, 'Test AZP projection map')
        self.append(self._test_skymap_car, 'Test CAR projection map')
        self.append(self._test_skymap_gls, 'Test GLS projection map')
        self.append(self._test_skymap_mer, 'Test MER projection map')
        self.append(self._test_skymap_mol, 'Test MOL projection map')
        self.append(self._test_skymap_sfl, 'Test SFL projection map')
        self.append(self._test_skymap_sin, 'Test SIN projection map')
        self.append(self._test_skymap_stg, 'Test STG projection map')
        self.append(self._test_skymap_tan, 'Test TAN projection map')
        self.append(self._test_regions_access, 'Test GSkyRegions region access')
        self.append(self._test_regions_slicing, 'Test GSkyRegions slicing')
        self.append(self._test_pickeling, 'Test sky class pickeling')

        # Return
        return
