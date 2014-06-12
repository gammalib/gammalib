/***************************************************************************
 *                  test_GSky.cpp - Test sky module                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file test_GSky.cpp
 * @brief Testing of sky module.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception
#include <stdlib.h>
#include "test_GSky.hpp"
#include "GTools.hpp"


/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_WCS_FORTH_BACK_PIXEL_DEBUG
#define G_WCS_COPY_DEBUG


/***********************************************************************//**
 * @brief Set parameters and tests
 ***************************************************************************/    
void TestGSky::set(void){
  
    // Test name
    name("GSky");
  
    // Append tests
    append(static_cast<pfunction>(&TestGSky::test_GWcs),"Test GWcs");
    append(static_cast<pfunction>(&TestGSky::test_GSkyPixel),"Test GSkyPixel");
    append(static_cast<pfunction>(&TestGSky::test_GSkymap_healpix_construct),"Test Healpix GSkymap constructors");
    append(static_cast<pfunction>(&TestGSky::test_GSkymap_healpix_io),"Test Healpix GSkymap I/O");
    append(static_cast<pfunction>(&TestGSky::test_GSkymap_wcs_construct),"Test WCS GSkymap constructors");
    append(static_cast<pfunction>(&TestGSky::test_GSkymap_wcs_io),"Test WCS GSkymap I/O");
    append(static_cast<pfunction>(&TestGSky::test_GSkymap),"Test GSkymap");
    append(static_cast<pfunction>(&TestGSky::test_GSkyRegions_io),"Test GSkyRegions");
    append(static_cast<pfunction>(&TestGSky::test_GSkyRegionCircle_construct),"Test GSkyRegionCircle constructors");
    append(static_cast<pfunction>(&TestGSky::test_GSkyRegionCircle_logic),"Test GSkyRegionCircle logic");
    append(static_cast<pfunction>(&TestGSky::test_GHorizDir),"Test GHorizDir");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGSky* TestGSky::clone(void) const
{
    // Clone test suite
    return new TestGSky(*this);
}


/***********************************************************************//**
 * @brief Test consistency of forward and background transformations
 *
 * @param[in] wcs WCS object
 * @param[in] nx Number of points in X
 * @param[in] ny Number of points in Y
 * @param[in] crpix1 Reference pixel in X
 * @param[in] crpix2 Reference pixel in Y
 ***************************************************************************/
double TestGSky::wcs_forth_back_pixel(GWcs* wcs, int nx, int ny, double& crpix1, double& crpix2)
{
    // Initialise maximal distance
    double dist_max = 0.0;
    
    // Loop over x pixels
    for (int ix = -nx; ix <= nx; ++ix) {
        
        // Set x value
        double x = double(ix) + crpix1;

        // Skip pixels outside valid range (they are not expected to
        // transform bijectively)
        if (x < 0 || x >= nx) {
            continue;
        }

        // Loop over y pixels
        for (int iy = -ny; iy <= ny; ++iy) {

            // Set y value
            double    y = double(iy) + crpix2;

            // Skip pixels outside valid range (they are not expected to
            // transform bijectively)
            if (y < 0 || y >= ny) {
                continue;
            }

            // Set sky pixel
            GSkyPixel pix_in(x,y);

            // Forth: Transform to world
            GSkyDir dir = wcs->pix2dir(pix_in);

            // Back: Transform to pixel
            GSkyPixel pix_out = wcs->dir2pix(dir);

            // Compute distance
            double dx      = pix_in.x()-pix_out.x();
            double dy      = pix_in.y()-pix_out.y();
            double dist    = sqrt(dx*dx+dy*dy);

            // Debug option: Dump discrepancies
            #if defined(G_WCS_FORTH_BACK_PIXEL_DEBUG)
            if (dist > 0.001) {
                std::cout << std::endl;
                std::cout << "dist=" << dist;
                std::cout << " dx=" << dx;
                std::cout << " dy=" << dy;
                std::cout << " pix_in=" << pix_in;
                std::cout << " pix_out=" << pix_out;
                std::cout << " (x,y)=(" << x << "," << y << ")";
                std::cout << " dir=" << dir;
            }
            #endif

            // Store maximum distance
            if (dist > dist_max) {
                dist_max = dist;
            }

        } // endfor: y pixels
    } // endfor: x pixels

    // Return
    return dist_max;
}


/***********************************************************************//**
 * @brief Test consistency of copy
 *
 * @param[in] wcs WCS object
 * @param[in] nx Number of points in X
 * @param[in] ny Number of points in Y
 * @param[in] crpix1 Reference pixel in X
 * @param[in] crpix2 Reference pixel in Y
 ***************************************************************************/
double TestGSky::wcs_copy(GWcs* wcs, int nx, int ny, double& crpix1, double& crpix2)
{
    // Make a copy using clone
    GWcs* cpy = wcs->clone();

    // Initialise maximal angle and distance
    double angle_max = 0.0;
    double dist_max  = 0.0;

    // Loop over x pixels
    for (int ix = -nx; ix <= nx; ++ix) {

        // Set x value
        double x = double(ix) + crpix1;

        // Skip pixels outside valid range (they are not expected to
        // transform bijectively)
        if (x < 0 || x >= nx) {
            continue;
        }

        // Loop over y pixels
        for (int iy = -ny; iy <= ny; ++iy) {

            // Set y value
            double y = double(iy) + crpix2;

            // Skip pixels outside valid range (they are not expected to
            // transform bijectively)
            if (y < 0 || y >= ny) {
                continue;
            }

            // Set sky pixel
            GSkyPixel pix_in(x,y);

            // Transform to world
            GSkyDir dir1 = wcs->pix2dir(pix_in);
            GSkyDir dir2 = cpy->pix2dir(pix_in);
            double angle = dir1.dist_deg(dir2);
            if (angle > angle_max) {
                angle_max = angle;
            }

            // Debug option: Dump discrepancies
            #if defined(G_WCS_COPY_DEBUG)
            if (angle > 0.001) {
                std::cout << std::endl;
                std::cout << "angle=" << angle;
                std::cout << " dir1=" << dir1;
                std::cout << " dir2=" << dir2;
            }
            #endif

            // Transform to pixel
            GSkyPixel pix_out1 = wcs->dir2pix(dir1);
            GSkyPixel pix_out2 = cpy->dir2pix(dir2);

            // Compute distance
            double dx   = pix_out1.x()-pix_out2.x();
            double dy   = pix_out1.y()-pix_out2.y();
            double dist = std::sqrt(dx*dx+dy*dy);
            if (dist > dist_max) {
                dist_max = dist;
            }

            // Debug option: Dump discrepancies
            #if defined(G_WCS_COPY_DEBUG)
            if (dist > 0.001) {
                std::cout << std::endl;
                std::cout << "dist=" << dist;
                std::cout << " dx=" << dx;
                std::cout << " dy=" << dy;
                std::cout << " pix_out1=" << pix_out1;
                std::cout << " pix_out2=" << pix_out2;
                std::cout << " dir1=" << dir1;
                std::cout << " dir2=" << dir2;
            }
            #endif

        } // endfor: y pixels
    } // endfor: x pixels

    // Delete copy
    delete cpy;

    // Return
    return ((dist_max > angle_max) ? dist_max : angle_max);
}


/***********************************************************************//**
 * @brief Test GSkyPixel class
 *
 * Test GSkyPixel class.
 ***************************************************************************/
void TestGSky::test_GSkyPixel(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GSkyPixel pixel;
        test_value(pixel.size(), 0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 1D constructor
    test_try("Test 1D constructor (int version)");
    try {
        GSkyPixel pixel(41);
        test_assert(pixel.is_1D(), "Pixel is not 1D");
        test_assert(!pixel.is_2D(), "Pixel is 2D but it should be 1D");
        test_value(pixel.size(), 1);
        test_value(pixel.index(), 41.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 1D constructor
    test_try("Test 1D constructor (double version)");
    try {
        GSkyPixel pixel(41.001);
        test_assert(pixel.is_1D(), "Pixel is not 1D");
        test_assert(!pixel.is_2D(), "Pixel is 2D but it should be 1D");
        test_value(pixel.size(), 1);
        test_value(pixel.index(), 41.001);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 2D constructor
    test_try("Test 2D constructor (int version)");
    try {
        GSkyPixel pixel(41,14);
        test_assert(!pixel.is_1D(), "Pixel is 1D but it should be 2D");
        test_assert(pixel.is_2D(), "Pixel is not 2D");
        test_value(pixel.size(), 2);
        test_value(pixel.x(), 41.0);
        test_value(pixel.y(), 14.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 2D constructor
    test_try("Test 2D constructor (double version)");
    try {
        GSkyPixel pixel(41.001,14.003);
        test_assert(!pixel.is_1D(), "Pixel is 1D but it should be 2D");
        test_assert(pixel.is_2D(), "Pixel is not 2D");
        test_value(pixel.size(), 2);
        test_value(pixel.x(), 41.001);
        test_value(pixel.y(), 14.003);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 1D indexing
    test_try("Test 1D indexing");
    try {
        GSkyPixel pixel;
        test_value(pixel.size(), 0);
        pixel = 41;
        int index = pixel;
        test_value(index, 41);
        pixel = 41.001;
        index = pixel;
        test_value(index, 41);
        double dindex = pixel;
        test_value(dindex, 41.001);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test cloning
    test_try("Test cloning");
    try {
        GSkyPixel pixel;
        test_value(pixel.size(), 0);
        GSkyPixel* clone = pixel.clone();
        test_value(clone->size(), 0);
        delete clone;
        pixel = 41;
        clone = pixel.clone();
        test_value(clone->size(), 1);
        int index = pixel;
        test_value(index, 41);
        delete clone;
        pixel.x(41.001);
        pixel.y(14.003);
        clone = pixel.clone();
        test_value(clone->size(), 2);
        test_value(clone->x(), 41.001);
        test_value(clone->y(), 14.003);
        delete clone;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GWcs projections
 *
 * This method tests all WCS projections that have been registered.
 ***************************************************************************/
void TestGSky::test_GWcs(void)
{
    // Allocate parameters for testing
    double crval1 = 83.63;
    double crval2 = 22.01;
    double crpix1 = 100.5;
    double crpix2 = 100.5;
    double cdelt1 =   0.5;
    double cdelt2 =   0.5;
    int    nx     =   200;
    int    ny     =   200;

    // Loop over all non-HealPix projections
    GWcsRegistry registry;
    for (int i = 0; i < registry.size(); ++i) {

        // Set header
        std::string test_class = "GWcs" + registry.code(i);

        // Perform tests
        test_try("Test "+test_class);
        try {

            // Allocate projection CEL and GAL
            GWcs *cel = dynamic_cast<GWcs*>(registry.alloc(registry.code(i)));
            GWcs *gal = dynamic_cast<GWcs*>(registry.alloc(registry.code(i)));

            // Throw an error if allocation failed
            if (cel == NULL || gal == NULL) {
                std::cout << std::endl 
                        << "TEST ERROR: Could not allocate "+test_class+"."
                        << std::endl;
                throw exception_error("Could not allocate"+test_class);
            }

            // Set projections
            test_try("Set projections");
            try {
                cel->set("CEL", crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
                gal->set("GAL", crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
                test_try_success();
            }
            catch (std::exception &e) {
                test_try_failure(e);
            }

            // Test forth-back CEL conversion
            test_try("Test forth-back CEL conversion");
            try {
                double tol = 0.0;
                if ((tol = wcs_forth_back_pixel(cel, nx, ny, crpix1, crpix2)) > 1.0e-10) {
                    throw exception_failure("CEL forth-back transformation tolerance 1.0e-10 exceeded: "+gammalib::str(tol));
                }
                test_try_success(); 
            }
            catch (std::exception &e) {
                test_try_failure(e);
            }

            // Test forth-back GAL conversion
            test_try("Test forth-back GAL conversion");
            try {
                double tol = 0.0;
                if ((tol = wcs_forth_back_pixel(gal, nx, ny, crpix1, crpix2)) > 1.0e-10) {
                    throw exception_failure("GAL forth-back transformation tolerance 1.0e-10 exceeded: "+gammalib::str(tol));
                }
                test_try_success();
            }
            catch (std::exception &e) {
                test_try_failure(e);
            }

            // Test CEL copy
            test_try("Test CEL copy");
            try {
                double tol = 0.0;
                if ((tol = wcs_copy(cel, nx, ny, crpix1, crpix2)) > 1.0e-4) {
                    throw exception_failure("CEL copy tolerance 1.0e-4 exceeded: "+gammalib::str(tol));
                }
                test_try_success();
            }
            catch (std::exception &e) {
                test_try_failure(e);
            }

            // Test GAL copy
            test_try("Test GAL copy");
            try {
                double tol = 0.0;
                if ((tol = wcs_copy(gal, nx, ny, crpix1, crpix2)) > 1.0e-4) {
                    throw exception_failure("GAL copy tolerance 1.0e-4 exceeded: "+gammalib::str(tol));
                }
                test_try_success();
            }
            catch (std::exception &e) {
                test_try_failure(e);
            }

            // Free memory
            delete cel;
            delete gal;

            // Signal success
            test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }

    } // endfor: looped over all non-HealPix projections

    // Exit test
    return;
}


/***************************************************************************
 * @brief GSkymap_healpix_construct
 ***************************************************************************/
void TestGSky::test_GSkymap_healpix_construct(void)
{

    // Test void constructor
    test_try("Test void constructor");
    try {
        GSkymap map;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test correct Healpix constructors
    test_try("Test correct Healpix constructors");
    try {
        GSkymap ring1("GAL", 1, "RING", 1);
        GSkymap ring2("GAL", 2, "RING", 1);
        GSkymap ring4("GAL", 4, "RING", 1);
        GSkymap ring8("GAL", 8, "RING", 1);
        GSkymap ring16("GAL", 16, "RING", 1);
        GSkymap ring32("GAL", 32, "RING", 1);
        GSkymap ring64("GAL", 64, "RING", 1);
        GSkymap ring128("GAL", 128, "RING", 1);
        GSkymap ring256("GAL", 256, "RING", 1);
        GSkymap ring512("GAL", 512, "RING", 1);
        GSkymap nest1("GAL", 1, "NEST", 1);
        GSkymap nest2("GAL", 2, "NEST", 1);
        GSkymap nest4("GAL", 4, "NEST", 1);
        GSkymap nest8("GAL", 8, "NEST", 1);
        GSkymap nest16("GAL", 16, "NEST", 1);
        GSkymap nest32("GAL", 32, "NEST", 1);
        GSkymap nest64("GAL", 64, "NEST", 1);
        GSkymap nest128("GAL", 128, "NEST", 1);
        GSkymap nest256("GAL", 256, "NEST", 1);
        GSkymap nest512("GAL", 512, "NEST", 1);
        GSkymap map1("CEL", 1, "RING", 1);
        GSkymap map2("CEL", 1, "RING", 2);
        GSkymap map3("EQU", 1, "RING", 1);
        GSkymap map4("EQU", 1, "NESTED", 1);

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test Healpix copy constructor
    test_try("Test Healpix copy constructor");
    try {
        GSkymap ring1("GAL", 1, "RING", 1);
        GSkymap ring2 = ring1;

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test invalid coordsys in constructor
    test_try("Test invalid coordsys in constructor");
    try {
        GSkymap map("HOR", 1, "RING", 1);
        test_try_failure();
    }
    catch (GException::wcs_bad_coords &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test invalid nside in constructor
    test_try("Test invalid nside in constructor");
    try {
        GSkymap map("GAL", 3, "RING", 1);
        test_try_failure();
    }
    catch (GException::wcs_hpx_bad_nside &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test invalid ordering in constructor
    test_try("Test invalid ordering in constructor");
    try {
        GSkymap map("GAL", 2, "SPHERICAL", 1);
        test_try_failure();
    }
    catch (GException::wcs_hpx_bad_ordering &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }


    // Test invalid nmaps in constructor
    test_try("Test invalid nmaps in constructor");
    try {
        GSkymap map("GAL", 2, "NEST", 0);
        test_try_failure();
    }
    catch (GException::skymap_bad_par &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;

}


/***************************************************************************
 * @brief GSkymap_healpix_io
 ***************************************************************************/
void TestGSky::test_GSkymap_healpix_io(void)
{
    // Set filenames
    const std::string file1 = "test_skymap_hpx_1.fits";

    // Define Healpix map for comparison
    GSkymap refmap("GAL", 4, "RING", 1);

    // Test Healpix map saving
    test_try("Test Healpix map saving");
    try {
        for (int pix = 0; pix < refmap.npix(); ++pix)
            refmap(pix) = pix+1;
        refmap.save(file1, true);

        test_try_success(); 
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test Healpix map loading
    test_try("Test Healpix map loading");
    try {
        GSkymap map;
        map.load(file1);
        int diff = 0;
        for (int pix = 0; pix < refmap.npix(); ++pix) {
            if (map(pix) != refmap(pix))
                diff++;
        }
        if (diff > 0) {
            throw exception_error("Loaded file differs from saved file.");
        }

        test_try_success(); 
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test Healpix map instantiation
    test_try("Test Healpix map instantiation");
    try {
        GSkymap map(file1);
        int diff = 0;
        for (int pix = 0; pix < refmap.npix(); ++pix) {
            if (map(pix) != refmap(pix))
                diff++;
        }
        if (diff > 0) {
            throw exception_error("Loaded file differs from saved file.");
        }
        test_try_success(); 
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***************************************************************************
 * @brief GSkymap_wcs_construct
 ***************************************************************************/
void TestGSky::test_GSkymap_wcs_construct(void)
{
    // Set precision
    double eps = 1.0e-5;

    // Test void constructor
    test_try("Test void constructor");
    try {
        GSkymap map;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test non-Healpix constructors
    test_try("Test non-Healpix constructors");
    try {
        GSkymap map1("CAR", "GAL", 0.0, 0.0, 1.0, 1.0, 100, 100);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test CAR projection
    test_try("Test CAR projection");
    try {
        GSkymap map1("CAR", "GAL", 138.817, 37.293, 0.521, 0.931, 100, 100);
        GSkyDir dir;
        for (int l = -180; l < 180; ++l) {
            for (int b = -90; b < 90; ++b) {
                dir.lb_deg(double(l),double(b));
                GSkyPixel pixel    = map1.dir2pix(dir);
                GSkyDir   dir_back = map1.pix2dir(pixel);
                double    dist     = dir.dist_deg(dir_back);
                if (dist > eps) {
                    throw exception_failure("Sky direction differs: dir="+dir.print()+" pixel="+pixel.print()+" dir_back"+ dir_back.print()+" dist="+gammalib::str(dist)+" deg");
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;

}


/***************************************************************************
 * @brief GSkymap_wcs_io
 ***************************************************************************/
void TestGSky::test_GSkymap_wcs_io(void)
{
    // Set filenames
    const std::string file1 = "test_skymap_wcs_1.fits";
    const std::string file2 = "test_skymap_wcs_2.fits";

    // Define WCS map for comparison
    GSkymap refmap1("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 100, 100);
    GSkymap refmap2("CAR", "CEL", 0.0, 0.0, -0.3, 0.3, 200, 200);

    // Test WCS map saving
    test_try("Test WCS map saving");
    try {
        for (int pix = 0; pix < refmap1.npix(); ++pix) {
            refmap1(pix) = pix+1;
        }
        refmap1.save(file1, 1);
        for (int pix = 0; pix < refmap2.npix(); ++pix) {
            refmap2(pix) = pix+1;
        }
        refmap2.save(file2, 1);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test WCS map loading (1)
    test_try("Test WCS map loading (1)");
    try {
        GSkymap map;
        map.load(file1);
        int diff = 0;
        for (int pix = 0; pix < refmap1.npix(); ++pix) {
            if (map(pix) != refmap1(pix)) {
                diff++;
            }
        }
        if (diff > 0) {
            throw exception_error("Loaded file "+file1+" differs from saved file .");
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test WCS map loading (2)
    test_try("Test WCS map loading (2)");
    try {
        GSkymap map;
        map.load(file2);
        int diff = 0;
        for (int pix = 0; pix < refmap2.npix(); ++pix) {
            if (map(pix) != refmap2(pix)) {
                diff++;
            }
        }
        if (diff > 0) {
            throw exception_error("Loaded file "+file2+" differs from saved file .");
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***************************************************************************
 * @brief GSkymap
 ***************************************************************************/
void TestGSky::test_GSkymap(void)
{
    // Define several maps for comparison
    GSkymap map_src("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 10, 10, 2);
    GSkymap map_dst("CAR", "GAL", 0.0, 0.0, -0.1, 0.1, 100, 100, 2);

    // Fill map pixels
    double total_src = 0.0;
    for (int pix = 0; pix < map_src.npix(); ++pix) {
        map_src(pix) = pix+1;
        total_src   += map_src(pix);
    }

    // Add pixels to destination map
    map_dst += map_src;

    // Check total in destination map. Note that the total in the destination
    // map should be 100 times the total in the source maps as the operator
    // is expected to perform strict bi-linear interpolation
    double total_dst = 0.0;
    for (int pix = 0; pix < map_dst.npix(); ++pix) {
        total_dst   += map_dst(pix);
    }
    total_dst /= 100.0;
	test_value(total_dst, total_src, 1.0e-3, "Test operator+");

    // Subtract pixels from destination map
    map_dst -= map_src;

    // Check total in destination map. Note that the total in the destination
    // map should be zero.
    total_dst = 0.0;
    for (int pix = 0; pix < map_dst.npix(); ++pix) {
        total_dst   += map_dst(pix);
    }
	test_value(total_dst, 0.0, 1.0e-3, "Test operator+");

    // Save maps
    map_src.save("test_map_src.fits", true);
    map_dst.save("test_map_dst.fits", true);

    // Exit test
    return;
}


/***************************************************************************
 * @brief GSkyRegionCircle_construct
 ***************************************************************************/
void TestGSky::test_GSkyRegionCircle_construct(void)
{
	// Define region for comparison
	GSkyDir refdir_radeczerozero = GSkyDir();
	refdir_radeczerozero.radec_deg(0,0);
	double refradius = 10.;

    // Test constructing:
    test_try("Test constructor");
    try {
        GSkyRegionCircle circle(refdir_radeczerozero,refradius);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test constructing with radius 0
    test_try("Test constructor2");
    try {
		GSkyRegionCircle circle(refdir_radeczerozero,-1);
		test_try_failure();
	}
	catch (GException::invalid_argument &e) {
        test_try_success();
	}
	catch (std::exception &e) {
        test_try_failure(e);
	}

    // Test constructing with radius 0
    test_try("Test radius assignment after");
    try {
		GSkyRegionCircle circle(refdir_radeczerozero,refradius);
		circle.radius(-1.0);
		test_try_failure();
	}
	catch (GException::invalid_argument &e) {
        test_try_success();
	}
	catch (std::exception &e) {
        test_try_failure(e);
	}

	// Check radius assignment
	GSkyRegionCircle refregion(refdir_radeczerozero,refradius);
	double refradius_check = refregion.radius();
	test_value(refradius,refradius_check,1.0e-10, "Test radius assignment");

	// Check solid angle assignment
	double solidangle_check = refregion.solidangle();
	double solidangle = 2*gammalib::pi*(1- std::cos(refradius /180 * gammalib::pi));
	test_value(solidangle_check,solidangle,1.0e-10, "Test solid angle assignment");

	//exit test
    return;
}

/***************************************************************************
 * @brief GSkyRegions_logic
 ***************************************************************************/
void TestGSky::test_GSkyRegionCircle_logic(void)
{

	// set input files /filenames
//    const std::string  file1 = "test_GSkyRegionCircle1.reg";
//    const std::string  file2 = "test_GSkyRegionCircle2.reg";

    //set reference values

    GSkyDir refdir_radeczerozero = GSkyDir();
    refdir_radeczerozero.radec_deg(0,0);

    GSkyDir refdir_lbzerozero = GSkyDir();
    refdir_lbzerozero.lb_deg(0,0);

    GSkyDir refdir_outside_refregion = GSkyDir();
    refdir_outside_refregion.radec_deg(100,80);

    GSkyDir refdir_raoffset = GSkyDir();
    refdir_raoffset.radec_deg(5,0);

    GSkyDir refdir_decpole = GSkyDir();
    refdir_decpole.radec_deg(0,89);

    GSkyDir refdir_ndecpole = GSkyDir();
	refdir_ndecpole.radec_deg(100,89);

    GSkyDir refdir_rapole = GSkyDir();
    refdir_rapole.radec_deg(359,0);

    GSkyRegionCircle refregion(refdir_radeczerozero,10);
    GSkyDir centre = refregion.centre();

    GSkyRegionCircle refregion_smaller(refdir_radeczerozero,5);
    GSkyRegionCircle refregion_larger(refdir_radeczerozero,20);
    GSkyRegionCircle refregion_raoffset(refdir_raoffset,10);
    GSkyRegionCircle refregion_rapole(refdir_rapole,3);
    GSkyRegionCircle refregion_decpole(refdir_decpole,3);

    // Test contain dirs
	test_assert(refregion.contains(refdir_radeczerozero),"test for containment");
    test_assert(refregion.contains(centre), "test if centre in circle");
	test_assert(!refregion.contains(refdir_outside_refregion), "test2 for containment");

	// Test contain regions
	test_assert(refregion.contains(refregion_smaller),"test for containment region");
	test_assert(!refregion.contains(refregion_larger), "test for containment region2 ");
	test_assert(refregion.contains(refregion), "test3 for containment region");

	test_assert(refregion.contains(refdir_rapole), "rapole for containment region");
	test_assert(refregion_decpole.contains(refdir_ndecpole), "rapole for containment region");

    // Test overlaps
	test_assert(refregion.overlaps(refregion_smaller),"test for overlap");
	test_assert(refregion.overlaps(refregion_larger),"test2 for overlap");
	test_assert(refregion.overlaps(refregion_raoffset),"test3 for overlap");
	test_assert(!refregion.overlaps(refregion_decpole),"test4 for overlap");

	// Exit test
    return;
}


/***************************************************************************
 * @brief Test GSkyRegions input and output
 ***************************************************************************/
 void TestGSky::test_GSkyRegions_io(void)
{
	// Set filenames
	const std::string filename = "data/test_circle_region.reg";

    // Allocate regions
    GSkyRegions regions;

	// Test regions loading
	test_try("Test regions loading");
	try {
		regions.load(filename);
		test_try_success();
	}
	catch (std::exception &e) {
		test_try_failure(e);
	}

    // Check if region was loaded correctly
    GSkyRegionCircle* circle = dynamic_cast<GSkyRegionCircle*>(regions[0]);
    test_assert(circle->type() == "Circle", "Region is not a circle");
    test_value(circle->radius(), 10.0, 1.0e-10);
    test_value(circle->ra(), 0.1, 1.0e-10);
    test_value(circle->dec(), -35.6, 1.0e-10);

	// Test regions saving
	test_try("Test regions saving");
	try {
		regions.save("region.reg");
		test_try_success();
	}
	catch (std::exception &e) {
		test_try_failure(e);
	}

	// Test regions reloading
	test_try("Test regions reloading");
	try {
		regions.load("region.reg");
		test_try_success();
	}
	catch (std::exception &e) {
		test_try_failure(e);
	}

    // Check if region was loaded correctly
    circle = dynamic_cast<GSkyRegionCircle*>(regions[0]);
    test_assert(circle->type() == "Circle", "Region is not a circle");
    test_value(circle->radius(), 10.0, 1.0e-10);
    test_value(circle->ra(), 0.1, 1.0e-10);
    test_value(circle->dec(), -35.6, 1.0e-10);

    // Exit test
    return;
        
}

/***************************************************************************
 * @brief Test GHorizDir class
 ***************************************************************************/
void TestGSky::test_GHorizDir(void){

    // Empty horizontal direction
    GHorizDir nulldir;

    // Set position
    GHorizDir altaz;
    double testalt =  45.0;
    double testaz  = 180.0;
    altaz.altaz_deg(testalt, testaz);

    // Retrieve position
    double newalt = altaz.alt_deg();
    double newaz  = altaz.az_deg();

    // Test position
    test_value(newalt, testalt, 1e-10, "Altitude" );
    test_value(newaz,  testaz,  1e-10, "Azimuth" );
    test_value(altaz.zenith_deg(), 90.0-testalt, 1e-10, "Zenith angle");

    // Return
    return;
}          
          

/***************************************************************************
 * @brief Main test function
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuites("GSky");
           
    bool was_successful = true;
            
    // Create a test suite
    TestGSky test;
            
    // Append to the container
    testsuites.append(test);

    // Run
    was_successful=testsuites.run();

    // Save xml report
    testsuites.save("reports/GSky.xml");

    // Return
    return was_successful ? 0:1;
}
