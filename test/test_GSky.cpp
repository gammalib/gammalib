/***************************************************************************
 *                  test_GSky.cpp - Test sky module                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
#include <cstdlib>                            // for system
#include "test_GSky.hpp"
#include "GTools.hpp"


/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_WCS_FORTH_BACK_PIXEL_DEBUG
#define G_WCS_COPY_DEBUG

/* __ Constants __________________________________________________________ */
const std::string datadir         = std::getenv("TEST_DATA");
const std::string sky_region      = datadir + "/test_circle_region.reg";
const std::string sky_region_rect = datadir + "/test_rect_region.reg";
const std::string sky_region_map  = datadir + "/test_map_region.fits";
const std::string sky_map         = datadir + "/cena_lobes_parkes.fits";


/***********************************************************************//**
 * @brief Set parameters and tests
 ***************************************************************************/
void TestGSky::set(void){

    // Test name
    name("GSky");

    // Append tests
    append(static_cast<pfunction>(&TestGSky::test_GWcs),
           "Test GWcs");
    append(static_cast<pfunction>(&TestGSky::test_GSkyDir),
           "Test GSkyDir");
    append(static_cast<pfunction>(&TestGSky::test_GSkyDirs),
           "Test GSkyDirs");
    append(static_cast<pfunction>(&TestGSky::test_GSkyPixel),
           "Test GSkyPixel");
    append(static_cast<pfunction>(&TestGSky::test_GSkyMap_healpix_construct),
           "Test Healpix GSkyMap constructors");
    append(static_cast<pfunction>(&TestGSky::test_GSkyMap_healpix_io),
           "Test Healpix GSkyMap I/O");
    append(static_cast<pfunction>(&TestGSky::test_GSkyMap_wcs_construct),
           "Test WCS GSkyMap constructors");
    append(static_cast<pfunction>(&TestGSky::test_GSkyMap_wcs_io),
           "Test WCS GSkyMap I/O");
    append(static_cast<pfunction>(&TestGSky::test_GSkyMap),
           "Test GSkyMap");
    append(static_cast<pfunction>(&TestGSky::test_GSkyMap_io),
           "Test GSkyMap I/O");
    append(static_cast<pfunction>(&TestGSky::test_GSkyRegions),
           "Test GSkyRegions");
    append(static_cast<pfunction>(&TestGSky::test_GSkyRegionCircle),
           "Test GSkyRegionCircle");
    append(static_cast<pfunction>(&TestGSky::test_GSkyRegionRect),
           "Test GSkyRegionRect");
    append(static_cast<pfunction>(&TestGSky::test_GSkyRegionMap),
           "Test GSkyRegionMap");
    append(static_cast<pfunction>(&TestGSky::test_GHorizDir),
           "Test GHorizDir");

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
double TestGSky::wcs_forth_back_pixel(GWcs*   wcs,
                                      int     nx,
                                      int     ny,
                                      double& crpix1,
                                      double& crpix2)
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
            double y = double(iy) + crpix2;

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
                std::cout << wcs->code() << ":";
                std::cout << " dist=" << dist;
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
double TestGSky::wcs_copy(GWcs*   wcs,
                          int     nx,
                          int     ny,
                          double& crpix1,
                          double& crpix2)
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


/***************************************************************************
 * @brief Test GSkyDir class
 ***************************************************************************/
void TestGSky::test_GSkyDir(void)
{
    // Test void constructor
    GSkyDir dir1;
    test_value(dir1.ra_deg(),  0.0, "Right Ascension of empty sky direction");
    test_value(dir1.dec_deg(), 0.0, "Declination of empty sky direction");
    test_value(dir1.l_deg(),   0.0, "Longitude of empty sky direction");
    test_value(dir1.b_deg(),   0.0, "Latitude of empty sky direction");

    // Test Crab sky direction
    GSkyDir dir2;
    dir2.radec_deg(83.633083, 22.0145);
    test_value(dir2.ra_deg(),  83.633083, "Right Ascension of Crab");
    test_value(dir2.dec_deg(), 22.0145,   "Declination of Crab");
    test_value(dir2.l_deg(),   184.55745, "Longitude of Crab");
    test_value(dir2.b_deg(),    -5.78436, "Latitude of Crab");

    // Test precess() method
    dir2.precess(2000.0, 1950.0);
    test_value(dir2.ra_deg(),  82.88086799, "Right Ascension of Crab (epoch 1950)");
    test_value(dir2.dec_deg(), 21.98181049, "Declination of Crab (epoch 1950)");
    dir2.precess(1950.0, 2000.0);
    test_value(dir2.ra_deg(),  83.633083, "Right Ascension of Crab");
    test_value(dir2.dec_deg(), 22.0145,   "Declination of Crab");

    // Test posang() and posang_deg() methods
    dir1.radec_deg(0.0,  0.0);
    dir2.radec_deg(0.0, 10.0);
    test_value(dir1.posang(dir2), 0.0, "1. Position angle posang()");
    test_value(dir1.posang(dir2, "CEL"), 0.0, "2. Position angle posang(CEL)");
    test_value(dir1.posang(dir2, "GAL"), 0.4098000, "3. Position angle posang(GAL)");
    test_value(dir1.posang_deg(dir2), 0.0, "4. Position angle posang_deg()");
    test_value(dir1.posang_deg(dir2, "CEL"), 0.0, "5. Position angle posang_deg(CEL)");
    test_value(dir1.posang_deg(dir2, "GAL"), 23.47981, "6. Position angle posang_deg(GAL)");
    dir2.radec_deg(90.0, 0.0);
    test_value(dir1.posang(dir2), 1.5707963, "7. Position angle posang()");
    test_value(dir1.posang(dir2, "CEL"), 1.5707963, "8. Position angle posang(CEL)");
    test_value(dir1.posang(dir2, "GAL"), 1.9805963, "9. Position angle posang(GAL)");
    test_value(dir1.posang_deg(dir2), 90.0, "10. Position angle posang_deg()");
    test_value(dir1.posang_deg(dir2, "CEL"), 90.0, "11. Position angle posang_deg(CEL)");
    test_value(dir1.posang_deg(dir2, "GAL"), 113.47981, "12. Position angle posang_deg(GAL)");
    dir1.lb_deg(0.0,  0.0);
    dir2.lb_deg(0.0, 10.0);
    test_value(dir1.posang(dir2), -1.0227397, "13. Position angle posang()");
    test_value(dir1.posang(dir2, "CEL"), -1.0227397, "14. Position angle posang(CEL)");
    test_value(dir1.posang(dir2, "GAL"), 0.0, "15. Position angle posang(GAL)");
    test_value(dir1.posang_deg(dir2), -58.5986663, "16. Position angle posang_deg()");
    test_value(dir1.posang_deg(dir2, "CEL"), -58.5986663, "17. Position angle posang_deg(CEL)");
    test_value(dir1.posang_deg(dir2, "GAL"), 0.0, "18. Position angle posang_deg(GAL)");
    dir2.lb_deg(90.0, 0.0);
    test_value(dir1.posang(dir2), 0.5480567, "19. Position angle posang()");
    test_value(dir1.posang(dir2, "CEL"), 0.5480567, "20. Position angle posang(CEL)");
    test_value(dir1.posang(dir2, "GAL"), 1.5707963, "21. Position angle posang(GAL)");
    test_value(dir1.posang_deg(dir2), 31.4013337, "22. Position angle posang_deg()");
    test_value(dir1.posang_deg(dir2, "CEL"), 31.4013337, "23. Position angle posang_deg(CEL)");
    test_value(dir1.posang_deg(dir2, "GAL"), 90.0, "24. Position angle posang_deg(GAL)");

    // Return
    return;
}


/***************************************************************************
 * @brief Test GSkyDirs class
 ***************************************************************************/
void TestGSky::test_GSkyDirs(void)
{
    // Manipulate GSkyDirs starting from an empty object
    GSkyDirs dirs;
    test_value(dirs.size(), 0, "GSkyDirs should have zero size.");
    test_assert(dirs.is_empty(), "GSkyDirs should be empty.");

    // Add a sky direction
    dirs.append(GSkyDir());
    test_value(dirs.size(), 1, "GSkyDirs should have 1 sky direction.");
    test_assert(!dirs.is_empty(), "GSkyDirs should not be empty.");

    // Remove sky direction
    dirs.remove(0);
    test_value(dirs.size(), 0, "GSkyDirs should have zero size.");
    test_assert(dirs.is_empty(), "GSkyDirs should be empty.");

    // Append two sky directions
    dirs.append(GSkyDir());
    dirs.append(GSkyDir());
    test_value(dirs.size(), 2, "GSkyDirs should have 2 sky directions.");
    test_assert(!dirs.is_empty(), "GSkyDirs should not be empty.");

    // Clear object
    dirs.clear();
    test_value(dirs.size(), 0, "GSkyDirs should have zero size.");
    test_assert(dirs.is_empty(), "GSkyDirs should be empty.");

    // Insert two sky directions
    dirs.insert(0, GSkyDir());
    dirs.insert(0, GSkyDir());
    test_value(dirs.size(), 2, "GSkyDirs should have 2 sky directions.");
    test_assert(!dirs.is_empty(), "GSkyDirs should not be empty.");

    // Extend sky directios
    dirs.extend(dirs);
    test_value(dirs.size(), 4, "GSkyDirs should have 4 sky directions.");
    test_assert(!dirs.is_empty(), "GSkyDirs should not be empty.");

    // Create 4 sky directions
    dirs.clear();
    for (int i = 0; i < 4; ++i) {
        GSkyDir dir;
        dir.radec_deg(double(i), double(i));
        dirs.append(dir);
    }
    for (int i = 0; i < 4; ++i) {
        test_value(dirs[i].ra_deg(), double(i));
        test_value(dirs[i].dec_deg(), double(i));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSkyPixel class
 *
 * Test GSkyPixel class.
 ***************************************************************************/
void TestGSky::test_GSkyPixel(void)
{
    // Test void constructor
    GSkyPixel pixel1;
    test_value(pixel1.size(), 0);

    // Test 1D constructor
    GSkyPixel pixel2(41);
    test_assert(pixel2.is_1D(), "Pixel is not 1D");
    test_assert(!pixel2.is_2D(), "Pixel is 2D but it should be 1D");
    test_value(pixel2.size(), 1);
    test_value(pixel2.index(), 41.0);

    // Test 1D constructor
    GSkyPixel pixel3(41.001);
    test_assert(pixel3.is_1D(), "Pixel is not 1D");
    test_assert(!pixel3.is_2D(), "Pixel is 2D but it should be 1D");
    test_value(pixel3.size(), 1);
    test_value(pixel3.index(), 41.001);

    // Test 2D constructor
    GSkyPixel pixel4(41,14);
    test_assert(!pixel4.is_1D(), "Pixel is 1D but it should be 2D");
    test_assert(pixel4.is_2D(), "Pixel is not 2D");
    test_value(pixel4.size(), 2);
    test_value(pixel4.x(), 41.0);
    test_value(pixel4.y(), 14.0);

    // Test 2D constructor
    GSkyPixel pixel5(41.001,14.003);
    test_assert(!pixel5.is_1D(), "Pixel is 1D but it should be 2D");
    test_assert(pixel5.is_2D(), "Pixel is not 2D");
    test_value(pixel5.size(), 2);
    test_value(pixel5.x(), 41.001);
    test_value(pixel5.y(), 14.003);

    // Test 1D indexing
    GSkyPixel pixel6;
    test_value(pixel6.size(), 0);
    pixel6 = 41;
    int index6 = pixel6;
    test_value(index6, 41);
    pixel6 = 41.001;
    index6 = pixel6;
    test_value(index6, 41);
    double dindex6 = pixel6;
    test_value(dindex6, 41.001);

    // Test cloning
    GSkyPixel pixel7;
    test_value(pixel7.size(), 0);
    GSkyPixel* clone = pixel7.clone();
    test_value(clone->size(), 0);
    delete clone;
    pixel7 = 41;
    clone = pixel7.clone();
    test_value(clone->size(), 1);
    int index7 = pixel7;
    test_value(index7, 41);
    delete clone;
    pixel7.x(41.001);
    pixel7.y(14.003);
    clone = pixel7.clone();
    test_value(clone->size(), 2);
    test_value(clone->x(), 41.001);
    test_value(clone->y(), 14.003);
    delete clone;

    // Test equality and inequality operators using pixels from above
    test_assert(pixel1 == pixel1, "Test equality operator for empty pixel");
    test_assert(!(pixel1 != pixel1), "Test inequality operator for empty pixel");
    test_assert(pixel2 == pixel2, "Test equality operator for 1D pixel");
    test_assert(!(pixel2 != pixel2), "Test inequality operator for 1D pixel");
    test_assert(pixel4 == pixel4, "Test equality operator for 2D pixel");
    test_assert(!(pixel4 != pixel4), "Test inequality operator for 2D pixel");
    test_assert(!(pixel2 == pixel4), "Test equality operator for 1D and 2D pixels");
    test_assert(pixel2 != pixel4, "Test inequality operator for 1D and 2D pixels");

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
    double cdelt1 =   0.1;
    double cdelt2 =   0.1;
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
                    throw exception_failure("CEL forth-back transformation"
                                            " tolerance 1.0e-10 exceeded: "+
                                            gammalib::str(tol));
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
                    throw exception_failure("GAL forth-back transformation"
                                            " tolerance 1.0e-10 exceeded: "+
                                            gammalib::str(tol));
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
                    throw exception_failure("CEL copy tolerance 1.0e-4"
                                            " exceeded: "+gammalib::str(tol));
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
                    throw exception_failure("GAL copy tolerance 1.0e-4"
                                            " exceeded: "+gammalib::str(tol));
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
 * @brief GSkyMap_healpix_construct
 ***************************************************************************/
void TestGSky::test_GSkyMap_healpix_construct(void)
{

    // Test void constructor
    GSkyMap empty_map;
    test_assert(empty_map.is_empty(), "Check for empty sky map");
    test_value(empty_map.npix(), 0, "Check for no sky map pixels");
    test_value(empty_map.nx(), 0, "Check for no X-direction map pixels");
    test_value(empty_map.ny(), 0, "Check for no Y-direction map pixels");
    test_value(empty_map.nmaps(), 0, "Check for no sky maps");

    // Test Healpix constructors
    GSkyMap ring1("GAL", 1, "RING", 1);
    test_assert(!ring1.is_empty(), "Check for non-empty sky map");
    test_value(ring1.npix(), 12, "Check for 12 sky map pixels");
    test_value(ring1.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring1.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring1.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring2("GAL", 2, "RING", 1);
    test_assert(!ring2.is_empty(), "Check for non-empty sky map");
    test_value(ring2.npix(), 48, "Check for 48 sky map pixels");
    test_value(ring2.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring2.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring2.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring4("GAL", 4, "RING", 1);
    test_assert(!ring4.is_empty(), "Check for non-empty sky map");
    test_value(ring4.npix(), 192, "Check for 192 sky map pixels");
    test_value(ring4.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring4.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring4.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring8("GAL", 8, "RING", 1);
    test_assert(!ring8.is_empty(), "Check for non-empty sky map");
    test_value(ring8.npix(), 768, "Check for 768 sky map pixels");
    test_value(ring8.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring8.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring8.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring16("GAL", 16, "RING", 1);
    test_assert(!ring16.is_empty(), "Check for non-empty sky map");
    test_value(ring16.npix(), 3072, "Check for 3072 sky map pixels");
    test_value(ring16.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring16.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring16.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring32("GAL", 32, "RING", 1);
    test_assert(!ring32.is_empty(), "Check for non-empty sky map");
    test_value(ring32.npix(), 12288, "Check for 12288 sky map pixels");
    test_value(ring32.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring32.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring32.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring64("GAL", 64, "RING", 1);
    test_assert(!ring64.is_empty(), "Check for non-empty sky map");
    test_value(ring64.npix(), 49152, "Check for 49152 sky map pixels");
    test_value(ring64.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring64.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring64.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring128("GAL", 128, "RING", 1);
    test_assert(!ring128.is_empty(), "Check for non-empty sky map");
    test_value(ring128.npix(), 196608, "Check for 196608 sky map pixels");
    test_value(ring128.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring128.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring128.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring256("GAL", 256, "RING", 1);
    test_assert(!ring256.is_empty(), "Check for non-empty sky map");
    test_value(ring256.npix(), 786432, "Check 786432 XX sky map pixels");
    test_value(ring256.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring256.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring256.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap ring512("GAL", 512, "RING", 1);
    test_assert(!ring512.is_empty(), "Check for non-empty sky map");
    test_value(ring512.npix(), 3145728, "Check for 3145728 sky map pixels");
    test_value(ring512.nx(), 0, "Check for no X-direction map pixels");
    test_value(ring512.ny(), 0, "Check for no Y-direction map pixels");
    test_value(ring512.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest1("GAL", 1, "NEST", 1);
    test_assert(!nest1.is_empty(), "Check for non-empty sky map");
    test_value(nest1.npix(), 12, "Check for 12 sky map pixels");
    test_value(nest1.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest1.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest1.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest2("GAL", 2, "NEST", 1);
    test_assert(!nest2.is_empty(), "Check for non-empty sky map");
    test_value(nest2.npix(), 48, "Check for 48 sky map pixels");
    test_value(nest2.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest2.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest2.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest4("GAL", 4, "NEST", 1);
    test_assert(!nest4.is_empty(), "Check for non-empty sky map");
    test_value(nest4.npix(), 192, "Check for 192 sky map pixels");
    test_value(nest4.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest4.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest4.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest8("GAL", 8, "NEST", 1);
    test_assert(!nest8.is_empty(), "Check for non-empty sky map");
    test_value(nest8.npix(), 768, "Check for 768 sky map pixels");
    test_value(nest8.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest8.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest8.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest16("GAL", 16, "NEST", 1);
    test_assert(!nest16.is_empty(), "Check for non-empty sky map");
    test_value(nest16.npix(), 3072, "Check for 3072 sky map pixels");
    test_value(nest16.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest16.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest16.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest32("GAL", 32, "NEST", 1);
    test_assert(!nest32.is_empty(), "Check for non-empty sky map");
    test_value(nest32.npix(), 12288, "Check for 12288 sky map pixels");
    test_value(nest32.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest32.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest32.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest64("GAL", 64, "NEST", 1);
    test_assert(!nest64.is_empty(), "Check for non-empty sky map");
    test_value(nest64.npix(), 49152, "Check for 49152 sky map pixels");
    test_value(nest64.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest64.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest64.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest128("GAL", 128, "NEST", 1);
    test_assert(!nest128.is_empty(), "Check for non-empty sky map");
    test_value(nest128.npix(), 196608, "Check for 196608 sky map pixels");
    test_value(nest128.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest128.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest128.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest256("GAL", 256, "NEST", 1);
    test_assert(!nest256.is_empty(), "Check for non-empty sky map");
    test_value(nest256.npix(), 786432, "Check for 786432 sky map pixels");
    test_value(nest256.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest256.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest256.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap nest512("GAL", 512, "NEST", 1);
    test_assert(!nest512.is_empty(), "Check for non-empty sky map");
    test_value(nest512.npix(), 3145728, "Check for 3145728 sky map pixels");
    test_value(nest512.nx(), 0, "Check for no X-direction map pixels");
    test_value(nest512.ny(), 0, "Check for no Y-direction map pixels");
    test_value(nest512.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap map1("CEL", 1, "RING", 1);
    test_assert(!map1.is_empty(), "Check for non-empty sky map");
    test_value(map1.npix(), 12, "Check for 12 sky map pixels");
    test_value(map1.nx(), 0, "Check for no X-direction map pixels");
    test_value(map1.ny(), 0, "Check for no Y-direction map pixels");
    test_value(map1.nmaps(), 1, "Check for one sky map");
    //
    GSkyMap map2("CEL", 1, "NESTED", 2);
    test_assert(!map2.is_empty(), "Check for non-empty sky map");
    test_value(map2.npix(), 12, "Check for 12 sky map pixels");
    test_value(map2.nx(), 0, "Check for no X-direction map pixels");
    test_value(map2.ny(), 0, "Check for no Y-direction map pixels");
    test_value(map2.nmaps(), 2, "Check for one sky map");

    // Test Healpix copy constructor
    GSkyMap ring1cpy = ring1;
	test_assert(!ring1cpy.is_empty(), "Check for non-empty sky map");
	test_value(ring1cpy.npix(), 12, "Check for 12 sky map pixels");
	test_value(ring1cpy.nx(), 0, "Check for no X-direction map pixels");
	test_value(ring1cpy.ny(), 0, "Check for no Y-direction map pixels");
	test_value(ring1cpy.nmaps(), 1, "Check for one sky map");

    // Test invalid coordsys in constructor
    test_try("Test invalid coordsys in constructor");
    try {
        GSkyMap map("HOR", 1, "RING", 1);
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
        GSkyMap map("GAL", 3, "RING", 1);
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
        GSkyMap map("GAL", 2, "SPHERICAL", 1);
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
        GSkyMap map("GAL", 2, "NEST", 0);
        test_try_failure();
    }
    catch (GException::skymap_bad_par &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test region circle
    GSkyRegionCircle region = ring1.region_circle();
	test_value(region.centre().ra_deg(), 0.0, "Check R.A. for region_circle()");
	test_value(region.centre().dec_deg(), 0.0, "Check Dec. for region_circle()");
	test_value(region.radius(), 180.0, "Check radius for region_circle()");

    // Exit test
    return;

}


/***************************************************************************
 * @brief GSkyMap_healpix_io
 ***************************************************************************/
void TestGSky::test_GSkyMap_healpix_io(void)
{
    // Set filenames
    const std::string file1 = "test_skymap_hpx_1.fits";

    // Define Healpix map for comparison
    GSkyMap refmap("GAL", 4, "RING", 1);

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
        GSkyMap map;
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
        GSkyMap map(file1);
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
 * @brief GSkyMap_wcs_construct
 ***************************************************************************/
void TestGSky::test_GSkyMap_wcs_construct(void)
{
    // Set precision
    const double eps = 1.0e-5;

    // Test void constructor
    test_try("Test void constructor");
    try {
        GSkyMap map;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test non-Healpix constructors
    test_try("Test non-Healpix constructors");
    try {
        GSkyMap map1("CAR", "GAL", 0.0, 0.0, 1.0, 1.0, 100, 100);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test CAR projection
    test_try("Test CAR projection");
    try {
        GSkyMap map1("CAR", "GAL", 138.817, 37.293, 0.521, 0.931, 100, 100);
        GSkyDir dir;
        for (int l = -180; l < 180; ++l) {
            for (int b = -90; b < 90; ++b) {
                dir.lb_deg(double(l),double(b));
                GSkyPixel pixel    = map1.dir2pix(dir);
                GSkyDir   dir_back = map1.pix2dir(pixel);
                double    dist     = dir.dist_deg(dir_back);
                if (dist > eps) {
                    throw exception_failure("Sky direction differs: dir="+
                                            dir.print()+" pixel="+
                                            pixel.print()+" dir_back"+
                                            dir_back.print()+" dist="+
                                            gammalib::str(dist)+" deg");
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test region circle
    GSkyMap          map2("CAR", "CEL", 138.817, 37.293, 0.521, 0.931, 100, 100);
    GSkyRegionCircle region = map2.region_circle();
	test_value(region.centre().ra_deg(), 139.146485, "Check R.A. for region_circle()");
	test_value(region.centre().dec_deg(), 37.758046, "Check Dec. for region_circle()");
	test_value(region.radius(), 52.370103, "Check radius for region_circle()");

    // Exit test
    return;
}


/***************************************************************************
 * @brief GSkyMap_wcs_io
 ***************************************************************************/
void TestGSky::test_GSkyMap_wcs_io(void)
{
    // Set filenames
    const std::string file1 = "test_skymap_wcs_1.fits";
    const std::string file2 = "test_skymap_wcs_2.fits";

    // Define WCS map for comparison
    GSkyMap refmap1("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 100, 100);
    GSkyMap refmap2("CAR", "CEL", 0.0, 0.0, -0.3, 0.3, 200, 200);

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
        GSkyMap map;
        map.load(file1);
        int diff = 0;
        for (int pix = 0; pix < refmap1.npix(); ++pix) {
            if (map(pix) != refmap1(pix)) {
                diff++;
            }
        }
        if (diff > 0) {
            throw exception_error("Loaded file "+file1+
                                  " differs from saved file.");
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test WCS map loading (2)
    test_try("Test WCS map loading (2)");
    try {
        GSkyMap map;
        map.load(file2);
        int diff = 0;
        for (int pix = 0; pix < refmap2.npix(); ++pix) {
            if (map(pix) != refmap2(pix)) {
                diff++;
            }
        }
        if (diff > 0) {
            throw exception_error("Loaded file "+file2+
                                  " differs from saved file.");
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
 * @brief GSkyMap
 *
 * This method tests the GSkyMap class.
 ***************************************************************************/
void TestGSky::test_GSkyMap(void)
{
    // Remove result files
    system("rm -rf test_map_src.fits");
    system("rm -rf test_map_dst.fits");
    system("rm -rf test_map_shape*.fits");

    // Define empty sky map
    GSkyMap empty_map;

    // Test that empty map is indeed empty
	test_assert(empty_map.is_empty(), "Check for empty sky map");
	test_value(empty_map.nmaps(), 0, "Check for no sky maps");
	test_value(empty_map.npix(), 0, "Check number of empty sky map pixels");
	test_value(empty_map.nx(), 0, "Check number of empty sky map X pixels");
	test_value(empty_map.ny(), 0, "Check number of empty sky map X pixels");
	test_value(empty_map.ndim(), 0, "Check empty sky map dimension");
	test_value(empty_map.shape().size(), 0, "Check empty sky map shape");

    // Test that writing, publishing and printing of empty sky map does not
    // lead to a segmentation fault
    empty_map.save("test_empty_map.fits", true);
    empty_map.publish("Empty sky map");
    empty_map.print();

    // Define several maps for comparison
    GSkyMap map_src("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 10, 10, 2);
    GSkyMap map_dst("CAR", "GAL", 0.0, 0.0, -0.1, 0.1, 100, 100, 2);
    GSkyMap map_new("CAR", "GAL", 0.0, 0.0, -0.1, 0.1, 100, 100, 2);

    // Test map dimensions and shape
	test_value(map_src.nmaps(), 2, "Check that sky map contains 2 maps");
	test_value(map_src.ndim(), 1, "Check that sky map has one dimension");
	test_value(map_src.shape().size(), 1, "Check that sky map has a shape size of 1");
	test_value(map_src.shape()[0], 2, "Check that sky map has 2 maps");

    // Fill map pixels
    double total_src = 0.0;
    for (int pix = 0; pix < map_src.npix(); ++pix) {
        for (int k = 0; k < map_src.nmaps(); ++k) {
            map_src(pix,k) = pix+1;
            total_src     += map_src(pix,k);
        }
    }

    // Add pixels to destination map
    map_new  = map_dst + map_src;
    map_dst += map_src;

    // Check total in destination map. Note that the total in the destination
    // map should be 100 times the total in the source map as the operator
    // is expected to perform strict bi-linear interpolation
    double total_dst = 0.0;
    double total_new = 0.0;
    double total_ref = 0.0;
    for (int pix = 0; pix < map_dst.npix(); ++pix) {
        for (int k = 0; k < map_dst.nmaps(); ++k) {
            total_dst += map_dst(pix,k);
            total_new += map_new(pix,k);
            total_ref += map_src(map_dst.pix2dir(pix),k);
        }
    }
    total_dst /= 100.0;
    total_new /= 100.0;
    total_ref /= 100.0;
	test_value(total_dst, total_ref, 1.0e-3, "Test operator+=(GSkyMap)");
	test_value(total_new, total_ref, 1.0e-3, "Test operator+(GSkyMap)");

    // Subtract pixels from destination map
    map_new  = map_dst - map_src;
    map_dst -= map_src;

    // Check total in destination map. Note that the total in the destination
    // map should be zero.
    total_dst = 0.0;
    total_new = 0.0;
    for (int pix = 0; pix < map_dst.npix(); ++pix) {
        for (int k = 0; k < map_dst.nmaps(); ++k) {
            total_dst += map_dst(pix,k);
            total_new += map_dst(pix,k);
        }
    }
	test_value(total_dst, 0.0, 1.0e-3, "Test operator-=(GSkyMap)");
	test_value(total_new, 0.0, 1.0e-3, "Test operator-(GSkyMap)");

    // Check map multiplication
    GSkyMap test_map  = map_src;
    test_map         *= map_src;
    map_new           = map_src * map_src;
    double total_test = 0.0;
    total_ref         = 0.0;
    total_new         = 0.0;
    for (int pix = 0; pix < test_map.npix(); ++pix) {
        for (int k = 0; k < test_map.nmaps(); ++k) {
            total_test += test_map(pix,k);
            total_new  += map_new(pix,k);
            total_ref  += map_src(pix,k) * map_src(pix,k);
        }
    }
	test_value(total_test, total_ref, 1.0e-3, "Test operator*=(GSkyMap)");
	test_value(total_new,  total_ref, 1.0e-3, "Test operator*(GSkyMap)");

    // Check map division
    test_map   = map_src;
    test_map  /= map_src;
    map_new    = map_src / map_src;
    total_test = 0.0;
    total_ref  = 0.0;
    total_new  = 0.0;
    for (int pix = 0; pix < test_map.npix(); ++pix) {
        for (int k = 0; k < test_map.nmaps(); ++k) {
            total_test += test_map(pix,k);
            total_new  += map_new(pix,k);
            total_ref  += map_src(pix,k) / map_src(pix,k);
        }
    }
	test_value(total_test, total_ref, 1.0e-3, "Test operator/=(GSkyMap)");
	test_value(total_new,  total_ref, 1.0e-3, "Test operator/(GSkyMap)");

    // Check map scaling
    test_map   = map_src;
    test_map  *= 3.3;
    total_test = 0.0;
    for (int pix = 0; pix < test_map.npix(); ++pix) {
        for (int k = 0; k < test_map.nmaps(); ++k) {
            total_test += test_map(pix,k);
        }
    }
	test_value(total_test, total_src*3.3, 1.0e-3, "Test operator*=(double)");

    // Check map division
    test_map   = map_src;
    test_map  /= 3.3;
    total_test = 0.0;
    for (int pix = 0; pix < test_map.npix(); ++pix) {
        for (int k = 0; k < test_map.nmaps(); ++k) {
            total_test += test_map(pix,k);
        }
    }
	test_value(total_test, total_src/3.3, 1.0e-3, "Test operator/=(double)");

    // Check map value addition
    test_map   = map_src;
    test_map  += 3.3;
    total_test = 0.0;
    double ref = total_src + 3.3*test_map.npix()*test_map.nmaps();
    for (int pix = 0; pix < test_map.npix(); ++pix) {
        for (int k = 0; k < test_map.nmaps(); ++k) {
            total_test += test_map(pix,k);
        }
    }
	test_value(total_test, ref, 1.0e-3, "Test operator+=(double)");

    // Check map value subtraction
    test_map   = map_src;
    test_map  -= 3.3;
    total_test = 0.0;
    ref        = total_src - 3.3*test_map.npix()*test_map.nmaps();
    for (int pix = 0; pix < test_map.npix(); ++pix) {
        for (int k = 0; k < test_map.nmaps(); ++k) {
            total_test += test_map(pix,k);
        }
    }
	test_value(total_test, ref, 1.0e-3, "Test operator-=(double)");

    // Save maps
    map_src.save("test_map_src.fits", true);
    map_dst.save("test_map_dst.fits", true);

    // Test map stacking
    GSkyMap map_stacked = map_src;
    map_stacked.stack_maps();
    double total_stacked = 0.0;
    for (int pix = 0; pix < map_stacked.npix(); ++pix) {
        total_stacked += map_stacked(pix);
    }
	test_value(total_stacked, total_src, 1.0e-3, "Test stack_maps() method");
	test_value(map_stacked.nmaps(), 1, "Test stack_maps() method");

    // Test total counts computation
    GNdarray counts_spectrum = map_src.counts();
    double total_counts = 0.0;
    for (int i = 0; i < counts_spectrum.size(); ++i) {
        total_counts += counts_spectrum(i);
    }
	test_value(total_counts, total_src, 1.0e-3, "Test counts() method");
	test_value(counts_spectrum.size(), map_src.nmaps(), "Test counts() method");

    // Test total flux computation
    GNdarray flux_spectrum = map_src.flux();
	test_value(flux_spectrum.size(), map_src.nmaps(), "Test flux() method");

    // Test map number changing
    GSkyMap map_more = map_src;
    map_more.nmaps(4);
    double total_more = 0.0;
    for (int k = 0; k < map_more.nmaps(); ++k) {
        for (int pix = 0; pix < map_more.npix(); ++pix) {
            total_more += map_more(pix,k);
        }
    }
	test_value(total_more, total_src, 1.0e-3, "Test nmaps() method with more maps");
	test_value(map_more.nmaps(), 4, "Test nmaps() method with more maps");
    GSkyMap map_less = map_src;
    map_less.nmaps(1);
    double total_less = 0.0;
    for (int k = 0; k < map_less.nmaps(); ++k) {
        for (int pix = 0; pix < map_less.npix(); ++pix) {
            total_less += map_less(pix,k);
        }
    }
	test_value(total_less, 0.5*total_src, 1.0e-3, "Test nmaps() method with less maps");
	test_value(map_less.nmaps(), 1, "Test nmaps() method with less maps");

    // Test map extraction
    GSkyMap map_extract = map_src.extract(0);
    double total_extract = 0.0;
    for (int k = 0; k < map_extract.nmaps(); ++k) {
        for (int pix = 0; pix < map_extract.npix(); ++pix) {
            total_extract += map_extract(pix,k);
        }
    }
	test_value(total_extract, 0.5*total_src, 1.0e-3, "Test extract() method with 1 map");
	test_value(map_extract.nmaps(), 1, "Test extract() method with 1 map");
    map_extract = map_src.extract(0,2);
    total_extract = 0.0;
    for (int k = 0; k < map_extract.nmaps(); ++k) {
        for (int pix = 0; pix < map_extract.npix(); ++pix) {
            total_extract += map_extract(pix,k);
        }
    }
	test_value(total_extract, total_src, 1.0e-3, "Test extract() method with 2 maps");
	test_value(map_extract.nmaps(), 2, "Test extract() method with 2 maps");

    // Define one more map for shaping manipulation
    GSkyMap map_shape0("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 10, 10);

    // Test map dimensions and shape
	test_value(map_shape0.nmaps(), 1, "Check that sky map contains one map");
	test_value(map_shape0.ndim(), 1, "Check that sky map has one dimension");
	test_value(map_shape0.shape().size(), 1, "Check that sky map has a shape size of 1");
	test_value(map_shape0.shape()[0], 1, "Check that sky map has one map");

    // Save map
    map_shape0.save("test_map_shape0.fits", true);

    // Load map
    GSkyMap map_load_shape0("test_map_shape0.fits");

    // Test map dimensions and shape
	test_value(map_load_shape0.nmaps(), 1, "Check that sky map contains one map");
	test_value(map_load_shape0.ndim(), 1, "Check that sky map has one dimension");
	test_value(map_load_shape0.shape().size(), 1, "Check that sky map has a shape size of 1");
	test_value(map_load_shape0.shape()[0], 1, "Check that sky map has one map");

    // Define one more map for shaping manipulation
    GSkyMap map_shape("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 10, 10, 12);

    // Test initial map dimensions and shape
	test_value(map_shape.nmaps(), 12, "Check that sky map contains 12 maps");
	test_value(map_shape.ndim(), 1, "Check that sky map has one dimension");
	test_value(map_shape.shape().size(), 1, "Check that sky map has a shape size of 1");
	test_value(map_shape.shape()[0], 12, "Check that sky map has 12 maps");

    // Save map
    map_shape.save("test_map_shape1.fits", true);

    // Load map
    GSkyMap map_load_shape1("test_map_shape1.fits");

    // Test loaded map dimensions and shape
	test_value(map_load_shape1.nmaps(), 12, "Check that sky map contains 12 maps");
	test_value(map_load_shape1.ndim(), 1, "Check that sky map has one dimension");
	test_value(map_load_shape1.shape().size(), 1, "Check that sky map has a shape size of 1");
	test_value(map_load_shape1.shape()[0], 12, "Check that sky map has 12 maps");

    // Set new map shape
    map_shape.shape(3,4);

    // Test map dimensions and shape
	test_value(map_shape.nmaps(), 12, "Check that sky map contains 12 maps");
	test_value(map_shape.ndim(), 2, "Check that sky map has two dimensions");
	test_value(map_shape.shape().size(), 2, "Check that sky map has a shape size of 2");
	test_value(map_shape.shape()[0], 3, "Check that sky map has 3 maps in first dimension");
	test_value(map_shape.shape()[1], 4, "Check that sky map has 4 maps in second dimension");

    // Save map
    map_shape.save("test_map_shape2.fits", true);

    // Load map
    GSkyMap map_load_shape2("test_map_shape2.fits");

    // Test map dimensions and shape
	test_value(map_load_shape2.nmaps(), 12, "Check that sky map contains 12 maps");
	test_value(map_load_shape2.ndim(), 2, "Check that sky map has two dimensions");
	test_value(map_load_shape2.shape().size(), 2, "Check that sky map has a shape size of 2");
	test_value(map_load_shape2.shape()[0], 3, "Check that sky map has 3 maps in first dimension");
	test_value(map_load_shape2.shape()[1], 4, "Check that sky map has 4 maps in second dimension");

    // Set new map shape
    map_shape.shape(2,3,2);

    // Test map dimensions and shape
	test_value(map_shape.nmaps(), 12, "Check that sky map contains 12 maps");
	test_value(map_shape.ndim(), 3, "Check that sky map has three dimensions");
	test_value(map_shape.shape().size(), 3, "Check that sky map has a shape size of 3");
	test_value(map_shape.shape()[0], 2, "Check that sky map has 2 maps in first dimension");
	test_value(map_shape.shape()[1], 3, "Check that sky map has 3 maps in second dimension");
	test_value(map_shape.shape()[2], 2, "Check that sky map has 2 maps in third dimension");

    // Save map
    map_shape.save("test_map_shape3.fits", true);

    // Load map
    GSkyMap map_load_shape3("test_map_shape3.fits");

    // Test map dimensions and shape
	test_value(map_load_shape3.nmaps(), 12, "Check that sky map contains 12 maps");
	test_value(map_load_shape3.ndim(), 3, "Check that sky map has three dimensions");
	test_value(map_load_shape3.shape().size(), 3, "Check that sky map has a shape size of 3");
	test_value(map_load_shape3.shape()[0], 2, "Check that sky map has 2 maps in first dimension");
	test_value(map_load_shape3.shape()[1], 3, "Check that sky map has 3 maps in second dimension");
	test_value(map_load_shape3.shape()[2], 2, "Check that sky map has 2 maps in third dimension");

    // Load map for smoothing
    GSkyMap map_smooth(sky_map);

    // Compute sum of sky map for reference
    double ref_smooth = 0.0;
    for (int pix = 0; pix < map_smooth.npix(); ++pix) {
        ref_smooth += map_smooth(pix);
    }

    // Smooth map using DISK kernel
    GSkyMap map_smooth1 = map_smooth;
    map_smooth1.smooth("DISK", 0.2);
    map_smooth1.save("test_map_smooth_disk.fits", true);

    // Compute sum of smoothed sky map
    double sum_smooth1 = 0.0;
    for (int pix = 0; pix < map_smooth1.npix(); ++pix) {
        sum_smooth1 += map_smooth1(pix);
    }
	test_value(sum_smooth1, ref_smooth, "Check DISK smoothing");

    // Smooth map using GAUSSIAN kernel
    GSkyMap map_smooth2 = map_smooth;
    map_smooth2.smooth("GAUSSIAN", 0.2);
    map_smooth2.save("test_map_smooth_gaussian.fits", true);

    // Compute sum of smoothed sky map
    double sum_smooth2 = 0.0;
    for (int pix = 0; pix < map_smooth2.npix(); ++pix) {
        sum_smooth2 += map_smooth2(pix);
    }
	test_value(sum_smooth2, ref_smooth, "Check GAUSSIAN smoothing");

    // Trim the skymap based on pixels
    int startx = 2;
    int stopx  = 5;
    int starty = 2;
    int stopy  = 5;
    test_map   = map_src.extract(startx, stopx, starty, stopy);
    test_value(test_map.npix(), 16, "Check trimmed map pixel count (index)");
    test_value(test_map.nmaps(), map_src.nmaps(), "Check total maps in trimmed map (index)");

    // Check the pixel values
    total_test = 0.0;
    total_src  = 0.0;
    for (int inx = 0; inx < test_map.npix(); ++inx) {
        // Update pixel for source map
        GSkyPixel pix = test_map.inx2pix(inx);
        pix.x(startx + pix.x());
        pix.y(starty + pix.y());

        // Sum over maps
        for (int k = 0; k < test_map.nmaps(); ++k) {
            total_test += test_map(inx,k);
            total_src  += map_src(pix,k);
        }
    }
    test_value(total_test, total_src, "Check Trimmed Map sum (index)");

    // Create a regions map
    GSkyRegionMap exclmap(test_map.extract(1));
    GSkyRegions inclusions;
    inclusions.append(exclmap);

    // Trim the skymap based on the inclusion regions
    test_map = map_src.extract(inclusions);
    test_value(test_map.npix(), 36, "Check trimmed map pixel count (region)");
    test_value(test_map.nmaps(), map_src.nmaps(), "Check total maps in trimmed map (region)");

    // Check the pixel values
    total_test = 0.0;
    total_src  = 0.0;
    for (int inx = 0; inx < test_map.npix(); ++inx) {
        // Update pixel for source map
        GSkyPixel pix = test_map.inx2pix(inx);
        pix.x(startx - 1 + pix.x());
        pix.y(starty - 1 + pix.y());

        // Sum over maps
        for (int k = 0; k < test_map.nmaps(); ++k) {
            total_test += test_map(inx,k);
            total_src  += map_src(pix,k);
        }
    }
    test_value(total_test, total_src, "Check Trimmed Map sum (region)");

    // Define maps for operator testing
    GSkyMap map_op1("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 10, 10, 2);
    GSkyMap map_op2("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 10, 10, 2);
    GSkyMap map_op3("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 12, 12, 2);
    for (int layer = 0; layer < map_op1.nmaps(); ++layer) {
        for (int inx = 0; inx < map_op1.npix(); ++inx) {
            map_op1(inx,layer) = double(inx+map_op1.npix()*layer);
        }
    }
    for (int layer = 0; layer < map_op2.nmaps(); ++layer) {
        for (int inx = 0; inx < map_op2.npix(); ++inx) {
            map_op2(inx,layer) = double(inx+map_op2.npix()*layer);
        }
    }
    for (int layer = 0; layer < map_op3.nmaps(); ++layer) {
        for (int inx = 0; inx < map_op3.npix(); ++inx) {
            map_op3(inx,layer) = double(inx+map_op3.npix()*layer);
        }
    }

    // Test addition
    test_map = map_op1 + map_op2;
    test_value(test_map(45,1), 290.0, "Check addition of identical maps");
    test_map = map_op1 + map_op3;
    test_value(test_map(45,1), 355.0, "Check addition of different maps");

    // Test subtraction
    test_map = map_op1 - map_op2;
    test_value(test_map(45,1), 0.0, "Check subtraction of identical maps");
    test_map = map_op1 - map_op3;
    test_value(test_map(45,1), -65.0, "Check subtraction of different maps");

    // Test multiplication
    test_map = map_op1 * map_op2;
    test_value(test_map(45,1), 21025.0, "Check multiplication of identical maps");
    test_map = map_op1 * map_op3;
    test_value(test_map(45,1), 30450.0, "Check multiplication of different maps");


    // Test division
    test_map = map_op1 / map_op2;
    test_value(test_map(45,1), 1.0, "Check division of identical maps");
    test_map = map_op1 / map_op3;
    test_value(test_map(45,1), 0.6904761905, "Check division of different maps");

    // Exit test
    return;
}


/***************************************************************************
 * @brief GSkyMap
 *
 * This method tests the GSkyMap class input and output. Specifically, the
 * test verifies that saving of a HEALPix and WCS map in the same FITS file
 * is possible, and that upon opening the correct extension is selected.
 ***************************************************************************/
void TestGSky::test_GSkyMap_io(void)
{
    // Set filename
    GFilename filename("test_skymap.fits");
    filename.remove();

    // Create HEALPix and WCS maps and save them to FITS file
    GSkyMap healpix1("GAL", 4, "RING", 1);
    GSkyMap healpix2("GAL", 2, "RING", 1);
    GSkyMap wcs1("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 10, 10);
    GSkyMap wcs2("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 5, 5);
    GSkyMap wcs3("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 5, 5, 12);
    for (int i = 0; i < healpix1.npix(); ++i) {
        healpix1(i) = double(i+1);
    }
    for (int i = 0; i < healpix2.npix(); ++i) {
        healpix2(i) = 3.0*double(i+1);
    }
    for (int i = 0; i < wcs1.npix(); ++i) {
        wcs1(i) = double(i+1);
    }
    for (int i = 0; i < wcs2.npix(); ++i) {
        wcs2(i) = 3.0*double(i+1);
    }
    for (int k = 0; k < 12; ++k) {
        for (int i = 0; i < wcs3.npix(); ++i) {
            wcs3(i,k) = 3.0*double(i+k+1);
        }
    }
    GFits fits(filename, true);
    healpix1.write(fits, "HEALPIX1");
    wcs1.write(fits, "WCS1");
    healpix2.write(fits, "HEALPIX2");
    wcs2.write(fits, "WCS2");
    wcs3.shape(2,2,3);
    wcs3.write(fits, "WCS3");
    fits.save(true);

    // Load HEALPix map and check content
    GSkyMap map1(filename+"[HEALPIX1]");
	test_value(map1.npix(), healpix1.npix(), "Check number of HEALPix pixels");
    int failures1(0);
    for (int i = 0; i < map1.npix(); ++i) {
        if (std::abs(map1(i)-double(i+1)) > 1.0e-6) {
            failures1++;
        }
    }
	test_value(failures1, 0, "Check number of invalid HEALPix pixel values");

    // Load WCS1 map and check content
    GSkyMap map2(filename+"[WCS1]");
	test_value(map2.npix(), wcs1.npix(), "Check number of WCS pixels");
    int failures2(0);
    for (int i = 0; i < map2.npix(); ++i) {
        if (std::abs(map2(i)-double(i+1)) > 1.0e-6) {
            failures2++;
        }
    }
	test_value(failures2, 0, "Check number of invalid WCS pixel values");

    // Load WCS3 map and check content
    GSkyMap map3(filename+"[WCS3]");
	test_value(map3.npix(),  wcs3.npix(), "Check number of WCS pixels");
	test_value(map3.nmaps(), wcs3.nmaps(), "Check number of WCS maps");
	test_value(map3.ndim(),  wcs3.ndim(), "Check dimension of WCS maps");
    int failures3(0);
    for (int k = 0; k < 12; ++k) {
        for (int i = 0; i < map3.npix(); ++i) {
            if (std::abs(map3(i,k)-3.0*double(i+k+1)) > 1.0e-6) {
                failures3++;
            }
        }
    }
	test_value(failures3, 0, "Check number of invalid WCS pixel values");

    // Exit test
    return;
}


/***************************************************************************
 * @brief Test GSkyRegions class
 ***************************************************************************/
 void TestGSky::test_GSkyRegions(void)
{
    // Allocate regions
    GSkyRegions regions;

    // Test regions loading
    test_try("Test regions loading");
    try {
        regions.load(sky_region);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Check if region was loaded correctly
    test_value(regions.size(), 5, "5 sky regions expected");
    GSkyRegionCircle* circle = dynamic_cast<GSkyRegionCircle*>(regions[0]);
    test_assert(circle->type() == "Circle", "Region is not a circle");
    test_value(circle->radius(), 10.0, 1.0e-10);
    test_value(circle->ra(), 0.1, 1.0e-10);
    test_value(circle->dec(), -35.6, 1.0e-10);
    test_value(regions[0]->name(), "", "No name for first sky region");
    test_value(regions[1]->name(), "", "Empty name for second sky region");
    test_value(regions[2]->name(), "Test 1", "Name \"Test 1\" for third sky region");
    test_value(regions[3]->name(), "Test 2", "Name \"Test 2\" for fourth sky region");
    test_value(regions[4]->name(), "Test 3", "Name \"Test 3\" for fifth sky region");

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
    test_value(regions.size(), 5, "5 sky regions expected");
    circle = dynamic_cast<GSkyRegionCircle*>(regions[0]);
    test_assert(circle->type() == "Circle", "Region is not a circle");
    test_value(circle->radius(), 10.0, 1.0e-10);
    test_value(circle->ra(), 0.1, 1.0e-10);
    test_value(circle->dec(), -35.6, 1.0e-10);
    test_value(regions[0]->name(), "", "No name for first sky region");
    test_value(regions[1]->name(), "", "Empty name for second sky region");
    test_value(regions[2]->name(), "Test 1", "Name \"Test 1\" for third sky region");
    test_value(regions[3]->name(), "Test 2", "Name \"Test 2\" for fourth sky region");
    test_value(regions[4]->name(), "Test 3", "Name \"Test 3\" for fifth sky region");

    // Exit test
    return;
}


/***************************************************************************
 * @brief Test GSkyRegionCircle
 ***************************************************************************/
void TestGSky::test_GSkyRegionCircle(void)
{
	// Define region for comparison
	GSkyDir refdir_radeczerozero;
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

    // Test constructing with radius -1
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

    // Set reference sky directions
    //GSkyDir refdir_radeczerozero;
    GSkyDir refdir_lbzerozero;
    GSkyDir refdir_outside_refregion;
    GSkyDir refdir_raoffset;
    GSkyDir refdir_decpole;
    GSkyDir refdir_ndecpole;
    GSkyDir refdir_rapole;
    refdir_radeczerozero.radec_deg(0,0);
    refdir_lbzerozero.lb_deg(0,0);
    refdir_outside_refregion.radec_deg(100,80);
    refdir_raoffset.radec_deg(5,0);
    refdir_decpole.radec_deg(0,89);
	refdir_ndecpole.radec_deg(100,89);
    refdir_rapole.radec_deg(359,0);

    // Set reference sky regions
    refregion = GSkyRegionCircle(refdir_radeczerozero, 10.0);
    GSkyDir          centre = refregion.centre();
    GSkyRegionCircle refregion_smaller(refdir_radeczerozero, 5.0);
    GSkyRegionCircle refregion_larger(refdir_radeczerozero, 20.0);
    GSkyRegionCircle refregion_raoffset(refdir_raoffset, 10.0);
    GSkyRegionCircle refregion_rapole(refdir_rapole, 3.0);
    GSkyRegionCircle refregion_decpole(refdir_decpole, 3.0);

    // Test contain dirs
	test_assert(refregion.contains(refdir_radeczerozero),"Test for containment");
    test_assert(refregion.contains(centre), "Test if centre in circle");
	test_assert(!refregion.contains(refdir_outside_refregion), "Test2 for containment");

	// Test contain regions
	test_assert(refregion.contains(refregion_smaller),"Test for containment region");
	test_assert(!refregion.contains(refregion_larger), "Test for containment region2 ");
	test_assert(refregion.contains(refregion), "Test3 for containment region");

	test_assert(refregion.contains(refdir_rapole), "Test rapole for containment region");
	test_assert(refregion_decpole.contains(refdir_ndecpole), "Test rapole for containment region");

    // Test overlaps
	test_assert(refregion.overlaps(refregion_smaller),"Test for overlap");
	test_assert(refregion.overlaps(refregion_larger),"Test2 for overlap");
	test_assert(refregion.overlaps(refregion_raoffset),"Test3 for overlap");
	test_assert(!refregion.overlaps(refregion_decpole),"Test4 for overlap");

    // Test equality and non-equality operators
	test_assert(refregion_smaller == refregion_smaller,"Test for equality");
	test_assert(!(refregion_smaller == refregion_larger),"Test for equality");
	test_assert(!(refregion_smaller != refregion_smaller),"Test for non-equality");
	test_assert(refregion_smaller != refregion_larger,"Test for non-equality");

	// Exit test
    return;
}


/***************************************************************************
 * @brief GSkyRegionRect
 ***************************************************************************/
void TestGSky::test_GSkyRegionRect(void)
{
    // Define region for comparison
    GSkyDir refdir_radeczerozero = GSkyDir();
    refdir_radeczerozero.radec_deg(0,0);
    double ref_width  = 1.;
    double ref_height = 1.;
    double ref_pa_deg = 0.;

    // Test constructing
    test_try("Test constructor");
    try {
        GSkyRegionRect rect(refdir_radeczerozero, ref_width, ref_height, ref_pa_deg);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test constructing with negative width (should be forbidden)
    test_try("Test constructor2");
    try {
        GSkyRegionRect rect(refdir_radeczerozero, -1, ref_height, ref_pa_deg);
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test constructing with negative height (should be forbidden)
    test_try("Test constructor3");
    try {
        GSkyRegionRect rect(refdir_radeczerozero, ref_width, -1, ref_pa_deg);
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test constructing with negative posang (should be allowed)
    test_try("Test constructor4");
    try {
        GSkyRegionRect rect(refdir_radeczerozero, ref_width, ref_height, -1);
        test_try_success();
    }
    catch (GException::invalid_argument &e) {
        test_try_failure();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test assignment of neg width after construction
    test_try("Test width assignment after");
    try {
        GSkyRegionRect rect(refdir_radeczerozero, ref_width, ref_height, ref_pa_deg);
        rect.width(-1.0);
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test assignment of neg height after construction
    test_try("Test height assignment after");
    try {
        GSkyRegionRect rect(refdir_radeczerozero, ref_width, ref_height, ref_pa_deg);
        rect.height(-1.0);
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test assignment of neg posang after construction
    test_try("Test posang assignment after");
    try {
        GSkyRegionRect rect(refdir_radeczerozero, ref_width, ref_height, ref_pa_deg);
        rect.posang(-1.0);
        test_try_success();
    }
    catch (GException::invalid_argument &e) {
        test_try_failure();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Check radius assignment
    GSkyRegionRect rect(refdir_radeczerozero, ref_width, ref_height, ref_pa_deg);

    double check_width  = rect.width();
    test_value(ref_width, check_width, 1.0e-10, "Test width assignment");

    double check_height = rect.width();
    test_value(ref_height, check_height, 1.0e-10, "Test height assignment");

    double check_pa_deg = rect.posang();
    test_value(ref_pa_deg, check_pa_deg, 1.0e-10, "Test posang assignment");

    // Check solid angle assignment
    double solidangle_check = rect.solidangle();
    double solidangle = ref_width*gammalib::deg2rad * ref_height*gammalib::deg2rad;
    test_value(solidangle_check, solidangle, 1.0e-10, "Test solid angle");

    // Initialise objects
    refdir_radeczerozero = GSkyDir();
    GSkyDir test_dir_0 = GSkyDir();
    GSkyDir test_dir_1 = GSkyDir();
    GSkyDir test_dir_2 = GSkyDir();

    refdir_radeczerozero.radec_deg(0,0);
    GSkyRegionRect refregion(refdir_radeczerozero, 10, 6, 0);

    // Check GSkyDir containment
    // Upper right
    test_dir_0.radec_deg(6, 2);
    test_dir_1.radec_deg(4, 2);
    test_dir_2.radec_deg(4, 4);
    test_assert(!refregion.contains(test_dir_0),"test for containment");
    test_assert(refregion.contains(test_dir_1),"test for containment");
    test_assert(!refregion.contains(test_dir_2),"test for containment");

    // Upper left
    test_dir_0.radec_deg(-6, 2);
    test_dir_1.radec_deg(-4, 2);
    test_dir_2.radec_deg(-4, 4);
    test_assert(!refregion.contains(test_dir_0),"test for containment");
    test_assert(refregion.contains(test_dir_1),"test for containment");
    test_assert(!refregion.contains(test_dir_2),"test for containment");

    // Lower left
    test_dir_0.radec_deg(-6, -2);
    test_dir_1.radec_deg(-4, -2);
    test_dir_2.radec_deg(-4, -4);
    test_assert(!refregion.contains(test_dir_0),"test for containment");
    test_assert(refregion.contains(test_dir_1),"test for containment");
    test_assert(!refregion.contains(test_dir_2),"test for containment");

    // Lower right
    test_dir_0.radec_deg(6, -2);
    test_dir_1.radec_deg(4, -2);
    test_dir_2.radec_deg(4, -4);
    test_assert(!refregion.contains(test_dir_0),"test for containment");
    test_assert(refregion.contains(test_dir_1),"test for containment");
    test_assert(!refregion.contains(test_dir_2),"test for containment");

    // Corners
    test_assert(refregion.contains(refregion.get_corner(0)),"test for corner containment");
    test_assert(refregion.contains(refregion.get_corner(1)),"test for corner containment");
    test_assert(refregion.contains(refregion.get_corner(2)),"test for corner containment");
    test_assert(refregion.contains(refregion.get_corner(3)),"test for corner containment");


    // Check GSkyDir containment for rotated rectangle
    GSkyRegionRect refregion_rot(refdir_radeczerozero, 10, 10, 45);

    // Check corner containment
    test_assert(refregion_rot.contains(refregion_rot.get_corner(0)),"test for corner containment");
    test_assert(refregion_rot.contains(refregion_rot.get_corner(1)),"test for corner containment");
    test_assert(refregion_rot.contains(refregion_rot.get_corner(2)),"test for corner containment");
    test_assert(refregion_rot.contains(refregion_rot.get_corner(3)),"test for corner containment");

    // Check GSkyDir containment
    test_dir_0.radec_deg(refregion_rot.get_corner(0).ra_deg()+1, 0);
    test_dir_1.radec_deg(refregion_rot.get_corner(0).ra_deg()  , 0);
    test_dir_2.radec_deg(refregion_rot.get_corner(0).ra_deg()-1, 0);
    test_assert(!refregion_rot.contains(test_dir_0),"test for containment");
    test_assert(refregion_rot.contains(test_dir_1),"test for containment");
    test_assert(refregion_rot.contains(test_dir_2),"test for containment");

    // Check local<->global coordinate transformations
    test_dir_0.radec_deg(50,50);
    test_dir_1 = refregion_rot.transform_to_local(test_dir_0);
    test_dir_2 = refregion_rot.transform_to_global(test_dir_1);
    test_value(test_dir_0.ra(),  test_dir_2.ra(),  1.0e-10, "Test coord trafo");
    test_value(test_dir_0.dec(), test_dir_2.dec(), 1.0e-10, "Test coord trafo");

    // Check other rects containment
    GSkyDir refdir_raoffset = GSkyDir();
    refdir_raoffset.radec_deg(5,0);

    GSkyDir refdir_rapole = GSkyDir();
    refdir_rapole.radec_deg(359,0);

    GSkyDir refdir_decpole = GSkyDir();
    refdir_decpole.radec_deg(0,89);

    GSkyDir refdir_ndecpole = GSkyDir();
    refdir_ndecpole.radec_deg(100,89);

    GSkyRegionRect refregion_smaller(refdir_radeczerozero, 5,  3, 0);
    GSkyRegionRect refregion_larger(refdir_radeczerozero, 20, 12, 0);
    GSkyRegionRect refregion_raoffset(refdir_raoffset, 10, 6, 0);
    GSkyRegionRect refregion_rapole(refdir_rapole,   3, 3, 0);
    GSkyRegionRect refregion_decpole(refdir_decpole, 3, 3, 0);

    // Test contain dirs (same from circle)
    test_assert(refregion.contains(refdir_radeczerozero),"test for containment");
    test_assert(refregion.contains(refregion.centre()), "test if centre in rectangle");

    // Test contain regions (same from circle)
    test_assert(refregion.contains(refregion_smaller),"test for containment region");
    test_assert(!refregion.contains(refregion_larger), "test for containment region2 ");
    test_assert(refregion.contains(refregion), "test3 for containment region");

    test_assert(refregion.contains(refdir_rapole), "rapole for containment region");
    test_assert(refregion_decpole.contains(refdir_ndecpole), "rapole for containment region");

    // Test overlaps (same from circle)
    test_assert(refregion.overlaps(refregion_smaller),"test for overlap");
    test_assert(refregion.overlaps(refregion_larger),"test2 for overlap");
    test_assert(refregion.overlaps(refregion_raoffset),"test3 for overlap");
    test_assert(!refregion.overlaps(refregion_decpole),"test4 for overlap");

    // Exit test
    return;
}


/***************************************************************************
 * @brief GSkyRegionMap_construct
 ***************************************************************************/
void TestGSky::test_GSkyRegionMap(void)
{
    // Test file constructor:
    test_try("Test file constructor");
    try {
        GSkyRegionMap regmap(sky_region_map);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test skymap constructor
    GSkyMap skymap("TAN", "CEL", 0.0, 0.0, 0.1, 0.1, 100, 100);
    test_try("Test skymap constructor");
    try {
      GSkyRegionMap regmap(skymap);
      test_try_success();
    }
    catch (std::exception &e) {
      test_try_failure(e);
    }

    // Load reference region map (a mask centred on (RA,Dec)=(0,0) with radius 0.3deg)
    GSkyRegionMap regmap(sky_region_map);

    // Set reference directions and regions
    GSkyDir refdir_radeczerozero = GSkyDir();
    GSkyDir refdir_radeclargeshift = GSkyDir();
    GSkyDir refdir_radecsmallshift = GSkyDir();
    refdir_radeczerozero.radec_deg(0.0,0.0);
    refdir_radeclargeshift.radec_deg(2.0,0.0);
    refdir_radecsmallshift.radec_deg(0.3,0.0);

    // Test contains dirs
    test_assert(regmap.contains(refdir_radeczerozero),"test for direction containment");
    test_assert(!regmap.contains(refdir_radeclargeshift),"test2 for direction containment");


    // Test circles
    GSkyRegionCircle inregion(refdir_radeczerozero,0.2);
    GSkyRegionCircle outregion(refdir_radeclargeshift,0.2);
    GSkyRegionCircle overregion(refdir_radecsmallshift,0.3);

    // Test contains regions
    test_assert(regmap.contains(inregion),"test for region containment");
    test_assert(!regmap.contains(outregion),"test2 for region containment");

    // Test overlaps regions
    test_assert(regmap.overlaps(inregion),"test for region overlap");
    test_assert(regmap.overlaps(overregion),"test2 for region overlap");


    // Test rectangles
    GSkyRegionRect inrect(refdir_radeczerozero, 0.3, 0.3, 0);
    GSkyRegionRect outrect(refdir_radeclargeshift, 0.3, 0.3, 0);
    GSkyRegionRect overrect(refdir_radecsmallshift, 0.3, 0.3, 0);

    // Test contains regions
    test_assert(regmap.contains(inrect),"test for region containment");
    test_assert(!regmap.contains(outrect),"test2 for region containment");

    // Test overlaps regions
    test_assert(regmap.overlaps(inrect),"test for region overlap");
    test_assert(regmap.overlaps(overrect),"test2 for region overlap");


    // Exit test
    return;
}


/***************************************************************************
 * @brief Test GHorizDir class
 ***************************************************************************/
void TestGSky::test_GHorizDir(void)
{

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

    // Create a test suite
    TestGSky test;

    // Append to the container
    testsuites.append(test);

    // Run
    bool was_successful = testsuites.run();

    // Save xml report
    testsuites.save("reports/GSky.xml");

    // Return
    return was_successful ? 0:1;
}
