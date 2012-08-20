/***************************************************************************
 *                  test_GSky.cpp  -  test GSky classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include "GammaLib.hpp"
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_WCS_FORTH_BACK_PIXEL_DEBUG
#define G_WCS_COPY_DEBUG


/***********************************************************************//**
 * @brief Test consistency of forward and background transformations
 *
 * @param[in] nx Number of points in X
 * @param[in] ny Number of points in Y
 * @param[in] crpix1 Reference pixel in X
 * @param[in] crpix2 Reference pixel in Y
 ***************************************************************************/
double wcs_forth_back_pixel(GWcslib* wcs, int nx, int ny, double& crpix1, double& crpix2)
{
    // Initialise maximal distance
    double dist_max = 0.0;
    
    // Loop over x pixels
    for (int ix = -nx; ix <= nx; ++ix) {
    
        // Set x value
        double x = double(ix) + crpix1;
        
        // Skip pixels outside valid range (they are not expected to
        // transform bijectively)
        if (x < 0 || x >= nx)
            continue;
        
        // Loop over y pixels
        for (int iy = -ny; iy <= ny; ++iy) {
        
            // Set y value
            double    y = double(iy) + crpix2;
            
            // Skip pixels outside valid range (they are not expected to
            // transform bijectively)
            if (y < 0 || y >= ny)
                continue;
            
            // Set sky pixel
            GSkyPixel pix_in(x,y);
            
            // Forth: Transform to world
            GSkyDir dir = wcs->xy2dir(pix_in);
            
            // Back: Transform to pixel
            GSkyPixel pix_out = wcs->dir2xy(dir);
            
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
            if (dist > dist_max)
                dist_max = dist;
                
        } // endfor: y pixels
    } // endfor: x pixels

    // Return
    return dist_max;
}


/***********************************************************************//**
 * @brief Test consistency of copy
 *
 * @param[in] nx Number of points in X
 * @param[in] ny Number of points in Y
 * @param[in] crpix1 Reference pixel in X
 * @param[in] crpix2 Reference pixel in Y
 ***************************************************************************/
double wcs_copy(GWcslib* wcs, int nx, int ny, double& crpix1, double& crpix2)
{
    // Make a copy using clone
    GWcslib* cpy = wcs->clone();
    
    // Initialise maximal angle and distance
    double angle_max = 0.0;
    double dist_max  = 0.0;
    
    // Loop over x pixels
    for (int ix = -nx; ix <= nx; ++ix) {
    
        // Set x value
        double x = double(ix) + crpix1;
        
        // Skip pixels outside valid range (they are not expected to
        // transform bijectively)
        if (x < 0 || x >= nx)
            continue;
        
        // Loop over y pixels
        for (int iy = -ny; iy <= ny; ++iy) {
        
            // Set y value
            double    y = double(iy) + crpix2;
            
            // Skip pixels outside valid range (they are not expected to
            // transform bijectively)
            if (y < 0 || y >= ny)
                continue;
            
            // Set sky pixel
            GSkyPixel pix_in(x,y);
            
            // Transform to world
            GSkyDir dir1 = wcs->xy2dir(pix_in);
            GSkyDir dir2 = cpy->xy2dir(pix_in);
            double angle = dir1.dist_deg(dir2);
            if (angle > angle_max)
                angle_max = angle;

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
            GSkyPixel pix_out1 = wcs->dir2xy(dir1);
            GSkyPixel pix_out2 = cpy->dir2xy(dir2);
            
            // Compute distance
            double dx   = pix_out1.x()-pix_out2.x();
            double dy   = pix_out1.y()-pix_out2.y();
            double dist = sqrt(dx*dx+dy*dy);
            if (dist > dist_max)
                dist_max = dist;
            
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

    // Return
    return ((dist_max > angle_max) ? dist_max : angle_max);
}


/***********************************************************************//**
 * @brief Test GWcslib projections
 *
 * This function test all non-HealPix projections that have been registered.
 ***************************************************************************/
void test_GWcslib(void)
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

    // Main test protection
    try {
    
        // Loop over all non-HealPix projections
        GWcsRegistry registry;
        for (int i = 0; i < registry.size(); ++i) {
    
            // Skip HealPix
            if (registry.code(i) == "HPX")
                continue;
        
            // Dump header
            std::string test_class = "GWcs" + registry.code(i);
            std::cout << "Test "+test_class+": ";

            // Perform tests
            try {
        
                // Allocate projection CEL and GAL
                GWcslib *cel = (GWcslib*)registry.alloc(registry.code(i));
                GWcslib *gal = (GWcslib*)registry.alloc(registry.code(i));
                
                // Throw an error if allocation failed
                if (cel == NULL || gal == NULL) {
                    std::cout << std::endl 
                              << "TEST ERROR: Could not allocate "+test_class+"."
                              << std::endl;
                    throw;
                }
                
                // Set projections
                try {
                    cel->set("CEL", crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
                    gal->set("GAL", crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
                }
                catch (std::exception &e) {
                    std::cout << std::endl 
                              << "TEST ERROR: Error occured while setting up "+test_class+"."
                              << std::endl;
                    std::cout << e.what() << std::endl;
                    throw;
                }
                std::cout << ".";
                
                // Test forth-back CEL conversion
                try {
                    double tol = 0.0;
                    if ((tol = wcs_forth_back_pixel(cel, nx, ny, crpix1, crpix2)) > 1.0e-10) {
                        std::cout << std::endl 
                                  << "TEST ERROR: CEL forth-back transformation tolerance 1.0e-10 exceeded: "
                                  << tol << std::endl;
                        throw;
                    }
                }
                catch (std::exception &e) {
                    std::cout << std::endl 
                              << "TEST ERROR: Error while testing CEL forth-back transformations."
                              << std::endl;
                    std::cout << e.what() << std::endl;
                    throw;
                }
                std::cout << ".";

                // Test forth-back GAL conversion
                try {
                    double tol = 0.0;
                    if ((tol = wcs_forth_back_pixel(gal, nx, ny, crpix1, crpix2)) > 1.0e-10) {
                        std::cout << std::endl 
                                  << "TEST ERROR: GAL forth-back transformation tolerance 1.0e-10 exceeded: "
                                  << tol << std::endl;
                        throw;
                    }
                }
                catch (std::exception &e) {
                    std::cout << std::endl 
                              << "TEST ERROR: Error while testing GAL forth-back transformations."
                              << std::endl;
                    std::cout << e.what() << std::endl;
                    throw;
                }
                std::cout << ".";

                // Test CEL copy
                try {
                    double tol = 0.0;
                    if ((tol = wcs_copy(cel, nx, ny, crpix1, crpix2)) > 1.0e-4) {
                        std::cout << std::endl 
                                  << "TEST ERROR: CEL copy tolerance 1.0e-4 exceeded: "
                                  << tol << std::endl;
                        throw;
                    }
                }
                catch (std::exception &e) {
                    std::cout << std::endl 
                              << "TEST ERROR: Error while testing CEL copy."
                              << std::endl;
                    std::cout << e.what() << std::endl;
                    throw;
                }
                std::cout << ".";

                // Test GAL copy
                try {
                    double tol = 0.0;
                    if ((tol = wcs_copy(gal, nx, ny, crpix1, crpix2)) > 1.0e-4) {
                        std::cout << std::endl 
                                  << "TEST ERROR: GAL copy tolerance 1.0e-4 exceeded: "
                                  << tol << std::endl;
                        throw;
                    }
                }
                catch (std::exception &e) {
                    std::cout << std::endl 
                              << "TEST ERROR: Error while testing GAL copy."
                              << std::endl;
                    std::cout << e.what() << std::endl;
                    throw;
                }
                std::cout << ".";
                    
            } // endtry: main test block
            catch (std::exception &e) {
                std::cout << std::endl 
                          << "TEST ERROR: Error occured while testing "+test_class+"."
                          << std::endl;
                std::cout << e.what() << std::endl;
                throw;
            }
            std::cout << ".";
            
            // Signal final test success
            std::cout << " ok." << std::endl;

        } // endfor: looped over all non-HealPix projections

    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Error occured while testing GWcslib projections."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    
    // Exit test
    return;
}


/***************************************************************************
 *  Test: GSkymap_healpix_construct                                        *
 ***************************************************************************/
void test_GSkymap_healpix_construct(void)
{
    // Dump header
    std::cout << "Test Healpix GSkymap constructors: ";

    // Test void constructor
    try {
        GSkymap map;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to construct empty map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test correct Healpix constructors
    try {
        GSkymap ring1("HPX", "GAL", 1, "RING", 1);
        GSkymap ring2("HPX", "GAL", 2, "RING", 1);
        GSkymap ring4("HPX", "GAL", 4, "RING", 1);
        GSkymap ring8("HPX", "GAL", 8, "RING", 1);
        GSkymap ring16("HPX", "GAL", 16, "RING", 1);
        GSkymap ring32("HPX", "GAL", 32, "RING", 1);
        GSkymap ring64("HPX", "GAL", 64, "RING", 1);
        GSkymap ring128("HPX", "GAL", 128, "RING", 1);
        GSkymap ring256("HPX", "GAL", 256, "RING", 1);
        GSkymap ring512("HPX", "GAL", 512, "RING", 1);
/*
        GSkymap ring1024("HPX", "GAL", 1024, "RING", 1);
        GSkymap ring2048("HPX", "GAL", 2048, "RING", 1);
        GSkymap ring4096("HPX", "GAL", 4096, "RING", 1);
        GSkymap ring8192("HPX", "GAL", 8192, "RING", 1);
*/
        GSkymap nest1("HPX", "GAL", 1, "NEST", 1);
        GSkymap nest2("HPX", "GAL", 2, "NEST", 1);
        GSkymap nest4("HPX", "GAL", 4, "NEST", 1);
        GSkymap nest8("HPX", "GAL", 8, "NEST", 1);
        GSkymap nest16("HPX", "GAL", 16, "NEST", 1);
        GSkymap nest32("HPX", "GAL", 32, "NEST", 1);
        GSkymap nest64("HPX", "GAL", 64, "NEST", 1);
        GSkymap nest128("HPX", "GAL", 128, "NEST", 1);
        GSkymap nest256("HPX", "GAL", 256, "NEST", 1);
        GSkymap nest512("HPX", "GAL", 512, "NEST", 1);
/*
        GSkymap nest1024("HPX", "GAL", 1024, "NEST", 1);
        GSkymap nest2048("HPX", "GAL", 2048, "NEST", 1);
        GSkymap nest4096("HPX", "GAL", 4096, "NEST", 1);
        GSkymap nest8192("HPX", "GAL", 8192, "NEST", 1);
*/
        GSkymap map1("HPX", "CEL", 1, "RING", 1);
        GSkymap map2("HPX", "CEL", 1, "RING", 2);
        GSkymap map3("HPX", "EQU", 1, "RING", 1);
        GSkymap map4("HPX", "EQU", 1, "NESTED", 1);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to construct Healpix map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test Healpix copy constructor
    try {
        GSkymap ring1("HPX", "GAL", 1, "RING", 1);
        GSkymap ring2 = ring1;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to copy construct Healpix map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test invalid wcs in constructor
    try {
        GSkymap map("CAR", "GAL", 1, "RING", 1);
    }
    catch (GException::wcs_invalid &e) {
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Did not detect invalid WCS."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test invalid coordsys in constructor
    try {
        GSkymap map("HPX", "HOR", 1, "RING", 1);
    }
    catch (GException::wcs_bad_coords &e) {
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Did not detect invalid coordinate system."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test invalid nside in constructor
    try {
        GSkymap map("HPX", "GAL", 3, "RING", 1);
    }
    catch (GException::wcs_hpx_bad_nside &e) {
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Did not detect invalid nside parameter."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test invalid ordering in constructor
    try {
        GSkymap map("HPX", "GAL", 2, "SPHERICAL", 1);
    }
    catch (GException::wcs_hpx_bad_ordering &e) {
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Did not detect invalid ordering parameter."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test invalid nmaps in constructor
    try {
        GSkymap map("HPX", "GAL", 2, "NEST", 0);
    }
    catch (GException::skymap_bad_par &e) {
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Did not detect invalid ordering parameter."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 *  Test: GSkymap_healpix_io                                               *
 ***************************************************************************/
void test_GSkymap_healpix_io(void)
{
    // Set filenames
    const std::string file1 = "test_skymap_hpx_1.fits";
    const std::string file2 = "test_skymap_hpx_2.fits";

    // Dump header
    std::cout << "Test Healpix GSkymap I/O: ";

    // Define Healpix map for comparison
    GSkymap refmap("HPX", "GAL", 4, "RING", 1);

    // Test Healpix map saving
    try {
        for (int pix = 0; pix < refmap.npix(); ++pix)
            refmap(pix) = pix+1;
        refmap.save(file1, true);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to save Healpix map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test Healpix map loading
    try {
        GSkymap map;
        map.load(file1);
        int diff = 0;
        for (int pix = 0; pix < refmap.npix(); ++pix) {
            if (map(pix) != refmap(pix))
                diff++;
        }
        if (diff > 0) {
            std::cout << std::endl
                      << "TEST ERROR: Loaded file differs from saved file ."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to load Healpix map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test Healpix map instantiation
    try {
        GSkymap map(file1);
        int diff = 0;
        for (int pix = 0; pix < refmap.npix(); ++pix) {
            if (map(pix) != refmap(pix))
                diff++;
        }
        if (diff > 0) {
            std::cout << std::endl
                      << "TEST ERROR: Loaded file differs from saved file ."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to create Healpix instance"
                     " from FITS file." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 *  Test: GSkymap_wcs_construct                                            *
 ***************************************************************************/
void test_GSkymap_wcs_construct(void)
{
    // Set precision
    double eps = 1.0e-5;

    // Dump header
    std::cout << "Test WCS GSkymap constructors: ";

    // Test void constructor
    try {
        GSkymap map;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to construct empty map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test non-Healpix constructors
    try {
        GSkymap map1("CAR", "GAL", 0.0, 0.0, 1.0, 1.0, 100, 100);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct sky map(s)."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test CAR projection
    try {
        GSkymap map1("CAR", "GAL", 138.817, 37.293, 0.521, 0.931, 100, 100);
        GSkyDir dir;
        for (int l = -180; l < 180; ++l) {
            for (int b = -90; b < 90; ++b) {
                dir.lb_deg(double(l),double(b));
                GSkyPixel pixel    = map1.dir2xy(dir);
                GSkyDir   dir_back = map1.xy2dir(pixel);
                double    dist     = dir.dist_deg(dir_back);
                if (dist > eps) {
                    std::cout << std::endl
                      << "TEST ERROR: Sky direction differs:"
                      << " dir=" << dir
                      << " pixel=" << pixel
                      << " dir_back=" << dir_back
                      << " dist=" << dist << " deg" << std::endl;
                    throw;
                }
            }
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct sky map(s)."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 *  Test: GSkymap_wcs_io                                                   *
 ***************************************************************************/
void test_GSkymap_wcs_io(void)
{
    // Set filenames
    const std::string file1 = "test_skymap_wcs_1.fits";
    const std::string file2 = "test_skymap_wcs_2.fits";

    // Dump header
    std::cout << "Test WCS GSkymap I/O: ";

    // Define WCS map for comparison
    GSkymap refmap1("CAR", "GAL", 0.0, 0.0, -1.0, 1.0, 100, 100);
    GSkymap refmap2("CAR", "CEL", 0.0, 0.0, -0.3, 0.3, 200, 200);

    // Test WCS map saving
    try {
        for (int pix = 0; pix < refmap1.npix(); ++pix)
            refmap1(pix) = pix+1;
        refmap1.save(file1, 1);
        for (int pix = 0; pix < refmap2.npix(); ++pix)
            refmap2(pix) = pix+1;
        refmap2.save(file2, 1);
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to save WCS map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test WCS map loading (1)
    try {
        GSkymap map;
        map.load(file1);
        int diff = 0;
        for (int pix = 0; pix < refmap1.npix(); ++pix) {
            if (map(pix) != refmap1(pix))
                diff++;
        }
        if (diff > 0) {
            std::cout << std::endl
                      << "TEST ERROR: Loaded file "+file1+
                         " differs from saved file ."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to load WCS map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test WCS map loading (2)
    try {
        GSkymap map;
        map.load(file2);
        int diff = 0;
        for (int pix = 0; pix < refmap2.npix(); ++pix) {
            if (map(pix) != refmap2(pix))
                diff++;
        }
        if (diff > 0) {
            std::cout << std::endl
                      << "TEST ERROR: Loaded file "+file2+
                         " differs from saved file ."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to load WCS map."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "************************" << std::endl;
    std::cout << "* GSky classes testing *" << std::endl;
    std::cout << "************************" << std::endl;

    // Execute Healpix tests
    test_GWcslib();
    test_GSkymap_healpix_construct();
    test_GSkymap_healpix_io();
    test_GSkymap_wcs_construct();
    test_GSkymap_wcs_io();

    // Return
    return 0;
}
