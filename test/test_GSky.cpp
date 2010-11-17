/***************************************************************************
 *                  test_GSky.cpp  -  test GSky classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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
                      << "TEST ERROR: Sky direction differs: dir="
                      << dir << " dir_back=" << dir_back << " dist="
                      << dist << " deg" << std::endl;
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
    test_GSkymap_healpix_construct();
    test_GSkymap_healpix_io();
    test_GSkymap_wcs_construct();
    test_GSkymap_wcs_io();

    // Return
    return 0;
}
