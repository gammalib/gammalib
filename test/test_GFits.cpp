/***************************************************************************
 *                  test_GFits.cpp  -  test FITS classes                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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

/* __ Includes ___________________________________________________________ */
#include <iostream>                           // std::cout, cerr
#include <stdexcept>                          // std::exception
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "GammaLib.hpp"
#include "GTools.hpp"
#include "test_GFits.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGFits::set(void){

    // Test name
    name("GFits");

    //add tests
    add_test(static_cast<pfunction>(&TestGFits::test_create),"Test create");
    add_test(static_cast<pfunction>(&TestGFits::test_image_byte),"Test image byte");
    add_test(static_cast<pfunction>(&TestGFits::test_image_ushort),"Test image ushort");
    add_test(static_cast<pfunction>(&TestGFits::test_image_short),"Test image short");
    add_test(static_cast<pfunction>(&TestGFits::test_image_ulong),"Test image ulong");
    add_test(static_cast<pfunction>(&TestGFits::test_image_long),"Test image long");
    add_test(static_cast<pfunction>(&TestGFits::test_image_longlong),"Test image longlong");
    add_test(static_cast<pfunction>(&TestGFits::test_image_float),"Test image float");
    add_test(static_cast<pfunction>(&TestGFits::test_image_double),"Test image double");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_bit),"Test bintable bit");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_logical),"Test bintable logical");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_string),"Test bintable string");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_double),"Test bintable double");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_float),"Test bintable float");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_ushort),"Test bintable ushort");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_short),"Test bintable short");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_ulong),"Test bintable ulong");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_long),"Test bintable long");
    add_test(static_cast<pfunction>(&TestGFits::test_bintable_longlong),"Test bintable longlong");

    return;
}


/***********************************************************************//**
 * @brief Verification routines
 **************************************************************************/
int TestGFits::equal(double val, double ref, double eps)
{
    return (fabs(val-ref) < eps) ? 1 : 0;
}
int TestGFits::fequal(double val, double ref)
{
    return (fabs(val-ref) < 1.0e-5) ? 1 : 0;
}
int TestGFits::dequal(double val, double ref)
{
    return (fabs(val-ref) < 1.0e-10) ? 1 : 0;
}


/***************************************************************************
 * @brief Test FITS file creation
 *
 * @todo Add checks that verify that the file content has been saved corretly
 ***************************************************************************/
void TestGFits::test_create(void)
{
    // Remove FITS file
    int rc = 0;
    rc = system("rm -rf test_empty.fits");
    rc = system("rm -rf test_empty_image.fits");
    rc = system("rm -rf test.fits");
    rc = system("rm -rf test_create_bintable.fits");

    // Create empty FITS file
    test_try("Create empty FITS file");
    try {
        GFits fits("test_empty.fits", true);
        fits.save();
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure();
    }

    // Create FITS file with empty double precision image
    test_try("Create FITS file with empty double precision image");
    try {
        GFits fits;
        fits.open("test_empty_image.fits", true);
        GFitsImageDouble image;
        fits.append(image);
        fits.save();
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure();
    }

    // Attach double precision image
    double sum = 0.0;
    test_try("Attach double precision image");
    try {
        GFits fits;
        fits.open("test_empty_image.fits");
        int naxis       = 2;
        int nx          = 10;
        int ny          = 20;
        int naxes[]     = {nx,ny};
        GFitsImageDouble image(naxis, naxes);
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                image(ix,iy) = 0.01*ix + 0.01*iy;
                sum += image(ix,iy);
            }
        }
        fits.append(image);
        fits.saveto("test.fits");
        
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure();
    }

    // Re-open double precision image
    test_try("Re-open double precision image");
    {
        double total = 0.0;
     try{
            GFits fits;
            fits.open("test.fits");
            GFitsHDU*         hdu   = fits.hdu(0);
            GFitsImageDouble* image = (GFitsImageDouble*)fits.hdu(1);
            int nx = image->naxes(0);
            int ny = image->naxes(1);
            
            for (int ix = 0; ix < nx; ++ix) {
                for (int iy = 0; iy < ny; ++iy) {
                    total += (*image)(ix,iy);
                }
            }
            test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure();
        }
        test_assert(dequal(total, sum),"Bad values in loaded image.");
    }

    // Attach binary table (save variant)
    test_try("Attach binary table (save variant)");
    try {
        // Re-open FITS file
        GFits fits;
        fits.open("test.fits");

        // Create binary Table with 10 rows
        int nrows = 10;
        GFitsBinTable       table  = GFitsBinTable(nrows);
        GFitsTableDoubleCol first  = GFitsTableDoubleCol("First", nrows);
        GFitsTableDoubleCol second = GFitsTableDoubleCol("Second", nrows);
        first(0)  =  1.0;
        first(1)  =  2.0;
        first(2)  =  3.0;
        second(0) = -99.0;
        //std::cout << std::endl << first << std::endl;
        //std::cout << second << std::endl;
        table.append_column(first);
        table.append_column(second);

        // Append table to FILE file
        fits.append(table);

        // Save FITS file
        fits.save(true);

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure();
    }


    // Create binary table (saveto variant)
    test_try("Create binary table (saveto variant)");
    try {
        // Allocate FITS file
        GFits fits;

        // Create binary Table with 10 rows
        int nrows = 10;
        GFitsBinTable       table  = GFitsBinTable(nrows);
        GFitsTableDoubleCol first  = GFitsTableDoubleCol("First", nrows);
        GFitsTableDoubleCol second = GFitsTableDoubleCol("Second", nrows);
        first(0)  =  1.0;
        first(1)  =  2.0;
        first(2)  =  3.0;
        second(0) = -99.0;
        table.append_column(first);
        table.append_column(second);

        // Append table to FILE file
        fits.append(table);

        // Save FITS file
        fits.saveto("test_create_bintable.fits");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure();
    }

    // Return
    return;
}


/***************************************************************************
 * @brief Test GFitsImageByte class
 ***************************************************************************/
void TestGFits::test_image_byte(void)
{
    // Set filename
    std::string filename = "test_image_byte.fits";
    remove(filename.c_str());

    // Create pixel array
    unsigned char* pixels = new unsigned char[256];
    for (int i = 0; i < 256; ++i)
        pixels[i] = (unsigned char)(i);

    // Test pixel access 1D
    test_try("Test pixel access 1D");
    try {
        // Test 1D image
        GFitsImageByte image(256, pixels);
        for (int ix = 0, i = 0; ix < 256; ++ix, ++i) {

            if(!dequal(image(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix))+", expected "+str(pixels[i])+" (operator access).");
            if(!dequal(image.at(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix))+", expected "+str(pixels[i])+" (at access).");
            if(!dequal(image.pixel(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix))+", expected "+str(pixels[i])+" (pixel access).");
        }
  
        test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }

        // Test pixel access 2D
        test_try("Test pixel access 2D");
        try {
            GFitsImageByte image(16, 16, pixels);
            for (int iy = 0, i = 0; iy < 16; ++iy) {
            for (int ix = 0; ix < 16; ++ix, ++i) {
                if(!dequal(image(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy))+", expected "+str(pixels[i])+" (operator access).");

                if(!dequal(image.at(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy))+", expected "+str(pixels[i])+" (at access).");

                if(!dequal(image.pixel(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy))+", expected "+str(pixels[i])+" (pixel access).");
            }
            }
        test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }

        // Test pixel access 3D
        test_try("Test pixel access 3D");
        try {
            GFitsImageByte image(16, 4, 4, pixels);
            for (int iz = 0, i = 0; iz < 4; ++iz) {
            for (int iy = 0; iy < 4; ++iy) {
            for (int ix = 0; ix < 16; ++ix, ++i) {

                if (!dequal(image(ix,iy,iz), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz))+", expected "+str(pixels[i])+" (operator access).");

                if (!dequal(image.at(ix,iy,iz), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz))+", expected "+str(pixels[i])+" (at access).");

                if (!dequal(image.pixel(ix,iy,iz), pixels[i])) {
                    throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz))+", expected "+str(pixels[i])+" (pixel access).");
                }
            }
            }
            }
        test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }

        // Test 4D image
        test_try("Test pixel access 4D");
        try {
            GFitsImageByte image (4, 4, 4, 4, pixels);
            for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
            for (int iy = 0; iy < 4; ++iy) {
            for (int ix = 0; ix < 4; ++ix, ++i) {
                if (!dequal(image(ix,iy,iz,it), pixels[i])) {
                    throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz,it))+", expected "+str(pixels[i])+" (operator access).");
                }
                if (!dequal(image.at(ix,iy,iz,it), pixels[i])) {
                    throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz,it))+", expected "+str(pixels[i])+" (at access).");
                }
                if (!dequal(image.pixel(ix,iy,iz,it), pixels[i])) {
                    throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                }
            }
            }
            }
            }
            test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }



    // Test image I/O with 4D image
    test_try("Test image I/O with 4D image");
    try{
        // Create 4D image
        int naxes[] = {4,4,4,4};
        GFitsImageByte image(4, naxes, pixels);

        // Save image
        GFits fits(filename, true);
        fits.append(image);
        fits.save();

        // Open FITS image
        GFits infile(filename);
        GFitsImage* ptr = infile.image(0);

        // Test 4D image
        for (int it = 0, i = 0; it < 4; ++it) {
        for (int iz = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 4; ++ix, ++i) {
            if (!dequal(ptr->pixel(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(ptr->pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        }

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***************************************************************************
 * @brief Test GFitsImageUShort class
 ***************************************************************************/
void TestGFits::test_image_ushort(void)
{

    // Set filename
    std::string filename = "test_image_ushort.fits";
    remove(filename.c_str());

    // Create pixel array
    unsigned short* pixels = new unsigned short[256];
    for (int i = 0; i < 256; ++i)
        pixels[i] = (unsigned short)(i);

    // Test pixel access (1D to 4D)
    test_try("Test pixel access 1D");
    try {

        // Test 1D image
        GFitsImageUShort image(256, pixels);
        for (int ix = 0, i = 0; ix < 256; ++ix, ++i) {
            if(!dequal(image(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix))+", expected "+str(pixels[i])+" (operator access).");
            if(!dequal(image.at(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix))+", expected "+str(pixels[i])+" (at access).");
            if(!dequal(image.pixel(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix))+", expected "+str(pixels[i])+" (pixel access).");
        }

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
     // Test 2D image
    test_try("Test pixel access 2D");
    try {
        GFitsImageUShort image(16, 16, pixels);
        for (int iy = 0, i = 0; iy < 16; ++iy) {
        for (int ix = 0; ix < 16; ++ix, ++i) {
            if(!dequal(image(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy))+", expected "+str(pixels[i])+" (operator access).");

            if(!dequal(image.at(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy))+", expected "+str(pixels[i])+" (at access).");

            if(!dequal(image.pixel(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy))+", expected "+str(pixels[i])+" (pixel access).");
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 3D image
    test_try("Test pixel access 3D");
    try {
        GFitsImageUShort image(16, 4, 4, pixels);
        for (int iz = 0, i = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 16; ++ix, ++i) {
            if (!dequal(image(ix,iy,iz), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz))+", expected "+str(pixels[i])+" (operator access).");

            if (!dequal(image.at(ix,iy,iz), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz))+", expected "+str(pixels[i])+" (at access).");

            if (!dequal(image.pixel(ix,iy,iz), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 4D image
    test_try("Test pixel access 4D");
    try {
        GFitsImageUShort image(4, 4, 4, 4, pixels);
        for (int it = 0, i = 0; it < 4; ++it) {
        for (int iz = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 4; ++ix, ++i) {
            if (!dequal(image(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz,it))+", expected "+str(pixels[i])+" (operator access).");
            }
            if (!dequal(image.at(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz,it))+", expected "+str(pixels[i])+" (at access).");
            }
            if (!dequal(image.pixel(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test image I/O with 4D image
    test_try("Test image I/O with 4D image");
    try{
        // Create 4D image
        int naxes[] = {4,4,4,4};
        GFitsImageUShort image(4, naxes, pixels);

        // Save image
        GFits fits(filename, true);
        fits.append(image);
        fits.save();

        // Open FITS image
        GFits infile(filename);
        GFitsImage* ptr = infile.image(0);

        // Test 4D image
        for (int it = 0, i = 0; it < 4; ++it) {
        for (int iz = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 4; ++ix, ++i) {
            if (!dequal(ptr->pixel(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(ptr->pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Return
    return;
}


/***************************************************************************
 * @brief Test GFitsImageShort class
 ***************************************************************************/
void TestGFits::test_image_short(void)
{
    // Set filename
    std::string filename = "test_image_short.fits";
    remove(filename.c_str());

    // Create pixel array
    short* pixels = new short[256];
    for (int i = 0; i < 256; ++i)
        pixels[i] = (short)(i);

    // Test pixel access 1D
    test_try("Test pixel access 1D");
    try {

        // Test 1D image
        GFitsImageShort image(256, pixels);
        for (int ix = 0, i = 0; ix < 256; ++ix, ++i) {
            if(!dequal(image(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix))+", expected "+str(pixels[i])+" (operator access).");
            if(!dequal(image.at(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix))+", expected "+str(pixels[i])+" (at access).");
            if(!dequal(image.pixel(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix))+", expected "+str(pixels[i])+" (pixel access).");
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test 2D image
    test_try("Test pixel access 2D");
    try {
        GFitsImageShort image(16, 16, pixels);
        for (int iy = 0, i = 0; iy < 16; ++iy) {
        for (int ix = 0; ix < 16; ++ix, ++i) {
            if(!dequal(image(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy))+", expected "+str(pixels[i])+" (operator access).");

            if(!dequal(image.at(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy))+", expected "+str(pixels[i])+" (at access).");

            if(!dequal(image.pixel(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy))+", expected "+str(pixels[i])+" (pixel access).");
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 3D image
    test_try("Test pixel access 3D");
    try {
        GFitsImageShort image(16, 4, 4, pixels);
        for (int iz = 0, i = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 16; ++ix, ++i) {
            if (!dequal(image(ix,iy,iz), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz))+", expected "+str(pixels[i])+" (operator access).");

            if (!dequal(image.at(ix,iy,iz), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz))+", expected "+str(pixels[i])+" (at access).");

            if (!dequal(image.pixel(ix,iy,iz), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 4D image
    test_try("Test pixel access 4D");
    try {
        GFitsImageShort image(4, 4, 4, 4, pixels);
        for (int it = 0, i = 0; it < 4; ++it) {
        for (int iz = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 4; ++ix, ++i) {
            if (!dequal(image(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz,it))+", expected "+str(pixels[i])+" (operator access).");
            }
            if (!dequal(image.at(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz,it))+", expected "+str(pixels[i])+" (at access).");
            }
            if (!dequal(image.pixel(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test image I/O with 4D image
    test_try("Test image I/O with 4D image");
    try{
        // Create 4D image
        int naxes[] = {4,4,4,4};
        GFitsImageShort image(4, naxes, pixels);

        // Save image
        GFits fits(filename, true);
        fits.append(image);
        fits.save();

        // Open FITS image
        GFits infile(filename);
        GFitsImage* ptr = infile.image(0);

        // Test 4D image
        for (int it = 0, i = 0; it < 4; ++it) {
        for (int iz = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 4; ++ix, ++i) {
            if (!dequal(ptr->pixel(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(ptr->pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***************************************************************************
 * @brief Test GFitsImageULong class
 ***************************************************************************/
void TestGFits::test_image_ulong(void)
{
    // Set filename
    std::string filename = "test_image_ulong.fits";
    remove(filename.c_str());

    // Create pixel array
    unsigned long* pixels = new unsigned long[256];
    for (int i = 0; i < 256; ++i)
        pixels[i] = (unsigned long)(i);

    // Test pixel access (1D to 4D)
    test_try("Test pixel access 1D");
    try {

        // Test 1D image
        GFitsImageULong image(256, pixels);
        for (int ix = 0, i = 0; ix < 256; ++ix, ++i) {
            if(!dequal(image(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix))+", expected "+str(pixels[i])+" (operator access).");
            if(!dequal(image.at(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix))+", expected "+str(pixels[i])+" (at access).");
            if(!dequal(image.pixel(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix))+", expected "+str(pixels[i])+" (pixel access).");
        }

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 2D image
    test_try("Test pixel access 2D");
    try {
        GFitsImageULong image(16, 16, pixels);
        for (int iy = 0, i = 0; iy < 16; ++iy) {
        for (int ix = 0; ix < 16; ++ix, ++i) {
            if(!dequal(image(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy))+", expected "+str(pixels[i])+" (operator access).");

            if(!dequal(image.at(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy))+", expected "+str(pixels[i])+" (at access).");

            if(!dequal(image.pixel(ix,iy), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy))+", expected "+str(pixels[i])+" (pixel access).");
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

           // Test 3D image
    test_try("Test pixel access 3D");
    try{
        GFitsImageULong image(16, 4, 4, pixels);
        for (int iz = 0, i = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 16; ++ix, ++i) {
            if (!dequal(image(ix,iy,iz), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz))+", expected "+str(pixels[i])+" (operator access).");

            if (!dequal(image.at(ix,iy,iz), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz))+", expected "+str(pixels[i])+" (at access).");

            if (!dequal(image.pixel(ix,iy,iz), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 4D image
    test_try("Test pixel access 4D");
    try {
        GFitsImageULong image(4, 4, 4, 4, pixels);
        for (int it = 0, i = 0; it < 4; ++it) {
        for (int iz = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 4; ++ix, ++i) {
            if (!dequal(image(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz,it))+", expected "+str(pixels[i])+" (operator access).");
            }
            if (!dequal(image.at(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz,it))+", expected "+str(pixels[i])+" (at access).");
            }
            if (!dequal(image.pixel(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }


    // Test image I/O with 4D image
    test_try("Test image I/O with 4D image");
    try{
        // Create 4D image
        int naxes[] = {4,4,4,4};
        GFitsImageULong image(4, naxes, pixels);

        // Save image
        GFits fits(filename, true);
        fits.append(image);
        fits.save();

        // Open FITS image
        GFits infile(filename);
        GFitsImage* ptr = infile.image(0);
        std::cout << ".";

        // Test 4D image
        for (int it = 0, i = 0; it < 4; ++it) {
        for (int iz = 0; iz < 4; ++iz) {
        for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 4; ++ix, ++i) {
            if (!dequal(ptr->pixel(ix,iy,iz,it), pixels[i])) {
                throw exception_failure("Unexpected pixel content (has "+str(ptr->pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        }
        }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***************************************************************************
 * @brief Test GFitsImageLong class
 ***************************************************************************/
void TestGFits::test_image_long(void)
{
    // Set filename
    std::string filename = "test_image_long.fits";
    remove(filename.c_str());

    // Create pixel array
    long* pixels = new long[256];
    for (int i = 0; i < 256; ++i)
        pixels[i] = (long)(i);

    // Test pixel access (1D to 4D)
    test_try("Test pixel access 1D");
    try {

        // Test 1D image
        GFitsImageLong image(256, pixels);
        for (int ix = 0, i = 0; ix < 256; ++ix, ++i) {
            if(!dequal(image(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix))+", expected "+str(pixels[i])+" (operator access).");
            if(!dequal(image.at(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix))+", expected "+str(pixels[i])+" (at access).");
            if(!dequal(image.pixel(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix))+", expected "+str(pixels[i])+" (pixel access).");
        }

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 2D image
    test_try("Test pixel access 2D");
    try {
        GFitsImageLong image(16, 16, pixels);
        for (int iy = 0, i = 0; iy < 16; ++iy) {
            for (int ix = 0; ix < 16; ++ix, ++i) {
                if(!dequal(image(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy))+", expected "+str(pixels[i])+" (operator access).");

                if(!dequal(image.at(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy))+", expected "+str(pixels[i])+" (at access).");

                if(!dequal(image.pixel(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

           // Test 3D image
    test_try("Test pixel access 3D");
    try{
        GFitsImageLong image(16, 4, 4, pixels);
        for (int iz = 0, i = 0; iz < 4; ++iz) {
            for (int iy = 0; iy < 4; ++iy) {
                for (int ix = 0; ix < 16; ++ix, ++i) {
                    if (!dequal(image(ix,iy,iz), pixels[i]))
                        throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz))+", expected "+str(pixels[i])+" (operator access).");
    
                    if (!dequal(image.at(ix,iy,iz), pixels[i]))
                        throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz))+", expected "+str(pixels[i])+" (at access).");
    
                    if (!dequal(image.pixel(ix,iy,iz), pixels[i])) {
                        throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz))+", expected "+str(pixels[i])+" (pixel access).");
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 4D image
    test_try("Test pixel access 4D");
    try {
        GFitsImageLong image(4, 4, 4, 4, pixels);
        for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
                for (int iy = 0; iy < 4; ++iy) {
                    for (int ix = 0; ix < 4; ++ix, ++i) {
                        if (!dequal(image(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz,it))+", expected "+str(pixels[i])+" (operator access).");
                        }
                        if (!dequal(image.at(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz,it))+", expected "+str(pixels[i])+" (at access).");
                        }
                        if (!dequal(image.pixel(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                        }
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }


    // Test image I/O with 4D image
    test_try("Test image I/O with 4D image");
    try{
        // Create 4D image
        int naxes[] = {4,4,4,4};
        GFitsImageLong image(4, naxes, pixels);

        // Save image
        GFits fits(filename, true);
        fits.append(image);
        fits.save();

        // Open FITS image
        GFits infile(filename);
        GFitsImage* ptr = infile.image(0);
        std::cout << ".";

        // Test 4D image
        for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
                for (int iy = 0; iy < 4; ++iy) {
                    for (int ix = 0; ix < 4; ++ix, ++i) {
                        if (!dequal(ptr->pixel(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(ptr->pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                        }
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***************************************************************************
 * @brief Test GFitsImageLongLong class
 ***************************************************************************/

void TestGFits::test_image_longlong(void)
{
    // Set filename
    std::string filename = "test_image_longlong.fits";
    remove(filename.c_str());

    // Create pixel array
    long long* pixels = new long long[256];
    for (int i = 0; i < 256; ++i)
        pixels[i] = (long long)(i);

    // Test pixel access (1D to 4D)
    test_try("Test pixel access 1D");
    try {

        // Test 1D image
        GFitsImageLongLong image(256, pixels);
        for (int ix = 0, i = 0; ix < 256; ++ix, ++i) {
            if(!dequal(image(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix))+", expected "+str(pixels[i])+" (operator access).");
            if(!dequal(image.at(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix))+", expected "+str(pixels[i])+" (at access).");
            if(!dequal(image.pixel(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix))+", expected "+str(pixels[i])+" (pixel access).");
        }

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 2D image
    test_try("Test pixel access 2D");
    try {
        GFitsImageLongLong image(16, 16, pixels);
        for (int iy = 0, i = 0; iy < 16; ++iy) {
            for (int ix = 0; ix < 16; ++ix, ++i) {
                if(!dequal(image(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy))+", expected "+str(pixels[i])+" (operator access).");

                if(!dequal(image.at(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy))+", expected "+str(pixels[i])+" (at access).");

                if(!dequal(image.pixel(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

           // Test 3D image
    test_try("Test pixel access 3D");
    try{
        GFitsImageLongLong image(16, 4, 4, pixels);
        for (int iz = 0, i = 0; iz < 4; ++iz) {
            for (int iy = 0; iy < 4; ++iy) {
                for (int ix = 0; ix < 16; ++ix, ++i) {
                    if (!dequal(image(ix,iy,iz), pixels[i]))
                        throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz))+", expected "+str(pixels[i])+" (operator access).");
    
                    if (!dequal(image.at(ix,iy,iz), pixels[i]))
                        throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz))+", expected "+str(pixels[i])+" (at access).");
    
                    if (!dequal(image.pixel(ix,iy,iz), pixels[i])) {
                        throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz))+", expected "+str(pixels[i])+" (pixel access).");
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 4D image
    test_try("Test pixel access 4D");
    try {
        GFitsImageLongLong image(4, 4, 4, 4, pixels);
        for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
                for (int iy = 0; iy < 4; ++iy) {
                    for (int ix = 0; ix < 4; ++ix, ++i) {
                        if (!dequal(image(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz,it))+", expected "+str(pixels[i])+" (operator access).");
                        }
                        if (!dequal(image.at(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz,it))+", expected "+str(pixels[i])+" (at access).");
                        }
                        if (!dequal(image.pixel(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                        }
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }


    // Test image I/O with 4D image
    test_try("Test image I/O with 4D image");
    try{
        // Create 4D image
        int naxes[] = {4,4,4,4};
        GFitsImageLongLong image(4, naxes, pixels);

        // Save image
        GFits fits(filename, true);
        fits.append(image);
        fits.save();

        // Open FITS image
        GFits infile(filename);
        GFitsImage* ptr = infile.image(0);
        std::cout << ".";

        // Test 4D image
        for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
                for (int iy = 0; iy < 4; ++iy) {
                    for (int ix = 0; ix < 4; ++ix, ++i) {
                        if (!dequal(ptr->pixel(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(ptr->pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                        }
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}

/***************************************************************************
 * @brief Test GFitsImageFloat class
 ***************************************************************************/
void TestGFits::test_image_float(void)
{
    // Set filename
    std::string filename = "test_image_float.fits";
    remove(filename.c_str());

    // Create pixel array
    float* pixels = new float[256];
    for (int i = 0; i < 256; ++i)
        pixels[i] = (float)(i);

    // Test pixel access (1D to 4D)
    test_try("Test pixel access 1D");
    try {

        // Test 1D image
        GFitsImageFloat image(256, pixels);
        for (int ix = 0, i = 0; ix < 256; ++ix, ++i) {
            if(!dequal(image(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix))+", expected "+str(pixels[i])+" (operator access).");
            if(!dequal(image.at(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix))+", expected "+str(pixels[i])+" (at access).");
            if(!dequal(image.pixel(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix))+", expected "+str(pixels[i])+" (pixel access).");
        }

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 2D image
    test_try("Test pixel access 2D");
    try {
        GFitsImageFloat image(16, 16, pixels);
        for (int iy = 0, i = 0; iy < 16; ++iy) {
            for (int ix = 0; ix < 16; ++ix, ++i) {
                if(!dequal(image(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy))+", expected "+str(pixels[i])+" (operator access).");

                if(!dequal(image.at(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy))+", expected "+str(pixels[i])+" (at access).");

                if(!dequal(image.pixel(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

           // Test 3D image
    test_try("Test pixel access 3D");
    try{
        GFitsImageFloat image(16, 4, 4, pixels);
        for (int iz = 0, i = 0; iz < 4; ++iz) {
            for (int iy = 0; iy < 4; ++iy) {
                for (int ix = 0; ix < 16; ++ix, ++i) {
                    if (!dequal(image(ix,iy,iz), pixels[i]))
                        throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz))+", expected "+str(pixels[i])+" (operator access).");
    
                    if (!dequal(image.at(ix,iy,iz), pixels[i]))
                        throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz))+", expected "+str(pixels[i])+" (at access).");
    
                    if (!dequal(image.pixel(ix,iy,iz), pixels[i])) {
                        throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz))+", expected "+str(pixels[i])+" (pixel access).");
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 4D image
    test_try("Test pixel access 4D");
    try {
        GFitsImageFloat image(4, 4, 4, 4, pixels);
        for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
                for (int iy = 0; iy < 4; ++iy) {
                    for (int ix = 0; ix < 4; ++ix, ++i) {
                        if (!dequal(image(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz,it))+", expected "+str(pixels[i])+" (operator access).");
                        }
                        if (!dequal(image.at(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz,it))+", expected "+str(pixels[i])+" (at access).");
                        }
                        if (!dequal(image.pixel(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                        }
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }


    // Test image I/O with 4D image
    test_try("Test image I/O with 4D image");
    try{
        // Create 4D image
        int naxes[] = {4,4,4,4};
        GFitsImageFloat image(4, naxes, pixels);

        // Save image
        GFits fits(filename, true);
        fits.append(image);
        fits.save();

        // Open FITS image
        GFits infile(filename);
        GFitsImage* ptr = infile.image(0);
        std::cout << ".";

        // Test 4D image
        for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
                for (int iy = 0; iy < 4; ++iy) {
                    for (int ix = 0; ix < 4; ++ix, ++i) {
                        if (!dequal(ptr->pixel(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(ptr->pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                        }
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***************************************************************************
 * @brief Test GFitsImageDouble class
 ***************************************************************************/

void TestGFits::test_image_double(void)
{
    // Set filename
    std::string filename = "test_image_double.fits";
    remove(filename.c_str());

    // Create pixel array
    double* pixels = new double[256];
    for (int i = 0; i < 256; ++i)
        pixels[i] = (double)(i);

    // Test pixel access (1D to 4D)
    test_try("Test pixel access 1D");
    try {

        // Test 1D image
        GFitsImageDouble image(256, pixels);
        for (int ix = 0, i = 0; ix < 256; ++ix, ++i) {
            if(!dequal(image(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image(ix))+", expected "+str(pixels[i])+" (operator access).");
            if(!dequal(image.at(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.at(ix))+", expected "+str(pixels[i])+" (at access).");
            if(!dequal(image.pixel(ix), pixels[i]))
                throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix))+", expected "+str(pixels[i])+" (pixel access).");
        }

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 2D image
    test_try("Test pixel access 2D");
    try {
        GFitsImageDouble image(16, 16, pixels);
        for (int iy = 0, i = 0; iy < 16; ++iy) {
            for (int ix = 0; ix < 16; ++ix, ++i) {
                if(!dequal(image(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy))+", expected "+str(pixels[i])+" (operator access).");

                if(!dequal(image.at(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy))+", expected "+str(pixels[i])+" (at access).");

                if(!dequal(image.pixel(ix,iy), pixels[i]))
                    throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy))+", expected "+str(pixels[i])+" (pixel access).");
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

           // Test 3D image
    test_try("Test pixel access 3D");
    try{
        GFitsImageDouble image(16, 4, 4, pixels);
        for (int iz = 0, i = 0; iz < 4; ++iz) {
            for (int iy = 0; iy < 4; ++iy) {
                for (int ix = 0; ix < 16; ++ix, ++i) {
                    if (!dequal(image(ix,iy,iz), pixels[i]))
                        throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz))+", expected "+str(pixels[i])+" (operator access).");
    
                    if (!dequal(image.at(ix,iy,iz), pixels[i]))
                        throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz))+", expected "+str(pixels[i])+" (at access).");
    
                    if (!dequal(image.pixel(ix,iy,iz), pixels[i])) {
                        throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz))+", expected "+str(pixels[i])+" (pixel access).");
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test 4D image
    test_try("Test pixel access 4D");
    try {
        GFitsImageDouble image(4, 4, 4, 4, pixels);
        for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
                for (int iy = 0; iy < 4; ++iy) {
                    for (int ix = 0; ix < 4; ++ix, ++i) {
                        if (!dequal(image(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image(ix,iy,iz,it))+", expected "+str(pixels[i])+" (operator access).");
                        }
                        if (!dequal(image.at(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image.at(ix,iy,iz,it))+", expected "+str(pixels[i])+" (at access).");
                        }
                        if (!dequal(image.pixel(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(image.pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                        }
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }


    // Test image I/O with 4D image
    test_try("Test image I/O with 4D image");
    try{
        // Create 4D image
        int naxes[] = {4,4,4,4};
        GFitsImageDouble image(4, naxes, pixels);

        // Save image
        GFits fits(filename, true);
        fits.append(image);
        fits.save();

        // Open FITS image
        GFits infile(filename);
        GFitsImage* ptr = infile.image(0);
        std::cout << ".";

        // Test 4D image
        for (int it = 0, i = 0; it < 4; ++it) {
            for (int iz = 0; iz < 4; ++iz) {
                for (int iy = 0; iy < 4; ++iy) {
                    for (int ix = 0; ix < 4; ++ix, ++i) {
                        if (!dequal(ptr->pixel(ix,iy,iz,it), pixels[i])) {
                            throw exception_failure("Unexpected pixel content (has "+str(ptr->pixel(ix,iy,iz,it))+", expected "+str(pixels[i])+" (pixel access).");
                        }
                    }
                }
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}

/***************************************************************************
 * @brief Test double precision FITS binary table
 ***************************************************************************/

void TestGFits::test_bintable_double(void)
{

    // Set filename
    std::string filename = "test_bintable_double.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    double      sum_dbl       = 0.0;
    double      sum_dbl10     = 0.0;
    int         sum_dbl_int   = 0;
    int         sum_dbl10_int = 0;
    std::string sum_dbl_str;
    std::string sum_dbl10_str;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;
    
        //
        // ===== D O U B L E =====
        //
    
        GFitsTableDoubleCol col_dbl;
    
        test_try("Double");
        try{
             // Set string table
            test_try("Set table");
            try{
                col_dbl = GFitsTableDoubleCol("DOUBLE", nrows);
                for (int i = 0; i < nrows; ++i) {
                    double val_dbl = cos(double(i));
                    col_dbl(i)   = val_dbl;
                    sum_dbl     += val_dbl;
                    sum_dbl_int += int(val_dbl);
                    sum_dbl_str += ":"+str(val_dbl);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_dbl(i);
                if (!dequal(tot_dbl, sum_dbl)) {
                    throw exception_failure("GFitsTableDoubleCol::operator() : Reference sum: "+str(sum_dbl)+"  Derived sum:   "+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_dbl.real(i);
                if (!dequal(tot_dbl, sum_dbl)) {
                    throw exception_failure("GFitsTableDoubleCol::real() : Reference sum: "+str(sum_dbl)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_dbl.integer(i);
                if (tot_int != sum_dbl_int) {
                    throw exception_failure("GFitsTableDoubleCol::integer() : Reference sum: "+str(sum_dbl_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_dbl.string(i);
                if (tot_str != sum_dbl_str) {
                    throw exception_failure("GFitsTableDoubleCol::string() : Reference sum: "+sum_dbl_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try ulong
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== D O U B L E 1 0 =====
        //
    
        GFitsTableDoubleCol col_dbl10;
            
        test_try("Double 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col_dbl10 = GFitsTableDoubleCol("DOUBLE10", nrows,nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        double val_dbl = cos(double(i))*cos(0.33*double(j));
                        col_dbl10(i,j) = val_dbl;
                        sum_dbl10     += val_dbl;
                        sum_dbl10_int += int(val_dbl);
                        sum_dbl10_str += ":"+str(val_dbl);
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_dbl10(i,j);
                }
                if (!dequal(tot_dbl, sum_dbl10)) {
                    throw exception_failure("GFitsTableDoubleCol::operator() - 10 : Reference sum: "+str(sum_dbl10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_dbl10.real(i,j);
                }
                if (!dequal(tot_dbl, sum_dbl10)) {
                    throw exception_failure("GFitsTableDoubleCol::real() - 10 : Reference sum: "+str(sum_dbl10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_dbl10.integer(i,j);
                }
                if (tot_int != sum_dbl10_int) {
                    throw exception_failure("GFitsTableDoubleCol::integer() - 10 : Reference sum: "+str(sum_dbl10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_dbl10.string(i,j);
                }
                if (tot_str != sum_dbl10_str) {
                    throw exception_failure("GFitsTableDoubleCol::string() - 10 : Reference sum: "+sum_dbl10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try double 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col_dbl);
            table.insert_column(1, col_dbl10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;
        
        //
        // ===== D O U B L E =====
        //
        GFitsTableDoubleCol* col_dbl = NULL;
        test_try("Double");
        try {
            // Get column
            col_dbl = (GFitsTableDoubleCol*)&(*fits.table(1))["DOUBLE"];

            // Check table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += (*col_dbl)(i);
                if (!dequal(tot_dbl, sum_dbl)) {
                    throw exception_failure("GFitsTableDoubleCol::operator() : Reference sum: "+str(sum_dbl)+"  Derived sum:   "+str(tot_dbl));
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += (*col_dbl).real(i);
                if (!dequal(tot_dbl, sum_dbl)) {
                    throw exception_failure("GFitsTableDoubleCol::real() : Reference sum: "+str(sum_dbl)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += (*col_dbl).integer(i);
                if (tot_int != sum_dbl_int) {
                    throw exception_failure("GFitsTableDoubleCol::integer() : Reference sum: "+str(sum_dbl_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+ (*col_dbl).string(i);
                if (tot_str != sum_dbl_str) {
                    throw exception_failure("GFitsTableDoubleCol::string() : Reference sum: "+sum_dbl_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== D O U B L E 1 0 =====
        //
        GFitsTableDoubleCol* col_dbl10 = NULL;
        test_try("Double 10");
        try{
            // Get column
            col_dbl10 = (GFitsTableDoubleCol*)&(*fits.table(1))["DOUBLE10"];
            
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += (*col_dbl10)(i,j);
                }
                if (!dequal(tot_dbl, sum_dbl10)) {
                    throw exception_failure("GFitsTableDoubleCol::operator() - 10 : Reference sum: "+str(sum_dbl10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_dbl10->real(i,j);
                }
                if (!dequal(tot_dbl, sum_dbl10)) {
                    throw exception_failure("GFitsTableDoubleCol::real() - 10 : Reference sum: "+str(sum_dbl10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_dbl10->integer(i,j);
                }
                if (tot_int != sum_dbl10_int) {
                    throw exception_failure("GFitsTableDoubleCol::integer() - 10 : Reference sum: "+str(sum_dbl10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_dbl10->string(i,j);
                }
                if (tot_str != sum_dbl10_str) {
                    throw exception_failure("GFitsTableDoubleCol::string() - 10 : Reference sum: "+sum_dbl10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try float 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }


}

/***************************************************************************
 * @brief Test single precision FITS binary table
 ***************************************************************************/

void TestGFits::test_bintable_float(void)
{

    // Set filename
    std::string filename = "test_bintable_float.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    float       sum_flt       = 0.0;
    float       sum_flt10     = 0.0;
    int         sum_flt_int   = 0;
    int         sum_flt10_int = 0;
    std::string sum_flt_str;
    std::string sum_flt10_str;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;
    
        //
        //
        // ===== F L O A T =====
        //
    
        GFitsTableFloatCol col_flt;
    
        test_try("Float");
        try{
             // Set string table
            test_try("Set table");
            try{
                col_flt = GFitsTableFloatCol("FLOAT", nrows);
                for (int i = 0; i < nrows; ++i) {
                    float val_flt = cos(0.1*float(i));
                    col_flt(i)   = val_flt;
                    sum_flt     += val_flt;
                    sum_flt_int += int(val_flt);
                    sum_flt_str += ":"+str(val_flt);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_flt = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_flt += col_flt(i);
                if (!fequal(tot_flt, sum_flt)) {
                    throw exception_failure("GFitsTableFloatCol::operator() : Reference sum: "+str(sum_flt)+"  Derived sum:   "+str(tot_flt));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_flt = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_flt += col_flt.real(i);
                if (!fequal(tot_flt, sum_flt)) {
                    throw exception_failure("GFitsTableFloatCol::real() : Reference sum: "+str(sum_flt)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_flt.integer(i);
                if (tot_int != sum_flt_int) {
                    throw exception_failure("GFitsTableFloatCol::integer() : Reference sum: "+str(sum_flt_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_flt.string(i);
                if (tot_str != sum_flt_str) {
                    throw exception_failure("GFitsTableFloatCol::string() : Reference sum: "+sum_flt_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try ulong
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== F L O A T  1 0 =====
        //
    
        GFitsTableFloatCol col_flt10;
            
        test_try("Float 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col_flt10 = GFitsTableFloatCol("FLOAT10", nrows,nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        float val_flt  = cos(0.1*float(i))*cos(0.33*float(j));
                        col_flt10(i,j) = val_flt;
                        sum_flt10     += val_flt;
                        sum_flt10_int += int(val_flt);
                        sum_flt10_str += ":"+str(val_flt);
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_flt = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_flt += col_flt10(i,j);
                }
                if (!fequal(tot_flt, sum_flt10)) {
                    throw exception_failure("GFitsTableFloatCol::operator() - 10 : Reference sum: "+str(sum_flt10)+"Derived sum:"+str(tot_flt));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_flt = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_flt += col_flt10.real(i,j);
                }
                if (!fequal(tot_flt, sum_flt10)) {
                    throw exception_failure("GFitsTableFloatCol::real() - 10 : Reference sum: "+str(sum_flt10)+"Derived sum:"+str(tot_flt));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_flt10.integer(i,j);
                }
                if (tot_int != sum_flt10_int) {
                    throw exception_failure("GFitsTableFloatCol::integer() - 10 : Reference sum: "+str(sum_flt10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_flt10.string(i,j);
                }
                if (tot_str != sum_flt10_str) {
                    throw exception_failure("GFitsTableFloatCol::string() - 10 : Reference sum: "+sum_flt10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try short 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col_flt);
            table.insert_column(0, col_flt10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;
        
        //
        // ===== U S H O R T =====
        //
        GFitsTableFloatCol* col_flt = NULL;
        test_try("Float");
        try {
            // Get column
            col_flt = (GFitsTableFloatCol*)&(*fits.table(1))["FLOAT"];

            // Check table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_flt = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_flt += (*col_flt)(i);
                if (!fequal(tot_flt, sum_flt)) {
                    throw exception_failure("GFitsTableFloatCol::operator() : Reference sum: "+str(sum_flt)+"  Derived sum:   "+str(tot_flt));
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_flt = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_flt += col_flt->real(i);
                if (!fequal(tot_flt, sum_flt)) {
                    throw exception_failure("GFitsTableFloatCol::real() : Reference sum: "+str(sum_flt)+"Derived sum:"+str(tot_flt));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_flt->integer(i);
                if (tot_int != sum_flt_int) {
                    throw exception_failure("GFitsTableFloatCol::integer() : Reference sum: "+str(sum_flt_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_flt->string(i);
                if (tot_str != sum_flt_str) {
                    throw exception_failure("GFitsTableFloatCol::string() : Reference sum: "+sum_flt_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== F L O A T  1 0 =====
        //
        GFitsTableFloatCol* col_flt10 = NULL;
        test_try("Float 10");
        try{
            // Get column
            col_flt10 = (GFitsTableFloatCol*)&(*fits.table(1))["FLOAT10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_flt = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_flt += (*col_flt10)(i,j);
                }
                if (!fequal(tot_flt, sum_flt10)) {
                    throw exception_failure("GFitsTableFloatCol::operator() - 10 : Reference sum: "+str(sum_flt10)+"Derived sum:"+str(tot_flt));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_flt = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_flt += col_flt10->real(i,j);
                }
                if (!fequal(tot_flt, sum_flt10)) {
                    throw exception_failure("GFitsTableFloatCol::real() - 10 : Reference sum: "+str(sum_flt10)+"Derived sum:"+str(tot_flt));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_flt10->integer(i,j);
                }
                if (tot_int != sum_flt10_int) {
                    throw exception_failure("GFitsTableFloatCol::integer() - 10 : Reference sum: "+str(sum_flt10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_flt10->string(i,j);
                }
                if (tot_str != sum_flt10_str) {
                    throw exception_failure("GFitsTableFloatCol::string() - 10 : Reference sum: "+sum_flt10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try float 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }


}


/***************************************************************************
 * @brief Test short FITS binary table
 ***************************************************************************/
 
void TestGFits::test_bintable_short(void)
{

    // Set filename
    std::string filename = "test_bintable_short.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    short       sum_sht       = 0;
    short       sum_sht10     = 0;
    short       sum_sht_int   = 0;
    short       sum_sht10_int = 0;
    std::string sum_sht_str;
    std::string sum_sht10_str;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;
    
        //
        //
        // ===== S H O R T =====
        //
    
        GFitsTableShortCol col_sht;
    
        test_try("Short");
        try{
             // Set string table
            test_try("Set table");
            try{
                col_sht = GFitsTableShortCol("SHORT", nrows);
                for (int i = 0; i < nrows; ++i) {
                    short val_sht = short(1000.0 * cos(0.1*float(i)));
                    col_sht(i)    = val_sht;
                    sum_sht      += val_sht;
                    sum_sht_int  += int(val_sht);
                    sum_sht_str += ":"+str(val_sht);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_sht = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_sht += col_sht(i);
                if (tot_sht != sum_sht) {
                    throw exception_failure("GFitsTableShortCol::operator() : Reference sum: "+str(sum_sht)+"  Derived sum:   "+str(tot_sht));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_sht.real(i);
                if (!dequal(tot_dbl, double(sum_sht))) {
                    throw exception_failure("GFitsTableShortCol::real() : Reference sum: "+str(sum_sht)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_sht.integer(i);
                if (tot_int != sum_sht_int) {
                    throw exception_failure("GFitsTableShortCol::integer() : Reference sum: "+str(sum_sht_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_sht.string(i);
                if (tot_str != sum_sht_str) {
                    throw exception_failure("GFitsTableShortCol::string() : Reference sum: "+sum_sht_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try ulong
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== S H O R T  1 0 =====
        //
    
        GFitsTableShortCol col_sht10;
            
        test_try("Short 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col_sht10 = GFitsTableShortCol("SHORT10", nrows,nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        short val_sht  = short(100.0*cos(0.1*float(i))*
                                cos(0.33*float(j)));
                        col_sht10(i,j) = val_sht;
                        sum_sht10     += val_sht;
                        sum_sht10_int += int(val_sht);
                        sum_sht10_str += ":"+str(val_sht);
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_sht = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_sht += col_sht10(i,j);
                }
                if (tot_sht != sum_sht10) {
                    throw exception_failure("GFitsTableShortCol::operator() - 10 : Reference sum: "+str(sum_sht10)+"Derived sum:"+str(tot_sht));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_sht10.real(i,j);
                }
                if (!dequal(tot_dbl, sum_sht10)) {
                    throw exception_failure("GFitsTableShortCol::real() - 10 : Reference sum: "+str(sum_sht10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_sht10.integer(i,j);
                }
                if (tot_int != sum_sht10_int) {
                    throw exception_failure("GFitsTableShortCol::integer() - 10 : Reference sum: "+str(sum_sht10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_sht10.string(i,j);
                }
                if (tot_str != sum_sht10_str) {
                    throw exception_failure("GFitsTableShortCol::string() - 10 : Reference sum: "+sum_sht10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try short 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col_sht);
            table.insert_column(99, col_sht10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        short       tot_sht = 0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;
        
        //
        // ===== U S H O R T =====
        //
        GFitsTableShortCol* col_sht = NULL;
        test_try("Short");
        try {
            // Get column
            col_sht = (GFitsTableShortCol*)&(*fits.table(1))["SHORT"];

            // Check table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_sht = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_sht += (*col_sht)(i);
                if (tot_sht != sum_sht) {
                    throw exception_failure("GFitsTableShortCol::operator() : Reference sum: "+str(sum_sht)+"  Derived sum:   "+str(tot_sht));
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_sht->real(i);
                if (!dequal(tot_dbl, double(sum_sht))) {
                    throw exception_failure("GFitsTableShortCol::real() : Reference sum: "+str(sum_sht)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_sht->integer(i);
                if (tot_int != sum_sht_int) {
                    throw exception_failure("GFitsTableShortCol::integer() : Reference sum: "+str(sum_sht_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_sht->string(i);
                if (tot_str != sum_sht_str) {
                    throw exception_failure("GFitsTableShortCol::string() : Reference sum: "+sum_sht_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== S H O R T  1 0 =====
        //
        GFitsTableShortCol* col_sht10 = NULL;
        test_try("Short 10");
        try{
            // Get column
            col_sht10 = (GFitsTableShortCol*)&(*fits.table(1))["SHORT10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_sht = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_sht += (*col_sht10)(i,j);
                }
                if (tot_sht != sum_sht10) {
                    throw exception_failure("GFitsTableShortCol::operator() - 10 : Reference sum: "+str(sum_sht10)+"Derived sum:"+str(tot_sht));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_sht10->real(i,j);
                }
                if (!dequal(tot_dbl, sum_sht10)) {
                    throw exception_failure("GFitsTableShortCol::real() - 10 : Reference sum: "+str(sum_sht10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_sht10->integer(i,j);
                }
                if (tot_int != sum_sht10_int) {
                    throw exception_failure("GFitsTableShortCol::integer() - 10 : Reference sum: "+str(sum_sht10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_sht10->string(i,j);
                }
                if (tot_str != sum_sht10_str) {
                    throw exception_failure("GFitsTableShortCol::string() - 10 : Reference sum: "+sum_sht10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try logical 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }


}

/***************************************************************************
 * @brief Test unsigned short FITS binary table
 ***************************************************************************/
 void TestGFits::test_bintable_ushort(void)
{

    // Set filename
    std::string filename = "test_bintable_ushort.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    unsigned short sum_sht       = 0;
    unsigned short sum_sht10     = 0;
    unsigned short sum_sht_int   = 0;
    unsigned short sum_sht10_int = 0;
    std::string    sum_sht_str;
    std::string    sum_sht10_str;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float          tot_flt = 0.0;
        double         tot_dbl = 0.0;
        int            tot_int = 0;
        long           tot_lng = 0;
        unsigned short tot_sht = 0;
        std::string    tot_str;
    
        //
        //
        // ===== U S H O R T =====
        //
    
        GFitsTableUShortCol col_sht;
    
        test_try("UShort");
        try{
             // Set string table
            test_try("Set table");
            try{
                col_sht = GFitsTableUShortCol("USHORT", nrows);
                for (int i = 0; i < nrows; ++i) {
                    unsigned short val_sht = (unsigned short)std::abs(1000.0 * cos(0.1*float(i)));
                    col_sht(i)    = val_sht;
                    sum_sht      += val_sht;
                    sum_sht_int  += int(val_sht);
                    sum_sht_str += ":"+str(val_sht);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_sht = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_sht += col_sht(i);
                if (tot_sht != sum_sht) {
                    throw exception_failure("GFitsTableUShortCol::operator() : Reference sum: "+str(sum_sht)+"  Derived sum:   "+str(tot_sht));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_sht.real(i);
                if (!dequal(tot_dbl, double(sum_sht))) {
                    throw exception_failure("GFitsTableUShortCol::real() : Reference sum: "+str(sum_sht)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_sht.integer(i);
                if (tot_int != sum_sht_int) {
                    throw exception_failure("GFitsTableUShortCol::integer() : Reference sum: "+str(sum_sht_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_sht.string(i);
                if (tot_str != sum_sht_str) {
                    throw exception_failure("GFitsTableUShortCol::string() : Reference sum: "+sum_sht_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try ulong
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== S H O R T  1 0 =====
        //
    
        GFitsTableUShortCol col_sht10;
            
        test_try("UShort 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col_sht10 = GFitsTableUShortCol("USHORT10", nrows,nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        unsigned short val_sht  =
                                (unsigned short)(std::abs(100.0*cos(0.1*float(i))*cos(0.33*float(j))));
                        col_sht10(i,j) = val_sht;
                        sum_sht10     += val_sht;
                        sum_sht10_int += int(val_sht);
                        sum_sht10_str += ":"+str(val_sht);
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_sht = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_sht += col_sht10(i,j);
                }
                if (tot_sht != sum_sht10) {
                    throw exception_failure("GFitsTableUShortCol::operator() - 10 : Reference sum: "+str(sum_sht10)+"Derived sum:"+str(tot_sht));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_sht10.real(i,j);
                }
                if (!dequal(tot_dbl, sum_sht10)) {
                    throw exception_failure("GFitsTableUShortCol::real() - 10 : Reference sum: "+str(sum_sht10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_sht10.integer(i,j);
                }
                if (tot_int != sum_sht10_int) {
                    throw exception_failure("GFitsTableUShortCol::integer() - 10 : Reference sum: "+str(sum_sht10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_sht10.string(i,j);
                }
                if (tot_str != sum_sht10_str) {
                    throw exception_failure("GFitsTableUShortCol::string() - 10 : Reference sum: "+sum_sht10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try longlong 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col_sht);
            table.insert_column(99, col_sht10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float          tot_flt = 0.0;
        unsigned short tot_sht = 0;
        double         tot_dbl = 0.0;
        int            tot_int = 0;
        std::string    tot_str;
        
        //
        // ===== U S H O R T =====
        //
        GFitsTableUShortCol* col_sht = NULL;
        test_try("UShort");
        try {
            // Get column
            col_sht = (GFitsTableUShortCol*)&(*fits.table(1))["USHORT"];

            // Check table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_sht = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_sht += (*col_sht)(i);
                if (tot_sht != sum_sht) {
                    throw exception_failure("GFitsTableUShortCol::operator() : Reference sum: "+str(sum_sht)+"  Derived sum:   "+str(tot_sht));
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_sht->real(i);
                if (!dequal(tot_dbl, double(sum_sht))) {
                    throw exception_failure("GFitsTableUShortCol::real() : Reference sum: "+str(sum_sht)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_sht->integer(i);
                if (tot_int != sum_sht_int) {
                    throw exception_failure("GFitsTableUShortCol::integer() : Reference sum: "+str(sum_sht_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_sht->string(i);
                if (tot_str != sum_sht_str) {
                    throw exception_failure("GFitsTableUShortCol::string() : Reference sum: "+sum_sht_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== U S H O R T  1 0 =====
        //
        GFitsTableUShortCol* col_sht10 = NULL;
        test_try("UShort 10");
        try{
            // Get column
            col_sht10 = (GFitsTableUShortCol*)&(*fits.table(1))["USHORT10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_sht = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_sht += (*col_sht10)(i,j);
                }
                if (tot_sht != sum_sht10) {
                    throw exception_failure("GFitsTableUShortCol::operator() - 10 : Reference sum: "+str(sum_sht10)+"Derived sum:"+str(tot_sht));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_sht10->real(i,j);
                }
                if (!dequal(tot_dbl, sum_sht10)) {
                    throw exception_failure("GFitsTableUShortCol::real() - 10 : Reference sum: "+str(sum_sht10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_sht10->integer(i,j);
                }
                if (tot_int != sum_sht10_int) {
                    throw exception_failure("GFitsTableUShortCol::integer() - 10 : Reference sum: "+str(sum_sht10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_sht10->string(i,j);
                }
                if (tot_str != sum_sht10_str) {
                    throw exception_failure("GFitsTableUShortCol::string() - 10 : Reference sum: "+sum_sht10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try logical 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }


}

/***************************************************************************
 * @brief Test long FITS binary table
 ***************************************************************************/

void TestGFits::test_bintable_long(void)
{

    // Set filename
    std::string filename = "test_bintable_long.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    long        sum_lng       = 0;
    long        sum_lng10     = 0;
    long        sum_lng_int   = 0;
    long        sum_lng10_int = 0;
    std::string sum_lng_str;
    std::string sum_lng10_str;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;
    
        //
        //
        // ===== L O N G  =====
        //
    
        GFitsTableLongCol col_lng;
    
        test_try("Long");
        try{
             // Set string table
            test_try("Set table");
            try{
                col_lng = GFitsTableLongCol("LONG", nrows);
                for (int i = 0; i < nrows; ++i) {
                    long val_lng  = long(100000.0 * cos(0.1*float(i)));
                    col_lng(i)    = val_lng;
                    sum_lng      += val_lng;
                    sum_lng_int  += int(val_lng);
                    sum_lng_str += ":"+str((int)val_lng);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_lng += col_lng(i);
                if (tot_lng != sum_lng) {
                    throw exception_failure("GFitsTableLongCol::operator() : Reference sum: "+str(sum_lng)+"  Derived sum:   "+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_lng.real(i);
                if (!dequal(tot_dbl, double(sum_lng))) {
                    throw exception_failure("GFitsTableLongCol::real() : Reference sum: "+str(sum_lng)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_lng.integer(i);
                if (tot_int != sum_lng_int) {
                    throw exception_failure("GFitsTableLongCol::integer() : Reference sum: "+str(sum_lng_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_lng.string(i);
                if (tot_str != sum_lng_str) {
                    throw exception_failure("GFitsTableLongCol::string() : Reference sum: "+sum_lng_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try ulong
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== L O N G  1 0 =====
        //
    
        GFitsTableLongCol col_lng10;
            
        test_try("Long 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col_lng10 = GFitsTableLongCol("LONG10", nrows,nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        long val_lng   = long(1000.0*cos(0.1*float(i))*
                                cos(0.33*float(j)));
                        col_lng10(i,j) = val_lng;
                        sum_lng10     += val_lng;
                        sum_lng10_int += int(val_lng);
                        sum_lng10_str += ":"+str((int)val_lng);
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_lng += col_lng10(i,j);
                }
                if (tot_lng != sum_lng10) {
                    throw exception_failure("GFitsTableLongCol::operator() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_lng10.real(i,j);
                }
                if (!dequal(tot_dbl, sum_lng10)) {
                    throw exception_failure("GFitsTableLongCol::real() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_lng10.integer(i,j);
                }
                if (tot_int != sum_lng10_int) {
                    throw exception_failure("GFitsTableLongCol::integer() - 10 : Reference sum: "+str(sum_lng10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_lng10.string(i,j);
                }
                if (tot_str != sum_lng10_str) {
                    throw exception_failure("GFitsTableLongCol::string() - 10 : Reference sum: "+sum_lng10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try longlong 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col_lng);
            table.insert_column(99, col_lng10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        std::string tot_str;
        
        //
        // ===== L O N G  =====
        //
        GFitsTableLongCol* col_lng = NULL;
        test_try("Long");
        try {
            // Get column
            col_lng = (GFitsTableLongCol*)&(*fits.table(1))["LONG"];

            // Check string table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_lng += (*col_lng)(i);
                if (tot_lng != sum_lng) {
                    throw exception_failure("GFitsTableLongCol::operator() : Reference sum: "+str(sum_lng)+"  Derived sum:   "+str(tot_lng));
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_lng->real(i);
                if (!dequal(tot_dbl, double(sum_lng))) {
                    throw exception_failure("GFitsTableLongCol::real() : Reference sum: "+str(sum_lng)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_lng->integer(i);
                if (tot_int != sum_lng_int) {
                    throw exception_failure("GFitsTableLongCol::integer() : Reference sum: "+str(sum_lng_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_lng->string(i);
                if (tot_str != sum_lng_str) {
                    throw exception_failure("GFitsTableLongCol::string() : Reference sum: "+sum_lng_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== L O N G 1 0 =====
        //
        GFitsTableLongCol* col_lng10 = NULL;
        test_try("Long 10");
        try{
            // Get column
            col_lng10 = (GFitsTableLongCol*)&(*fits.table(1))["LONG10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_lng += (*col_lng10)(i,j);
                }
                if (tot_lng != sum_lng10) {
                    throw exception_failure("GFitsTableLongCol::operator() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_lng10->real(i,j);
                }
                if (!dequal(tot_dbl, sum_lng10)) {
                    throw exception_failure("GFitsTableLongCol::real() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_lng10->integer(i,j);
                }
                if (tot_int != sum_lng10_int) {
                    throw exception_failure("GFitsTableLongCol::integer() - 10 : Reference sum: "+str(sum_lng10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_lng10->string(i,j);
                }
                if (tot_str != sum_lng10_str) {
                    throw exception_failure("GFitsTableLongCol::string() - 10 : Reference sum: "+sum_lng10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try logical 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }
    

}

/***************************************************************************
 * @brief Test long long FITS binary table
 ***************************************************************************/
 
void TestGFits::test_bintable_longlong(void)
{

    // Set filename
    std::string filename = "test_bintable_longlong.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    long long   sum_lng       = 0;
    long long   sum_lng10     = 0;
    long long   sum_lng_int   = 0;
    long long   sum_lng10_int = 0;
    std::string sum_lng_str;
    std::string sum_lng10_str;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long long   tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;
    
        //
        //
        // ===== L O N G  L O N G =====
        //
    
        GFitsTableLongLongCol col_lng;
    
        test_try("LongLong");
        try{
             // Set string table
            test_try("Set table");
            try{
                col_lng = GFitsTableLongLongCol("LONGLONG", nrows);
                for (int i = 0; i < nrows; ++i) {
                    long val_lng  = long(100000.0 * cos(0.1*float(i)));
                    col_lng(i)    = val_lng;
                    sum_lng      += val_lng;
                    sum_lng_int  += int(val_lng);
                    sum_lng_str += ":"+str((int)val_lng);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_lng += col_lng(i);
                if (tot_lng != sum_lng) {
                    throw exception_failure("GFitsTableLongLongCol::operator() : Reference sum: "+str(sum_lng)+"  Derived sum:   "+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_lng.real(i);
                if (!dequal(tot_dbl, double(sum_lng))) {
                    throw exception_failure("GFitsTableLongLongCol::real() : Reference sum: "+str(sum_lng)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_lng.integer(i);
                if (tot_int != sum_lng_int) {
                    throw exception_failure("GFitsTableLongLongCol::integer() : Reference sum: "+str(sum_lng_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_lng.string(i);
                if (tot_str != sum_lng_str) {
                    throw exception_failure("GFitsTableLongLongCol::string() : Reference sum: "+sum_lng_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try ulong
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== L O N G  L O N G 1 0 =====
        //
    
        GFitsTableLongLongCol col_lng10;
            
        test_try("LongLong 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col_lng10 = GFitsTableLongLongCol("LONGLONG10", nrows,nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        long val_lng   = long(1000.0*cos(0.1*float(i))*
                                cos(0.33*float(j)));
                        col_lng10(i,j) = val_lng;
                        sum_lng10     += val_lng;
                        sum_lng10_int += int(val_lng);
                        sum_lng10_str += ":"+str((int)val_lng);
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_lng += col_lng10(i,j);
                }
                if (tot_lng != sum_lng10) {
                    throw exception_failure("GFitsTableLongLongCol::operator() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_lng10.real(i,j);
                }
                if (!dequal(tot_dbl, sum_lng10)) {
                    throw exception_failure("GFitsTableLongLongCol::real() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_lng10.integer(i,j);
                }
                if (tot_int != sum_lng10_int) {
                    throw exception_failure("GFitsTableLongLongCol::integer() - 10 : Reference sum: "+str(sum_lng10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_lng10.string(i,j);
                }
                if (tot_str != sum_lng10_str) {
                    throw exception_failure("GFitsTableLongLongCol::string() - 10 : Reference sum: "+sum_lng10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try longlong 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col_lng);
            table.insert_column(99, col_lng10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long long   tot_lng = 0;
        std::string tot_str;
        
        //
        // ===== L O N G L O N G =====
        //
        GFitsTableLongLongCol* col_lng = NULL;
        test_try("LongLong");
        try {
            // Get column
            col_lng = (GFitsTableLongLongCol*)&(*fits.table(1))["LONGLONG"];

            // Check string table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_lng += (*col_lng)(i);
                if (tot_lng != sum_lng) {
                    throw exception_failure("GFitsTableLongLongCol::operator() : Reference sum: "+str(sum_lng)+"  Derived sum:   "+str(tot_lng));
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_lng->real(i);
                if (!dequal(tot_dbl, double(sum_lng))) {
                    throw exception_failure("GFitsTableLongLongCol::real() : Reference sum: "+str(sum_lng)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_lng->integer(i);
                if (tot_int != sum_lng_int) {
                    throw exception_failure("GFitsTableLongLongCol::integer() : Reference sum: "+str(sum_lng_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_lng->string(i);
                if (tot_str != sum_lng_str) {
                    throw exception_failure("GFitsTableLongLongCol::string() : Reference sum: "+sum_lng_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== L O N G  L O N G 1 0 =====
        //
        GFitsTableLongLongCol* col_lng10 = NULL;
        test_try("LongLong 10");
        try{
            // Get column
            col_lng10 = (GFitsTableLongLongCol*)&(*fits.table(1))["LONGLONG10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_lng += (*col_lng10)(i,j);
                }
                if (tot_lng != sum_lng10) {
                    throw exception_failure("GFitsTableLongLongCol::operator() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_lng10->real(i,j);
                }
                if (!dequal(tot_dbl, sum_lng10)) {
                    throw exception_failure("GFitsTableLongLongCol::real() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_lng10->integer(i,j);
                }
                if (tot_int != sum_lng10_int) {
                    throw exception_failure("GFitsTableLongLongCol::integer() - 10 : Reference sum: "+str(sum_lng10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_lng10->string(i,j);
                }
                if (tot_str != sum_lng10_str) {
                    throw exception_failure("GFitsTableLongLongCol::string() - 10 : Reference sum: "+sum_lng10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try logical 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }
    

}

/***************************************************************************
 * @brief Test unsigned long FITS binary table
 ***************************************************************************/

void TestGFits::test_bintable_ulong(void)
{

    // Set filename
    std::string filename = "test_bintable_ulong.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    unsigned long sum_lng       = 0;
    unsigned long sum_lng10     = 0;
    unsigned long sum_lng_int   = 0;
    unsigned long sum_lng10_int = 0;
    std::string   sum_lng_str;
    std::string   sum_lng10_str;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float         tot_flt = 0.0;
        double        tot_dbl = 0.0;
        int           tot_int = 0;
        unsigned long tot_lng = 0;
        short         tot_sht = 0;
        std::string   tot_str;
    
        //
        // ===== U L O N G=====
        //
    
        GFitsTableULongCol col_lng;
    
        test_try("Ulong");
        try{
             // Set string table
            test_try("Set table");
            try{
                col_lng = GFitsTableULongCol("ULONG", nrows);
                for (int i = 0; i < nrows; ++i) {
                    unsigned long val_lng  = (unsigned long)std::abs(100000.0 * cos(0.1*float(i)));
                    col_lng(i)    = val_lng;
                    sum_lng      += val_lng;
                    sum_lng_int  += int(val_lng);
                    sum_lng_str += ":"+str((int)val_lng);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_lng += col_lng(i);
                if (tot_lng != sum_lng) {
                    throw exception_failure("GFitsTableULongCol::operator() : Reference sum: "+str(sum_lng)+"  Derived sum:   "+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_lng.real(i);
                if (!dequal(tot_dbl, double(sum_lng))) {
                    throw exception_failure("GFitsTableULongCol::real() : Reference sum: "+str(sum_lng)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_lng.integer(i);
                if (tot_int != sum_lng_int) {
                    throw exception_failure("GFitsTableULongCol::integer() : Reference sum: "+str(sum_lng_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_lng.string(i);
                if (tot_str != sum_lng_str) {
                    throw exception_failure("GFitsTableULongCol::string() : Reference sum: "+sum_lng_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try ulong
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== U L O N G 1 0 =====
        //
    
        GFitsTableULongCol col_lng10;
            
        test_try("Ulong 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col_lng10 = GFitsTableULongCol("ULONG10", nrows,nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        unsigned long val_lng = (unsigned long)std::abs(1000.0*cos(0.1*float(i))*
                                cos(0.33*float(j)));
                        col_lng10(i,j) = val_lng;
                        sum_lng10     += val_lng;
                        sum_lng10_int += int(val_lng);
                        sum_lng10_str += ":"+str((int)val_lng);
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_lng += col_lng10(i,j);
                }
                if (tot_lng != sum_lng10) {
                    throw exception_failure("GFitsTableULongCol::operator() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_lng10.real(i,j);
                }
                if (!dequal(tot_dbl, sum_lng10)) {
                    throw exception_failure("GFitsTableULongCol::real() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_lng10.integer(i,j);
                }
                if (tot_int != sum_lng10_int) {
                    throw exception_failure("GFitsTableULongCol::integer() - 10 : Reference sum: "+str(sum_lng10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_lng10.string(i,j);
                }
                if (tot_str != sum_lng10_str) {
                    throw exception_failure("GFitsTableULongCol::string() - 10 : Reference sum: "+sum_lng10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try ulong 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col_lng);
            table.insert_column(99, col_lng10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float         tot_flt = 0.0;
        double        tot_dbl = 0.0;
        int           tot_int = 0;
        unsigned long tot_lng = 0;
        std::string   tot_str;
        
        //
        // ===== U L O N G =====
        //
        GFitsTableULongCol* col_lng = NULL;
        test_try("ULong");
        try {
            // Get column
            col_lng = (GFitsTableULongCol*)&(*fits.table(1))["ULONG"];

            // Check string table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_lng += (*col_lng)(i);
                if (tot_lng != sum_lng) {
                    throw exception_failure("GFitsTableULongCol::operator() : Reference sum: "+str(sum_lng)+"  Derived sum:   "+str(tot_lng));
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_lng->real(i);
                if (!dequal(tot_dbl, double(sum_lng))) {
                    throw exception_failure("GFitsTableULongCol::real() : Reference sum: "+str(sum_lng)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_lng->integer(i);
                if (tot_int != sum_lng_int) {
                    throw exception_failure("GFitsTableULongCol::integer() : Reference sum: "+str(sum_lng_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_lng->string(i);
                if (tot_str != sum_lng_str) {
                    throw exception_failure("GFitsTableULongCol::string() : Reference sum: "+sum_lng_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== U L O N G   1 0 =====
        //
        GFitsTableULongCol* col_lng10 = NULL;
        test_try("ULong 10");
        try{
            // Get column
            col_lng10 = (GFitsTableULongCol*)&(*fits.table(1))["ULONG10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_lng = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_lng += (*col_lng10)(i,j);
                }
                if (tot_lng != sum_lng10) {
                    throw exception_failure("GFitsTableULongCol::operator() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_lng));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_lng10->real(i,j);
                }
                if (!dequal(tot_dbl, sum_lng10)) {
                    throw exception_failure("GFitsTableULongCol::real() - 10 : Reference sum: "+str(sum_lng10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_lng10->integer(i,j);
                }
                if (tot_int != sum_lng10_int) {
                    throw exception_failure("GFitsTableULongCol::integer() - 10 : Reference sum: "+str(sum_lng10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_lng10->string(i,j);
                }
                if (tot_str != sum_lng10_str) {
                    throw exception_failure("GFitsTableULongCol::string() - 10 : Reference sum: "+sum_lng10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try logical 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }
    

}

/***************************************************************************
 * @brief Test string FITS binary table
 ***************************************************************************/
 void TestGFits::test_bintable_string(void)
{
    
    // Set filename
    std::string filename = "test_bintable_string.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    std::string sum_str;
    std::string sum_str10;
    double      sum_str_dbl   = 0;
    double      sum_str10_dbl = 0;
    double      sum_str_int   = 0;
    double      sum_str10_int = 0;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        std::string tot_str;
    
        //
        // ===== S T R I N G =====
        //
    
        GFitsTableStringCol col_str;
    
        test_try("String");
        try{
             // Set string table
            test_try("Set table");
            try{
                col_str = GFitsTableStringCol("STRING", nrows, 20);
                for (int i = 0; i < nrows; ++i) {
                    double val_dbl = 100.0 * cos(0.1*double(i));
                    col_str(i)   = str(val_dbl);
                    sum_str     += ":" + str(val_dbl);
                    sum_str_dbl += val_dbl;
                    sum_str_int += int(val_dbl);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":" + col_str(i);
                if (tot_str != sum_str) {
                    throw exception_failure("GFitsTableStringCol::operator() : Reference sum: "+sum_str+"  Derived sum:   "+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_str.real(i);
                if (!equal(tot_dbl, sum_str_dbl, 1.0e-2)) {
                    throw exception_failure("GFitsTableStringCol::real() : Reference sum: "+str(sum_str_dbl)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_str.integer(i);
                if (tot_int != sum_str_int) {
                    throw exception_failure("GFitsTableStringCol::integer() : Reference sum: "+str(sum_str_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_str.string(i);
                if (tot_str != sum_str) {
                    throw exception_failure("GFitsTableBoolCol::string() : Reference sum: "+sum_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try logical
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== S T R I N G   1 0 =====
        //
    
        GFitsTableStringCol col_str10;
            
        test_try("String 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col_str10 = GFitsTableStringCol("STRING10", nrows, 20, nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        double val_dbl = 100.0 * cos(0.1*double(i)) * cos(0.3*double(j));
                        col_str10(i,j)  = str(val_dbl);
                        sum_str10      += ":" + str(val_dbl);
                        sum_str10_dbl  += val_dbl;
                        sum_str10_int  += int(val_dbl);
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":" + col_str10(i,j);
                }
                if (tot_str != sum_str10) {
                    throw exception_failure("GFitsTableStringCol::operator() - 10 : Reference sum: "+sum_str10+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_str10.real(i,j);
                }
                if (!equal(tot_dbl, sum_str10_dbl, 1.0e-2)) {
                    throw exception_failure("GFitsTableStringCol::real() - 10 : Reference sum: "+str(sum_str10_dbl)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_str10.integer(i,j);
                }
                if (tot_int != sum_str10_int) {
                    throw exception_failure("GFitsTableStringCol::integer() - 10 : Reference sum: "+str(sum_str10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_str10.string(i,j);
                }
                if (tot_str != sum_str10) {
                    throw exception_failure("GFitsTableStringCol::string() - 10 : Reference sum: "+sum_str10+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try string 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col_str);
            table.insert_column(0, col_str10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;
        
        //
        // ===== L O G I C A L =====
        //
        GFitsTableStringCol* col_str = NULL;
        test_try("String");
        try {
            // Get column
            col_str = (GFitsTableStringCol*)&(*fits.table(1))["STRING"];

            // Check string table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":" + (*col_str)(i);
                if (tot_str != sum_str) {
                    throw exception_failure("GFitsTableStringCol::operator() : Reference sum: "+sum_str+"  Derived sum:   "+tot_str);
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col_str->real(i);
                if (!equal(tot_dbl, sum_str_dbl, 1.0e-2)) {
                    throw exception_failure("GFitsTableStringCol::real() : Reference sum: "+str(sum_str_dbl)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check string table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col_str->integer(i);
                if (tot_int != sum_str_int) {
                    throw exception_failure("GFitsTableStringCol::integer() : Reference sum: "+str(sum_str_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col_str->string(i);
                if (tot_str != sum_str) {
                    throw exception_failure("GFitsTableStringCol::string() : Reference sum: "+sum_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== S T R I N G   1 0 =====
        //
        GFitsTableStringCol* col_str10 = NULL;
        test_try("String 10");
        try{
            // Get column
            col_str10 = (GFitsTableStringCol*)&(*fits.table(1))["STRING10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":" + (*col_str10)(i,j);
                }
                if (tot_str != sum_str10) {
                    throw exception_failure("GFitsTableStringCol::operator() - 10 : Reference sum: "+sum_str10+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check string table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col_str10->real(i,j);
                }
                if (!equal(tot_dbl, sum_str10_dbl, 1.0e-2)) {
                    throw exception_failure("GFitsTableStringCol::real() - 10 : Reference sum: "+str(sum_str10_dbl)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col_str10->integer(i,j);
                }
                if (tot_int != sum_str10_int) {
                    throw exception_failure("GFitsTableStringCol::integer() - 10 : Reference sum: "+str(sum_str10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col_str10->string(i,j);
                }
                if (tot_str != sum_str10) {
                    throw exception_failure("GFitsTableStringCol::string() - 10 : Reference sum: "+sum_str10+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try logical 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }
    

}

/***************************************************************************
 * @brief Test boolean FITS binary table
 ***************************************************************************/
void TestGFits::test_bintable_logical(void)
{
    
    // Set filename
    std::string filename = "test_bintable_logical.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    int         sum       = 0;
    int         sum10     = 0;
    int         sum_int   = 0;
    int         sum10_int = 0;
    std::string sum_str;
    std::string sum10_str;
   
    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        std::string tot_str;
    
        //
        // ===== L O G I C A L =====
        //
    
        GFitsTableBoolCol col;
    
        test_try("Logical");
        try{
            // Set table
            test_try("Set table");
            try{
                col = GFitsTableBoolCol("LOGICAL", nrows);
                for (int i = 0; i < nrows; ++i) {
                    bool val  = (i % 2) ? true : false;
                    col(i)    = val;
                    sum      += val;
                    sum_int  += int(val);
                    if (val)
                        sum_str += ":T";
                    else
                        sum_str += ":F";
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col(i);
                if (tot_int != sum) {
                    throw exception_failure("GFitsTableBoolCol::operator() : Reference sum: "+str(sum)+"  Derived sum:   "+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col.real(i);
                if (!dequal(tot_dbl, double(sum))) {
                    throw exception_failure("GFitsTableBoolCol::real() : Reference sum: "+str(sum)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col.integer(i);
                if (tot_int != sum_int) {
                    throw exception_failure("GFitsTableBoolCol::integer() : Reference sum: "+str(sum_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col.string(i);
                if (tot_str != sum_str) {
                    throw exception_failure("GFitsTableBoolCol::string() : Reference sum: "+sum_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try logical
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
        //
        // ===== L O G I C A L  1 0 =====
        //
    
        GFitsTableBoolCol col10;
            
        test_try("Logical 10");
        try{
    
                // Set table
            test_try("Set table");
            try{
                col10 = GFitsTableBoolCol("LOGICAL10", nrows, nvec);
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j) {
                        char val  = (i % 2) * (j % 2);
                        col10(i,j) = val;
                        sum10     += val;
                        sum10_int += int(val);
                        if (val)
                            sum10_str += ":T";
                        else
                            sum10_str += ":F";
                    }
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col10(i,j);
                }
                if (tot_int != sum10) {
                    throw exception_failure("GFitsTableBoolCol::operator() - 10 : Reference sum: "+str(sum10)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col10.real(i,j);
                }
                if (!dequal(tot_dbl, sum10)) {
                    throw exception_failure("GFitsTableBoolCol::real() - 10 : Reference sum: "+str(sum10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col10.integer(i,j);
                }
                if (tot_int != sum10_int) {
                    throw exception_failure("GFitsTableBoolCol::integer() - 10 : Reference sum: "+str(sum10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col10.string(i,j);
                }
                if (tot_str != sum10_str) {
                    throw exception_failure("GFitsTableBoolCol::string() - 10 : Reference sum: "+sum10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // end try logical 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
    
        //
            // ===== W R I T E   T A B L E =====
        //
        test_try("Write Table");
        try{
                // Build binary table
            GFitsBinTable table = GFitsBinTable(nrows);
            table.append_column(col);
            table.insert_column(0, col10);
        
                // Append to file
            fits.append(table);
        
                // Save FITS file
            fits.save();
        
            test_try_success();
        }
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
            
        test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    
    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        short       tot_sht = 0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;
        
        //
        // ===== L O G I C A L =====
        //
        GFitsTableBoolCol* col = NULL;
        test_try("Logical");
        try {
            // Get column
            col = (GFitsTableBoolCol*)&(*fits.table(1))["LOGICAL"];

            // Check table (operator access)
            test_try("Check table (operator access)");
            try {
                
                tot_int = 0;
      
                for (int i = 0; i < nrows; ++i)
                    tot_int += (*col)(i);
                if (tot_int != sum) {
                    throw exception_failure("GFitsTableBoolCol::operator() : Reference sum: "+str(sum)+"  Derived sum:   "+str(tot_int));
                }
                
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col->real(i);
                if (!dequal(tot_dbl, double(sum))) {
                    throw exception_failure("GFitsTableBoolCol::real() : Reference sum: "+str(sum)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col->integer(i);
                if (tot_int != sum_int) {
                    throw exception_failure("GFitsTableBoolCol::integer() : Reference sum: "+str(sum_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col->string(i);
                if (tot_str != sum_str) {
                    throw exception_failure("GFitsTableBoolCol::string() : Reference sum: "+sum_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        }// End try logical
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        //
        // ===== L O G I C A L  1 0 =====
        //
        GFitsTableBoolCol* col10 = NULL;
        test_try("Logical 10");
        try{
            // Get column
            col10 = (GFitsTableBoolCol*)&(*fits.table(1))["LOGICAL10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += (*col10)(i,j);
                }
                if (tot_int != sum10) {
                    throw exception_failure("GFitsTableBoolCol::operator() - 10 : Reference sum: "+str(sum10)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_dbl += col10->real(i,j);
                }
                if (!dequal(tot_dbl, sum10)) {
                    throw exception_failure("GFitsTableBoolCol::real() - 10 : Reference sum: "+str(sum10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_int += col10->integer(i,j);
                }
                if (tot_int != sum10_int) {
                    throw exception_failure("GFitsTableBoolCol::integer() - 10 : Reference sum: "+str(sum10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col10->string(i,j);
                }
                if (tot_str != sum10_str) {
                    throw exception_failure("GFitsTableBoolCol::string() - 10 : Reference sum: "+sum10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try logical 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }
    

}

/***************************************************************************
 * @brief Test bit FITS binary table
 ***************************************************************************/
void TestGFits::test_bintable_bit(void)
{
    // Set filename
    std::string filename = "test_bintable_bit.fits";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    int rc = system(cmd.c_str());

    // Allocate reference sums
    int         sum       = 0;
    int         sum10     = 0;
    int         sum_int   = 0;
    int         sum10_int = 0;
    std::string sum_str;
    std::string sum10_str;

    // Build tables
    test_try("Build tables");
    try {
        // Create FITS file
        GFits fits;
        fits.open(filename, true);
    
        // Set number of rows
        int nrows = 10;
        int nvec  = 10;
    
        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        std::string tot_str;
    
        //
        // ===== B I T =====
        //
    
        GFitsTableBitCol col;
    
        test_try("BIT");
        try{
            // Set table
            test_try("Set table");
            try{
                col = GFitsTableBitCol("BIT", nrows);
                for (int i = 0; i < nrows; ++i) {
                    char val  = (i % 2) ? 1 : 0;
                    col(i)    = val;
                    sum      += val;
                    sum_int  += int(val);
                    if (val)
                        sum_str += ":T";
                    else
                        sum_str += ":F";
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col(i);
                if (tot_int != sum) {
                    throw exception_failure("GFitsTableBitCol::operator() : Reference sum: "+str(sum)+"  Derived sum:   "+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col.real(i);
                if (!dequal(tot_dbl, double(sum))) {
                    throw exception_failure("GFitsTableBitCol::real() : Reference sum: "+str(sum)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (int access)
            test_try("Check table (int access)");
            try{
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += col.integer(i);
                if (tot_int != sum_int) {
                    throw exception_failure("GFitsTableBitCol::integer() : Reference sum: "+str(sum_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i)
                    tot_str += ":"+col.string(i);
                if (tot_str != sum_str) {
                    throw exception_failure("GFitsTableBitCol::string() : Reference sum: "+sum_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        } // End try bit
        catch(std::exception &e)
        { 
            test_try_failure(e);
        }
    
            //
            // ===== B I T  1 0 =====
            //
    
            GFitsTableBitCol col10;
            
            test_try("Bit 10");
            try{
    
                // Set table
                test_try("Set table");
                try{
                    col10 = GFitsTableBitCol("BIT10", nrows, nvec);
                    for (int i = 0; i < nrows; ++i) {
                        for (int j = 0; j < nvec; ++j) {
                            char val  = (i % 2) * (j % 2);
                            col10(i,j) = val;
                            sum10     += val;
                            sum10_int += int(val);
                            if (val)
                                sum10_str += ":T";
                            else
                                sum10_str += ":F";
                        }
                    }
                    test_try_success();
                }
                catch(std::exception &e)
                {
                    test_try_failure(e);
                }
    
                // Check table (operator access)
                test_try("Check table (operator access)");
                try{
                    tot_int = 0;
                    for (int i = 0; i < nrows; ++i) {
                        for (int j = 0; j < nvec; ++j)
                            tot_int += col10(i,j);
                    }
                    if (tot_int != sum10) {
                        throw exception_failure("GFitsTableBitCol::operator() - 10 : Reference sum: "+str(sum10)+"Derived sum:"+str(tot_int));
                    }
                    test_try_success();
                }
                catch(std::exception &e)
                {
                    test_try_failure(e);
                }
        
                // Check table (real access)
                test_try("Check table (real access)");
                try{
                    tot_dbl = 0.0;
                    for (int i = 0; i < nrows; ++i) {
                        for (int j = 0; j < nvec; ++j)
                            tot_dbl += col10.real(i,j);
                    }
                    if (!dequal(tot_dbl, sum10)) {
                        throw exception_failure("GFitsTableBitCol::real() - 10 : Reference sum: "+str(sum10)+"Derived sum:"+str(tot_dbl));
                    }
                    test_try_success();
                }
                catch(std::exception &e)
                {
                    test_try_failure(e);
                }
        
                // Check table (int access)
                test_try("Check table (int access)");
                try{
                    tot_int = 0;
                    for (int i = 0; i < nrows; ++i) {
                        for (int j = 0; j < nvec; ++j)
                            tot_int += col10.integer(i,j);
                    }
                    if (tot_int != sum10_int) {
                        throw exception_failure("GFitsTableBitCol::integer() - 10 : Reference sum: "+str(sum10_int)+"Derived sum:"+str(tot_int));
                    }
                    test_try_success();
                }
                catch(std::exception &e)
                {
                    test_try_failure(e);
                }
        
                // Check table (string access)
                test_try("Check table (string access)");
                try{
                    tot_str.clear();
                    for (int i = 0; i < nrows; ++i) {
                        for (int j = 0; j < nvec; ++j)
                            tot_str += ":"+col10.string(i,j);
                    }
                    if (tot_str != sum10_str) {
                        throw exception_failure("GFitsTableBitCol::string() - 10 : Reference sum: "+sum10_str+"Derived sum:"+tot_str);
                    }
                    test_try_success();
                }
                catch(std::exception &e)
                {
                    test_try_failure(e);
                }
                test_try_success();
            } // end try bit 10
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
            //
            // ===== W R I T E   T A B L E =====
            //
            test_try("Write Table");
            try{
                // Build binary table
                GFitsBinTable table = GFitsBinTable(nrows);
                table.append_column(col);
                table.insert_column(0, col10);
        
                // Append to file
                fits.append(table);
        
                // Save FITS file
                fits.save();
        
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
    } //end try build table
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================

    test_try("Read Table");
    try{
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        short       tot_sht = 0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;

        //
        // ===== B I T =====
        //
        GFitsTableBitCol* col = NULL;
        test_try("Bit");
        try {
            // Get column
            col = (GFitsTableBitCol*)&(*fits.table(1))["BIT"];

            // Check table (operator access)
            test_try("Check table (operator access)");
            try {
                tot_int = 0;
                for (int i = 0; i < nrows; ++i)
                    tot_int += (*col)(i);
                if (tot_int != sum) {
                    throw exception_failure("GFitsTableBitCol::operator() : Reference sum: "+str(sum)+"  Derived sum:   "+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
    
                    // Check table (real access)
            test_try("Check table (real access)");
            try{
                tot_dbl = 0.0;
                for (int i = 0; i < nrows; ++i)
                    tot_dbl += col->real(i);
                if (!dequal(tot_dbl, double(sum))) {
                    throw exception_failure("GFitsTableBitCol::real() : Reference sum: "+str(sum)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                    tot_int = 0;
                    for (int i = 0; i < nrows; ++i)
                        tot_int += col->integer(i);
                    if (tot_int != sum_int) {
                    throw exception_failure("GFitsTableBitCol::integer() : Reference sum: "+str(sum_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                    tot_str.clear();
                    for (int i = 0; i < nrows; ++i)
                        tot_str += ":"+col->string(i);
                    if (tot_str != sum_str) {
                    throw exception_failure("GFitsTableBitCol::string() : Reference sum: "+sum_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            test_try_success();
        }// End try bit
        catch(std::exception &e)
        {
            test_try_failure(e);
        }

        //
        // ===== B I T  1 0 =====
        //
        GFitsTableBitCol* col10 = NULL;
        test_try("Bit 10");
        try{
            // Get column
            col10 = (GFitsTableBitCol*)&(*fits.table(1))["BIT10"];
            // Check table (operator access)
            test_try("Check table (operator access)");
            try{
                    tot_int = 0;
                    for (int i = 0; i < nrows; ++i) {
                        for (int j = 0; j < nvec; ++j)
                            tot_int += (*col10)(i,j);
                    }
                    if (tot_int != sum10) {
                    throw exception_failure("GFitsTableBitCol::operator() - 10 : Reference sum: "+str(sum10)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (real access)
            test_try("Check table (real access)");
            try{
                    tot_dbl = 0.0;
                    for (int i = 0; i < nrows; ++i) {
                        for (int j = 0; j < nvec; ++j)
                            tot_dbl += col10->real(i,j);
                    }
                    if (!dequal(tot_dbl, sum10)) {
                    throw exception_failure("GFitsTableBitCol::real() - 10 : Reference sum: "+str(sum10)+"Derived sum:"+str(tot_dbl));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (int access)
            test_try("Check table (int access)");
            try{
                    tot_int = 0;
                    for (int i = 0; i < nrows; ++i) {
                        for (int j = 0; j < nvec; ++j)
                            tot_int += col10->integer(i,j);
                    }
                    if (tot_int != sum10_int) {
                    throw exception_failure("GFitsTableBitCol::integer() - 10 : Reference sum: "+str(sum10_int)+"Derived sum:"+str(tot_int));
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
        
                // Check table (string access)
            test_try("Check table (string access)");
            try{
                tot_str.clear();
                for (int i = 0; i < nrows; ++i) {
                    for (int j = 0; j < nvec; ++j)
                        tot_str += ":"+col10->string(i,j);
                }
                if (tot_str != sum10_str) {
                    throw exception_failure("GFitsTableBitCol::string() - 10 : Reference sum: "+sum10_str+"Derived sum:"+tot_str);
                }
                test_try_success();
            }
            catch(std::exception &e)
            {
                test_try_failure(e);
            }
            
            test_try_success();
        } //end try bit 10
        catch(std::exception &e)
        {
            test_try_failure(e);
        }
        
        test_try_success();
    } //end try read tabe
    catch(std::exception &e)
    {
        test_try_failure(e);
    }

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuites("GFits");

    bool was_successful=true;

    //Create a test suite
    TestGFits test;

    //Append to the container
    testsuites.append(test);

    //Run
    was_successful=testsuites.run();

    //save xml report
    testsuites.save("reports/GFits.xml");

    // Return
    return was_successful ? 0:1;
}
