/***************************************************************************
 *                   test_GFits.cpp - test FITS classes                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file test_GFits.cpp
 * @brief Implementation of unit tests for FITS classes
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include <cstdlib>        // for system
#include "GTools.hpp"
#include "test_GFits.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Test macros ________________________________________________________ */
#define TEST_1D_ACCESS(NX) \
    for (int ix = 0, i = 0; ix < NX; ++ix, ++i) { \
        test_value(image1(ix), pixels[i], 1e-10, \
                   "Test pixel access operator for pixel "+gammalib::str(i)); \
        test_value(image1.at(ix), pixels[i], 1e-10, \
                   "Test at() method for pixel "+gammalib::str(i)); \
        test_value(image1.pixel(ix), pixels[i], 1e-10, \
                   "Test pixel() method for pixel "+gammalib::str(i)); \
    }
#define TEST_2D_ACCESS(NX,NY) \
    for (int iy = 0, i = 0; iy < NY; ++iy) { \
        for (int ix = 0; ix < NX; ++ix, ++i) { \
            test_value(image2(ix,iy), pixels[i], 1e-10, \
                       "Test pixel access operator for pixel "+gammalib::str(i)); \
            test_value(image2.at(ix,iy), pixels[i], 1e-10, \
                       "Test at() method for pixel "+gammalib::str(i)); \
            test_value(image2.pixel(ix,iy), pixels[i], 1e-10, \
                       "Test pixel() method for pixel "+gammalib::str(i)); \
        } \
    }
#define TEST_3D_ACCESS(NX,NY,NZ) \
    for (int iz = 0, i = 0; iz < NZ; ++iz) { \
        for (int iy = 0; iy < NY; ++iy) { \
            for (int ix = 0; ix < NX; ++ix, ++i) { \
                test_value(image3(ix,iy,iz), pixels[i], 1e-10, \
                           "Test pixel access operator for pixel "+gammalib::str(i)); \
                test_value(image3.at(ix,iy,iz), pixels[i], 1e-10, \
                           "Test at() method for pixel "+gammalib::str(i)); \
                test_value(image3.pixel(ix,iy,iz), pixels[i], 1e-10, \
                           "Test pixel() method for pixel "+gammalib::str(i)); \
            } \
        } \
    }
#define TEST_4D_ACCESS(NX,NY,NZ,NT) \
    for (int it = 0, i = 0; it < NT; ++it) { \
        for (int iz = 0; iz < NZ; ++iz) { \
            for (int iy = 0; iy < NY; ++iy) { \
                for (int ix = 0; ix < NX; ++ix, ++i) { \
                    test_value(image4(ix,iy,iz, it), pixels[i], 1e-10, \
                               "Test pixel access operator for pixel "+gammalib::str(i)); \
                    test_value(image4.at(ix,iy,iz, it), pixels[i], 1e-10, \
                               "Test at() method for pixel "+gammalib::str(i)); \
                    test_value(image4.pixel(ix,iy,iz, it), pixels[i], 1e-10, \
                               "Test pixel() method for pixel "+gammalib::str(i)); \
                } \
            } \
        } \
    }
#define TEST_4D_ACCESS_IO(NX,NY,NZ,NT) \
    for (int it = 0, i = 0; it < NT; ++it) { \
        for (int iz = 0; iz < NZ; ++iz) { \
            for (int iy = 0; iy < NY; ++iy) { \
                for (int ix = 0; ix < NX; ++ix, ++i) { \
                    test_value(ptr->pixel(ix,iy,iz, it), pixels[i], 1e-10, \
                               "Test pixel() method for pixel "+gammalib::str(i)); \
                } \
            } \
        } \
    }
#define TEST_TABLE1 \
    for (int i = 0; i < nrows; ++i) { \
        col1(i) = (double(i)*3.57+1.29); \
    } \
    for (int i = 0; i < nrows; ++i) { \
        double val = (double(i)*3.57+1.29); \
        test_value(col1(i), val, 1e-6, \
                   "Test access operator row "+gammalib::str(i)); \
        test_value(col1.real(i), val, 1e-6, \
                   "Test real() method for row "+gammalib::str(i)); \
        test_value(col1.integer(i), int(val), 1e-6, \
                   "Test integer() method for row "+gammalib::str(i)); \
        test_assert((col1.string(i) == gammalib::str(val)), \
                    "Test string() method for row "+gammalib::str(i), \
                    col1.string(i)+" is not "+gammalib::str(val)); \
    }
#define TEST_TABLE1_INT \
    for (int i = 0; i < nrows; ++i) { \
        col1(i) = i*3+1; \
    } \
    for (int i = 0; i < nrows; ++i) { \
        int val = i*3+1; \
        test_assert((col1(i) == val), \
                    "Test access operator row "+gammalib::str(i), \
                    gammalib::str(col1(i))+" is not "+gammalib::str(val)); \
        test_assert((col1.real(i) == double(val)), \
                    "Test real() method for row "+gammalib::str(i), \
                    gammalib::str(col1.real(i))+" is not "+gammalib::str(double(val))); \
        test_assert((col1.integer(i) == val), \
                    "Test integer() method for row "+gammalib::str(i), \
                    gammalib::str(col1.integer(i))+" is not "+gammalib::str(val)); \
        test_assert((col1.string(i) == gammalib::str(val)), \
                    "Test string() method for row "+gammalib::str(i), \
                    col1.string(i)+" is not "+gammalib::str(val)); \
    }
#define TEST_TABLE1_BOOL \
    for (int i = 0; i < nrows; ++i) { \
        col1(i) = (i % 2); \
    } \
    for (int i = 0; i < nrows; ++i) { \
        int         ival = (i % 2); \
        std::string sval = (i % 2) ? "T" : "F"; \
        test_assert((col1(i) == ival), \
                    "Test access operator row "+gammalib::str(i), \
                    gammalib::str(col1(i))+" is not "+gammalib::str(ival)); \
        test_assert((col1.real(i) == double(ival)), \
                    "Test real() method for row "+gammalib::str(i), \
                    gammalib::str(col1.real(i))+" is not "+gammalib::str(double(ival))); \
        test_assert((col1.integer(i) == int(ival)), \
                    "Test integer() method for row "+gammalib::str(i), \
                    gammalib::str(col1.integer(i))+" is not "+gammalib::str(int(ival))); \
        test_assert((col1.string(i) == gammalib::str(ival) || col1.string(i) == sval), \
                    "Test string() method for row "+gammalib::str(i), \
                    col1.string(i)+" is not "+gammalib::str(ival)+" or "+sval); \
    }
#define TEST_TABLE1_STRING \
    for (int i = 0; i < nrows; ++i) { \
        col1(i) = gammalib::str(double(i)*3.57+1.29); \
    } \
    for (int i = 0; i < nrows; ++i) { \
        double val = (double(i)*3.57+1.29); \
        test_value(gammalib::todouble(col1(i)), val, 1e-6, \
                   "Test access operator row "+gammalib::str(i)); \
        test_value(col1.real(i), val, 1e-6, \
                   "Test real() method for row "+gammalib::str(i)); \
        test_value(col1.integer(i), int(val), 1e-6, \
                   "Test integer() method for row "+gammalib::str(i)); \
        test_assert((col1.string(i) == gammalib::str(val)), \
                    "Test string() method for row "+gammalib::str(i), \
                    col1.string(i)+" is not "+gammalib::str(val)); \
    }
#define TEST_TABLE2 \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < nvec; ++j) { \
            col2(i,j) = double(i)*2.3 + double(j)*11.7 + 0.95; \
        } \
    } \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < nvec; ++j) { \
            double val = double(i)*2.3 + double(j)*11.7 + 0.95; \
            test_value(col2(i,j), val, 1e-6, \
                       "Test access operator row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_value(col2.real(i,j), val, 1e-6, \
                       "Test real() method for row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_value(col2.integer(i,j), int(val), 1e-6, \
                       "Test integer() method for row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_assert((col2.string(i,j) == gammalib::str(val)), \
                        "Test string() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        col2.string(i,j)+" is not "+gammalib::str(val)); \
        } \
    }
#define TEST_TABLE2_INT \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < nvec; ++j) { \
            col2(i,j) = i + j*11; \
        } \
    } \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < nvec; ++j) { \
            int val = i + j*11; \
            test_assert((col2(i,j) == val), \
                        "Test access operator row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col2(i,j))+" is not "+gammalib::str(val)); \
            test_assert((col2.real(i,j) == double(val)), \
                        "Test real() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col2.real(i,j))+" is not "+gammalib::str(double(val))); \
            test_assert((col2.integer(i,j) == val), \
                        "Test integer() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col2.integer(i,j))+" is not "+gammalib::str(val)); \
            test_assert((col2.string(i,j) == gammalib::str(val)), \
                        "Test string() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        col2.string(i,j)+" is not "+gammalib::str(val)); \
        } \
    }
#define TEST_TABLE2_BOOL \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < nvec; ++j) { \
            col2(i,j) = (i % 2) * (j % 2); \
        } \
    } \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < nvec; ++j) { \
            int         ival = ((i % 2) * (j % 2)); \
            std::string sval = ((i % 2) * (j % 2)) ? "T" : "F"; \
            test_assert((col2(i,j) == ival), \
                        "Test access operator row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col2(i,j))+" is not "+gammalib::str(ival)); \
            test_assert((col2.real(i,j) == double(ival)), \
                        "Test real() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col2.real(i,j))+" is not "+gammalib::str(double(ival))); \
            test_assert((col2.integer(i,j) == int(ival)), \
                        "Test integer() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col2.integer(i,j))+" is not "+gammalib::str(int(ival))); \
            test_assert((col2.string(i,j) == gammalib::str(ival) || col2.string(i,j) == sval), \
                        "Test string() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        col2.string(i,j)+" is not "+gammalib::str(ival)+" or "+sval); \
        } \
    }
#define TEST_TABLE2_STRING \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < nvec; ++j) { \
            col2(i,j) = gammalib::str(double(i)*2.3 + double(j)*11.7 + 0.95); \
        } \
    } \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < nvec; ++j) { \
            double val = double(i)*2.3 + double(j)*11.7 + 0.95; \
            test_value(gammalib::todouble(col2(i,j)), val, 1e-6, \
                       "Test access operator row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_value(col2.real(i,j), val, 1e-6, \
                       "Test real() method for row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_value(col2.integer(i,j), int(val), 1e-6, \
                       "Test integer() method for row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_assert((col2.string(i,j) == gammalib::str(val)), \
                        "Test string() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        col2.string(i,j)+" is not "+gammalib::str(val)); \
        } \
    }
#define TEST_VARTABLE \
    for (int i = 0; i < nrows; ++i) { \
        col3.elements(i,i+1); \
        for (int j = 0; j < i+1; ++j) { \
            col3(i,j) = double(i)*2.3 + double(j)*11.7 + 0.95; \
        } \
    } \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < i+1; ++j) { \
            double val = double(i)*2.3 + double(j)*11.7 + 0.95; \
            test_value(col3(i,j), val, 1e-6, \
                       "Test access operator row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_value(col3.real(i,j), val, 1e-6, \
                       "Test real() method for row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_value(col3.integer(i,j), int(val), 1e-6, \
                       "Test integer() method for row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_assert((col3.string(i,j) == gammalib::str(val)), \
                        "Test string() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        col3.string(i,j)+" is not "+gammalib::str(val)); \
        } \
    }
#define TEST_VARTABLE_INT \
    for (int i = 0; i < nrows; ++i) { \
        col3.elements(i,i+1); \
        for (int j = 0; j < i+1; ++j) { \
            col3(i,j) = i + j*11; \
        } \
    } \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < i+1; ++j) { \
            int val = i + j*11; \
            test_assert((col3(i,j) == val), \
                        "Test access operator row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col3(i,j))+" is not "+gammalib::str(val)); \
            test_assert((col3.real(i,j) == double(val)), \
                        "Test real() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col3.real(i,j))+" is not "+gammalib::str(double(val))); \
            test_assert((col3.integer(i,j) == val), \
                        "Test integer() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col3.integer(i,j))+" is not "+gammalib::str(val)); \
            test_assert((col3.string(i,j) == gammalib::str(val)), \
                        "Test string() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        col3.string(i,j)+" is not "+gammalib::str(val)); \
        } \
    }
#define TEST_VARTABLE_BOOL \
    for (int i = 0; i < nrows; ++i) { \
        col3.elements(i,i+1); \
        for (int j = 0; j < i+1; ++j) { \
            col3(i,j) = (i % 2) * (j % 2); \
        } \
    } \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < i+1; ++j) { \
            int         ival = ((i % 2) * (j % 2)); \
            std::string sval = ((i % 2) * (j % 2)) ? "T" : "F"; \
            test_assert((col3(i,j) == ival), \
                        "Test access operator row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col3(i,j))+" is not "+gammalib::str(ival)); \
            test_assert((col3.real(i,j) == double(ival)), \
                        "Test real() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col3.real(i,j))+" is not "+gammalib::str(double(ival))); \
            test_assert((col3.integer(i,j) == int(ival)), \
                        "Test integer() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        gammalib::str(col3.integer(i,j))+" is not "+gammalib::str(int(ival))); \
            test_assert((col3.string(i,j) == gammalib::str(ival) || col3.string(i,j) == sval), \
                        "Test string() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        col3.string(i,j)+" is not "+gammalib::str(ival)+" or "+sval); \
        } \
    }
#define TEST_VARTABLE_STRING \
    for (int i = 0; i < nrows; ++i) { \
        col3.elements(i,i+1); \
        for (int j = 0; j < i+1; ++j) { \
            col3(i,j) = gammalib::str(double(i)*2.3 + double(j)*11.7 + 0.95); \
        } \
    } \
    for (int i = 0; i < nrows; ++i) { \
        for (int j = 0; j < i+1; ++j) { \
            double val = double(i)*2.3 + double(j)*11.7 + 0.95; \
            test_value(gammalib::todouble(col3(i,j)), val, 1e-6, \
                       "Test access operator row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_value(col3.real(i,j), val, 1e-6, \
                       "Test real() method for row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_value(col3.integer(i,j), int(val), 1e-6, \
                       "Test integer() method for row "+gammalib::str(i)+" and element "+gammalib::str(j)); \
            test_assert((col3.string(i,j) == gammalib::str(val)), \
                        "Test string() method for row "+gammalib::str(i)+" and element "+gammalib::str(j), \
                        col3.string(i,j)+" is not "+gammalib::str(val)); \
        } \
    }
#define TEST_WRITE_TABLES \
    test_try("Write Table"); \
    try { \
        GFits fits; \
        fits.open(filename, true); \
        GFitsBinTable table(nrows); \
        table.append(col1); \
        table.insert(1, col2); \
        table.append(col3); \
        fits.append(table); \
        fits.save(); \
        fits.close(); \
        test_try_success(); \
    } \
    catch(std::exception &e) { \
        test_try_failure(e); \
    }


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGFits::set(void)
{
    // Set test name
    name("GFits");

    // Add tests
    append(static_cast<pfunction>(&TestGFits::test_create), "Test create");
    append(static_cast<pfunction>(&TestGFits::test_image_byte), "Test image byte");
    append(static_cast<pfunction>(&TestGFits::test_image_ushort), "Test image ushort");
    append(static_cast<pfunction>(&TestGFits::test_image_short), "Test image short");
    append(static_cast<pfunction>(&TestGFits::test_image_ulong), "Test image ulong");
    append(static_cast<pfunction>(&TestGFits::test_image_long), "Test image long");
    append(static_cast<pfunction>(&TestGFits::test_image_longlong), "Test image longlong");
    append(static_cast<pfunction>(&TestGFits::test_image_float), "Test image float");
    append(static_cast<pfunction>(&TestGFits::test_image_double), "Test image double");
    append(static_cast<pfunction>(&TestGFits::test_bintable_bit), "Test bintable bit");
    append(static_cast<pfunction>(&TestGFits::test_bintable_logical), "Test bintable logical");
    append(static_cast<pfunction>(&TestGFits::test_bintable_string), "Test bintable string");
    append(static_cast<pfunction>(&TestGFits::test_bintable_double), "Test bintable double");
    append(static_cast<pfunction>(&TestGFits::test_bintable_float), "Test bintable float");
    append(static_cast<pfunction>(&TestGFits::test_bintable_ushort), "Test bintable ushort");
    append(static_cast<pfunction>(&TestGFits::test_bintable_short), "Test bintable short");
    append(static_cast<pfunction>(&TestGFits::test_bintable_ulong), "Test bintable ulong");
    append(static_cast<pfunction>(&TestGFits::test_bintable_long), "Test bintable long");
    append(static_cast<pfunction>(&TestGFits::test_bintable_longlong), "Test bintable longlong");

    // Return
    return;
}


/***************************************************************************
 * @brief Test FITS file creation
 *
 * @todo Add checks that verify that the file content has been saved corretly
 ***************************************************************************/
void TestGFits::test_create(void)
{
    // Remove FITS file
    system("rm -rf test_empty.fits");
    system("rm -rf test_empty_image.fits");
    system("rm -rf test.fits");
    system("rm -rf test_create_bintable.fits");

    // Create empty FITS file
    test_try("Create empty FITS file");
    try {
        GFits fits("test_empty.fits", true);
        fits.save();
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
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
        test_try_failure(e);
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
    double total = 0.0;
    test_try("Re-open double precision image");
    try {
        GFits fits;
        fits.open("test.fits");
        GFitsHDU*         hdu   = fits.hdu(0);
        GFitsImageDouble* image = static_cast<GFitsImageDouble*>(fits.hdu(1));
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
        test_try_failure(e);
    }
    test_value(total, sum, 1e-10, "Test loading of double precision image");

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
        table.append(first);
        table.append(second);

        // Append table to FILE file
        fits.append(table);

        // Save FITS file
        fits.save(true);

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
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
        table.append(first);
        table.append(second);

        // Append table to FILE file
        fits.append(table);

        // Save FITS file
        fits.saveto("test_create_bintable.fits");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
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
    unsigned char* pixels = new unsigned char[16];
    for (int i = 0; i < 16; ++i) {
        pixels[i] = (unsigned char)(i);
    }
    
    // Test 1D pixel access
    GFitsImageByte image1(2, pixels);
    TEST_1D_ACCESS(2)

    // Test 2D pixel access
    GFitsImageByte image2(2, 2, pixels);
    TEST_2D_ACCESS(2,2)

    // Test 3D pixel access
    GFitsImageByte image3(2, 2, 2, pixels);
    TEST_3D_ACCESS(2,2,2)

    // Test 4D pixel access
    GFitsImageByte image4(2, 2, 2, 2, pixels);
    TEST_4D_ACCESS(2,2,2,2)

    // Test image I/O with 4D image
    int naxes[] = {2,2,2,2};
    GFitsImageByte image(4, naxes, pixels);

    // Save image
    GFits fits(filename, true);
    fits.append(image);
    fits.save();

    // Open FITS image
    GFits infile(filename);
    GFitsImage* ptr = infile.image(0);
    
    // Test 4D pixel access
    TEST_4D_ACCESS_IO(2,2,2,2)

    // Free pixels
    delete [] pixels;

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
    unsigned short* pixels = new unsigned short[16];
    for (int i = 0; i < 16; ++i) {
        pixels[i] = (unsigned short)(i);
    }

    // Test 1D pixel access
    GFitsImageUShort image1(2, pixels);
    TEST_1D_ACCESS(2)

    // Test 2D pixel access
    GFitsImageUShort image2(2, 2, pixels);
    TEST_2D_ACCESS(2,2)

    // Test 3D pixel access
    GFitsImageUShort image3(2, 2, 2, pixels);
    TEST_3D_ACCESS(2,2,2)

    // Test 4D pixel access
    GFitsImageUShort image4(2, 2, 2, 2, pixels);
    TEST_4D_ACCESS(2,2,2,2)

    // Test image I/O with 4D image
    int naxes[] = {2,2,2,2};
    GFitsImageUShort image(4, naxes, pixels);

    // Save image
    GFits fits(filename, true);
    fits.append(image);
    fits.save();

    // Open FITS image
    GFits infile(filename);
    GFitsImage* ptr = infile.image(0);
    
    // Test 4D pixel access
    TEST_4D_ACCESS_IO(2,2,2,2)

    // Free pixels
    delete [] pixels;
    
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
    short* pixels = new short[16];
    for (int i = 0; i < 16; ++i) {
        pixels[i] = (short)(i);
    }

    // Test 1D pixel access
    GFitsImageShort image1(2, pixels);
    TEST_1D_ACCESS(2)

    // Test 2D pixel access
    GFitsImageShort image2(2, 2, pixels);
    TEST_2D_ACCESS(2,2)

    // Test 3D pixel access
    GFitsImageShort image3(2, 2, 2, pixels);
    TEST_3D_ACCESS(2,2,2)

    // Test 4D pixel access
    GFitsImageShort image4(2, 2, 2, 2, pixels);
    TEST_4D_ACCESS(2,2,2,2)

    // Test image I/O with 4D image
    int naxes[] = {2,2,2,2};
    GFitsImageShort image(4, naxes, pixels);

    // Save image
    GFits fits(filename, true);
    fits.append(image);
    fits.save();

    // Open FITS image
    GFits infile(filename);
    GFitsImage* ptr = infile.image(0);
    
    // Test 4D pixel access
    TEST_4D_ACCESS_IO(2,2,2,2)

    // Free pixels
    delete [] pixels;

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
    unsigned long* pixels = new unsigned long[16];
    for (int i = 0; i < 16; ++i) {
        pixels[i] = (unsigned long)(i);
    }

    // Test 1D pixel access
    GFitsImageULong image1(2, pixels);
    TEST_1D_ACCESS(2)

    // Test 2D pixel access
    GFitsImageULong image2(2, 2, pixels);
    TEST_2D_ACCESS(2,2)

    // Test 3D pixel access
    GFitsImageULong image3(2, 2, 2, pixels);
    TEST_3D_ACCESS(2,2,2)

    // Test 4D pixel access
    GFitsImageULong image4(2, 2, 2, 2, pixels);
    TEST_4D_ACCESS(2,2,2,2)

    // Test image I/O with 4D image
    int naxes[] = {2,2,2,2};
    GFitsImageULong image(4, naxes, pixels);

    // Save image
    GFits fits(filename, true);
    fits.append(image);
    fits.save();

    // Open FITS image
    GFits infile(filename);
    GFitsImage* ptr = infile.image(0);
    
    // Test 4D pixel access
    TEST_4D_ACCESS_IO(2,2,2,2)

    // Free pixels
    delete [] pixels;

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
    long* pixels = new long[16];
    for (int i = 0; i < 16; ++i) {
        pixels[i] = (long)(i);
    }

    // Test 1D pixel access
    GFitsImageLong image1(2, pixels);
    TEST_1D_ACCESS(2)

    // Test 2D pixel access
    GFitsImageLong image2(2, 2, pixels);
    TEST_2D_ACCESS(2,2)

    // Test 3D pixel access
    GFitsImageLong image3(2, 2, 2, pixels);
    TEST_3D_ACCESS(2,2,2)

    // Test 4D pixel access
    GFitsImageLong image4(2, 2, 2, 2, pixels);
    TEST_4D_ACCESS(2,2,2,2)

    // Test image I/O with 4D image
    int naxes[] = {2,2,2,2};
    GFitsImageLong image(4, naxes, pixels);

    // Save image
    GFits fits(filename, true);
    fits.append(image);
    fits.save();

    // Open FITS image
    GFits infile(filename);
    GFitsImage* ptr = infile.image(0);
    
    // Test 4D pixel access
    TEST_4D_ACCESS_IO(2,2,2,2)

    // Free pixels
    delete [] pixels;

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
    long long* pixels = new long long[16];
    for (int i = 0; i < 16; ++i) {
        pixels[i] = (long long)(i);
    }

    // Test 1D pixel access
    GFitsImageLongLong image1(2, pixels);
    TEST_1D_ACCESS(2)

    // Test 2D pixel access
    GFitsImageLongLong image2(2, 2, pixels);
    TEST_2D_ACCESS(2,2)

    // Test 3D pixel access
    GFitsImageLongLong image3(2, 2, 2, pixels);
    TEST_3D_ACCESS(2,2,2)

    // Test 4D pixel access
    GFitsImageLongLong image4(2, 2, 2, 2, pixels);
    TEST_4D_ACCESS(2,2,2,2)

    // Test image I/O with 4D image
    int naxes[] = {2,2,2,2};
    GFitsImageLongLong image(4, naxes, pixels);

    // Save image
    GFits fits(filename, true);
    fits.append(image);
    fits.save();

    // Open FITS image
    GFits infile(filename);
    GFitsImage* ptr = infile.image(0);
    
    // Test 4D pixel access
    TEST_4D_ACCESS_IO(2,2,2,2)

    // Free pixels
    delete [] pixels;

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
    float* pixels = new float[16];
    for (int i = 0; i < 16; ++i) {
        pixels[i] = (float)(i);
    }

    // Test 1D pixel access
    GFitsImageFloat image1(2, pixels);
    TEST_1D_ACCESS(2)

    // Test 2D pixel access
    GFitsImageFloat image2(2, 2, pixels);
    TEST_2D_ACCESS(2,2)

    // Test 3D pixel access
    GFitsImageFloat image3(2, 2, 2, pixels);
    TEST_3D_ACCESS(2,2,2)

    // Test 4D pixel access
    GFitsImageFloat image4(2, 2, 2, 2, pixels);
    TEST_4D_ACCESS(2,2,2,2)

    // Test image I/O with 4D image
    int naxes[] = {2,2,2,2};
    GFitsImageFloat image(4, naxes, pixels);

    // Save image
    GFits fits(filename, true);
    fits.append(image);
    fits.save();

    // Open FITS image
    GFits infile(filename);
    GFitsImage* ptr = infile.image(0);
    
    // Test 4D pixel access
    TEST_4D_ACCESS_IO(2,2,2,2)

    // Free pixels
    delete [] pixels;

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
    double* pixels = new double[16];
    for (int i = 0; i < 16; ++i) {
        pixels[i] = (double)(i);
    }

    // Test 1D pixel access
    GFitsImageDouble image1(2, pixels);
    TEST_1D_ACCESS(2)

    // Test 2D pixel access
    GFitsImageDouble image2(2, 2, pixels);
    TEST_2D_ACCESS(2,2)

    // Test 3D pixel access
    GFitsImageDouble image3(2, 2, 2, pixels);
    TEST_3D_ACCESS(2,2,2)

    // Test 4D pixel access
    GFitsImageDouble image4(2, 2, 2, 2, pixels);
    TEST_4D_ACCESS(2,2,2,2)

    // Test image I/O with 4D image
    int naxes[] = {2,2,2,2};
    GFitsImageDouble image(4, naxes, pixels);

    // Save image
    GFits fits(filename, true);
    fits.append(image);
    fits.save();

    // Open FITS image
    GFits infile(filename);
    GFitsImage* ptr = infile.image(0);
    
    // Test 4D pixel access
    TEST_4D_ACCESS_IO(2,2,2,2)

    // Free pixels
    delete [] pixels;

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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableDoubleCol col1("DOUBLE", nrows);
    TEST_TABLE1;

    // Test multiple column table
    GFitsTableDoubleCol col2("DOUBLE10", nrows, nvec);
    TEST_TABLE2;

    // Test variable-length column table
    GFitsTableDoubleCol col3("DOUBLEVAR", nrows, -1);
    TEST_VARTABLE;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableDoubleCol&>(*(*fits.table(1))["DOUBLE"]);
        col2 = static_cast<GFitsTableDoubleCol&>(*(*fits.table(1))["DOUBLE10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1;

    // Test multiple column table
    TEST_TABLE2;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableFloatCol col1("FLOAT", nrows);
    TEST_TABLE1;

    // Test multiple column table
    GFitsTableFloatCol col2("FLOAT10", nrows, nvec);
    TEST_TABLE2;

    // Test variable-length column table
    GFitsTableFloatCol col3("FLOATVAR", nrows, -1);
    TEST_VARTABLE;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableFloatCol&>(*(*fits.table(1))["FLOAT"]);
        col2 = static_cast<GFitsTableFloatCol&>(*(*fits.table(1))["FLOAT10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1;

    // Test multiple column table
    TEST_TABLE2;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableShortCol col1("SHORT", nrows);
    TEST_TABLE1_INT;

    // Test multiple column table
    GFitsTableShortCol col2("SHORT10", nrows, nvec);
    TEST_TABLE2_INT;

    // Test variable-length column table
    GFitsTableShortCol col3("SHORTVAR", nrows, -1);
    TEST_VARTABLE_INT;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableShortCol&>(*(*fits.table(1))["SHORT"]);
        col2 = static_cast<GFitsTableShortCol&>(*(*fits.table(1))["SHORT10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1_INT;

    // Test multiple column table
    TEST_TABLE2_INT;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableUShortCol col1("USHORT", nrows);
    TEST_TABLE1_INT;

    // Test multiple column table
    GFitsTableUShortCol col2("USHORT10", nrows, nvec);
    TEST_TABLE2_INT;

    // Test variable-length column table
    GFitsTableUShortCol col3("USHORTVAR", nrows, -1);
    TEST_VARTABLE_INT;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableUShortCol&>(*(*fits.table(1))["USHORT"]);
        col2 = static_cast<GFitsTableUShortCol&>(*(*fits.table(1))["USHORT10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1_INT;

    // Test multiple column table
    TEST_TABLE2_INT;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableLongCol col1("LONG", nrows);
    TEST_TABLE1_INT;

    // Test multiple column table
    GFitsTableLongCol col2("LONG10", nrows, nvec);
    TEST_TABLE2_INT;

    // Test variable-length column table
    GFitsTableLongCol col3("LONGVAR", nrows, -1);
    TEST_VARTABLE_INT;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableLongCol&>(*(*fits.table(1))["LONG"]);
        col2 = static_cast<GFitsTableLongCol&>(*(*fits.table(1))["LONG10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1_INT;

    // Test multiple column table
    TEST_TABLE2_INT;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableLongLongCol col1("LONGLONG", nrows);
    TEST_TABLE1_INT;

    // Test multiple column table
    GFitsTableLongLongCol col2("LONGLONG10", nrows, nvec);
    TEST_TABLE2_INT;

    // Test variable-length column table
    GFitsTableLongLongCol col3("LONGLONGVAR", nrows, -1);
    TEST_VARTABLE_INT;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableLongLongCol&>(*(*fits.table(1))["LONGLONG"]);
        col2 = static_cast<GFitsTableLongLongCol&>(*(*fits.table(1))["LONGLONG10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1_INT;

    // Test multiple column table
    TEST_TABLE2_INT;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableULongCol col1("ULONG", nrows);
    TEST_TABLE1_INT;

    // Test multiple column table
    GFitsTableULongCol col2("ULONG10", nrows, nvec);
    TEST_TABLE2_INT;

    // Test variable-length column table
    GFitsTableULongCol col3("ULONGVAR", nrows, -1);
    TEST_VARTABLE_INT;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableULongCol&>(*(*fits.table(1))["ULONG"]);
        col2 = static_cast<GFitsTableULongCol&>(*(*fits.table(1))["ULONG10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1_INT;

    // Test multiple column table
    TEST_TABLE2_INT;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableStringCol col1("STRING", nrows, 20);
    TEST_TABLE1_STRING;

    // Test multiple column table
    GFitsTableStringCol col2("STRING10", nrows, 20, nvec);
    TEST_TABLE2_STRING;

    // Test variable-length column table
    GFitsTableStringCol col3("STRINGVAR", nrows, 1);
    //GFitsTableStringCol col3("STRINGVAR", nrows, -1);
    //TEST_VARTABLE_STRING;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableStringCol&>(*(*fits.table(1))["STRING"]);
        col2 = static_cast<GFitsTableStringCol&>(*(*fits.table(1))["STRING10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1_STRING;

    // Test multiple column table
    TEST_TABLE2_STRING;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableBoolCol col1("LOGICAL", nrows);
    TEST_TABLE1_BOOL;

    // Test multiple column table
    GFitsTableBoolCol col2("LOGICAL10", nrows, nvec);
    TEST_TABLE2_BOOL;

    // Test variable-length column table
    GFitsTableBoolCol col3("LOGICALVAR", nrows, -1);
    TEST_VARTABLE_BOOL;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = static_cast<GFitsTableBoolCol&>(*(*fits.table(1))["LOGICAL"]);
        col2 = static_cast<GFitsTableBoolCol&>(*(*fits.table(1))["LOGICAL10"]);
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }
    
    // Test single column table
    TEST_TABLE1_BOOL;

    // Test multiple column table
    TEST_TABLE2_BOOL;

    // Return
    return;
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
    system(cmd.c_str());

    // Set number of rows and vector columns
    int nrows = 3;
    int nvec  = 3;
    
    // Test single column table
    GFitsTableBitCol col1("BIT", nrows);
    TEST_TABLE1_BOOL;

    // Test multiple column table
    GFitsTableBitCol col2("BIT10", nrows, nvec);
    TEST_TABLE2_BOOL;

    // Test variable-length column table
    GFitsTableBitCol col3("BITVAR", nrows, -1);
    //TEST_VARTABLE_BOOL;

    // Write tables
    TEST_WRITE_TABLES;

    // Read tables back
    test_try("Read Tables");
    try {
        GFits fits(filename);
        col1 = *(static_cast<GFitsTableBitCol*>((*fits.table(1))["BIT"]));
        col2 = *(static_cast<GFitsTableBitCol*>((*fits.table(1))["BIT10"]));
        fits.close();
        test_try_success();
    }
    catch(std::exception &e) {
        test_try_failure(e);
    }

    // Test single column table
    TEST_TABLE1_BOOL;

    // Test multiple column table
    TEST_TABLE2_BOOL;

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("FITS module");

    // Initially assume that we pass all tests
    bool success = true;

    // Create a test suite
    TestGFits test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GFits.xml");

    // Return success status
    return (success ? 0 : 1);
}
