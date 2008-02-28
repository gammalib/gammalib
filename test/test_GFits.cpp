/***************************************************************************
 *                 test_GMatrix.cpp  -  test FITS classes                  *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "test_GFits.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *                              Test: Open                                 *
 ***************************************************************************/
int fequal(double val, double ref)
{
    return (fabs(val-ref) < 1.0e-5) ? 1 : 0;
}


/***************************************************************************
 *                              Test: Open                                 *
 ***************************************************************************/
void test_open(void)
{
    std::cout << "Test GFits: Open FITS file: ";
    try {
        GFits fits;
        fits.open("test_GFits.fits");
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to open FITS file." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;
}


/***************************************************************************
 *                             Test: Columns                               *
 ***************************************************************************/
void test_columns(void)
{
    GFits fits;
    try {
        // Open FITS file
        fits.open("test_GFits.fits");
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to open test file." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }

    // Test single floating point column
    std::cout << "Test GFits: Floating point column handling: ";
    try {
        GFitsTableFltCol flt = *((GFitsTableFltCol*)fits.hdu("BinTable")->column("TFLOAT"));
        if (!fequal(flt.real(0),1.0) > 1e-10 || !fequal(flt.real(1),-2.1)) {
          cout << endl << "TEST ERROR: Bad single TFLOAT values read (" << flt.real(0) << ", " << flt.real(1) << ")." << endl;
          throw;
        }
        else
          std::cout << ".";
        float nullval = -99.9;
        flt.set_nullval(&nullval);
        if (!fequal(flt.real(2),-99.9)) {
          cout << endl << "TEST ERROR: Bad single TFLOAT NULL value read (" << flt.real(2) << ")." << endl;
          throw;
        }
        else
          std::cout << ".";
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to handle single TFLOAT column." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "***********************" << std::endl;
    std::cout << "* GFIts class testing *" << std::endl;
    std::cout << "***********************" << std::endl;

    // Execute the tests
    test_open();
    test_columns();

    // Return
    return 0;
}
