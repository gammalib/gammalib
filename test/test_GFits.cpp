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
 *                                 Test: Open                              *
 ***************************************************************************/
void test_open(void)
{
    std::cout << "Test GFits: Open FITS file: " << endl;
    try {
        GFits fits;
        fits.open("srcid.fits");
        
        // Access data via conversion functions
        cout << fits.hdu("LAT_POINT_SOURCE_CATALOG")->column("ID_3EG_PROB_1")->string(0) << endl;
        cout << fits.hdu("LAT_POINT_SOURCE_CATALOG")->column("ID_3EG_PROB_1")->real(0) << endl;
        cout << fits.hdu("LAT_POINT_SOURCE_CATALOG")->column("ID_3EG_PROB_1")->integer(0) << endl;

        // Access data via pointer
        double* ptr = fits.hdu("LAT_POINT_SOURCE_CATALOG")->column("ID_3EG_PROB_1")->ptr_double();
        cout << ptr[0] << endl;
        
        GFitsTableDblCol column = *((GFitsTableDblCol*)fits.hdu("LAT_POINT_SOURCE_CATALOG")->column("ID_3EG_PROB_1"));
        cout << column.real(0) << endl;

    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to open FITS file." << std::endl;
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

    // Return
    return 0;
}
