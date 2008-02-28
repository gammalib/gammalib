/***************************************************************************
 *              test_GResponse.cpp  -  test Response classes               *
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
#include "test_GResponse.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *                              Test: Open                                 *
 ***************************************************************************/
void test_open(void)
{
    std::cout << "Test GResponse: Open PSF FITS file: ";
    try {
        // Open FITS file
        GFits fits;
        fits.open("irf/lat/aeff_Pass5_v0_front.fits");

        // Get pointer towards effective area HDU
        GFitsHDU* aeff = fits.hdu("EFFECTIVE AREA");
        if (aeff == NULL) {
            std::cout << std::endl << "TEST ERROR: Unable to find HDU <EFFECTIVE AREA>." << std::endl;
            throw;
        }

        // Get columns
        GFitsTableFltCol energ_lo  = *((GFitsTableFltCol*)aeff->column("ENERG_LO"));
        GFitsTableFltCol energ_hi  = *((GFitsTableFltCol*)aeff->column("ENERG_HI"));
        GFitsTableFltCol ctheta_lo = *((GFitsTableFltCol*)aeff->column("CTHETA_LO"));
        GFitsTableFltCol ctheta_hi = *((GFitsTableFltCol*)aeff->column("CTHETA_HI"));
        GFitsTableFltCol effarea   = *((GFitsTableFltCol*)aeff->column("EFFAREA"));

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
    std::cout << "***************************" << std::endl;
    std::cout << "* GResponse class testing *" << std::endl;
    std::cout << "**************************" << std::endl;

    // Execute the tests
    test_open();

    // Return
    return 0;
}
