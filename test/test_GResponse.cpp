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
#include <iostream>                           // cout, cerr

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *                          Test: LAT Response                             *
 ***************************************************************************/
void test_lat_response(void)
{
    cout << "Test GLATResponse: Open PSF FITS file: ";
    try {
        // Get HANDOFF Response
        GLATResponse rsp;
        rsp.set_caldb("irf/lat");
        rsp.load("Pass5_v0", "front");
/*
        // Open FITS file
        GFits fits;
        fits.open("irf/lat/aeff_.fits");

        // Get pointer towards effective area HDU
        GFitsHDU* aeff = fits.hdu("EFFECTIVE AREA");
        if (aeff == NULL) {
            cout << endl << "TEST ERROR: Unable to find HDU <EFFECTIVE AREA>." << endl;
            throw;
        }

        // Get columns
        GFitsTableFltCol energ_lo  = *((GFitsTableFltCol*)aeff->column("ENERG_LO"));
        GFitsTableFltCol energ_hi  = *((GFitsTableFltCol*)aeff->column("ENERG_HI"));
        GFitsTableFltCol ctheta_lo = *((GFitsTableFltCol*)aeff->column("CTHETA_LO"));
        GFitsTableFltCol ctheta_hi = *((GFitsTableFltCol*)aeff->column("CTHETA_HI"));
        GFitsTableFltCol effarea   = *((GFitsTableFltCol*)aeff->column("EFFAREA"));
*/
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to open FITS file." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;
}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Dump header
    cout << std::endl;
    cout << "***************************" << endl;
    cout << "* GResponse class testing *" << endl;
    cout << "**************************" << endl;

    // Execute the tests
    test_lat_response();

    // Return
    return 0;
}
