/***************************************************************************
 *                   test_GData.cpp  -  test GData class                   *
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
#include <stdlib.h>
#include "test_GData.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *                          Test: test_lat_data                            *
 ***************************************************************************/
void test_lat_data(void)
{
	// Write header
    std::cout << "Test LAT data handling: ";

	// Declare observations
	GData           data;
	GLATObservation obs;

	// Load LAT observation
    try {
		obs = GLATObservation("data/FT1_253582800.fits.gz", "data/FT2_253582800.fits.gz");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load LAT observation." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;
	std::cout << obs << std::endl;
	
	// Add observation to data
	data.add(&obs);
	std::cout << data << std::endl;
}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "***********************" << std::endl;
    std::cout << "* GData class testing *" << std::endl;
    std::cout << "***********************" << std::endl;

    // Execute the tests
    test_lat_data();

    // Return
    return 0;
}
