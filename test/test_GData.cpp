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
    std::cout << ".";
	
	// Add observation (twice) to data
    try {
        data.add(obs);
        data.add(obs);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to add LAT observation to data." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    
    // Loop over all events in GData using iterators
    try {
        int num = 0;
        for (GData::iterator event = data.begin(); event != data.end(); ++event) {
//          std::cout << *event << std::endl;
            num++;
        }
        if (num != 552724) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in GData::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GData." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Loop over all events in GEvents using iterator
    try {
        int num = 0;
        GLATEvents *ptr = (GLATEvents*)obs.events();
        for (GLATEvents::iterator event = ptr->begin(); event != ptr->end(); ++event) {
//          std::cout << *event << std::endl;
            num++;
        }
        if (num != 276362) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in GLATEvents::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GLATEvents." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
    std::cout << ". ok." << std::endl;

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
    std::cout << "***********************" << std::endl;
    std::cout << "* GData class testing *" << std::endl;
    std::cout << "***********************" << std::endl;

    // Execute the tests
    test_lat_data();

    // Return
    return 0;
}
