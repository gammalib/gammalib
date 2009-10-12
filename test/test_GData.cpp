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
 *  Test: LAT unbinned observation handling                                *
 ***************************************************************************/
void test_lat_unbinned(void)
{
	// Write header
    std::cout << "Test LAT unbinned data handling: ";

	// Declare observations
	GData           data;
	GLATObservation obs;

	// Load unbinned LAT observation
    try {
        obs.load_unbinned("data/lat/ft1.fits.gz", "data/lat/ft2.fits.gz", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load LAT observation." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
//std::cout << obs << std::endl;
	
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
//            std::cout << *event << std::endl;
//            std::cout << event->test() << std::endl;
            num++;
        }
        if (num != 4038) {
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

    // Loop over all events using iterator
    try {
        int num = 0;
        GLATEventList *ptr = (GLATEventList*)obs.events();
        for (GLATEventList::iterator event = ptr->begin(); event != ptr->end(); ++event) {
//            std::cout << *event << std::endl;
            num++;
        }
        if (num != 2019) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in GLATEventList::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GLATEventList." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    
    // Plot final test success
    std::cout << ". ok." << std::endl;

//GLATEventList *ptr = (GLATEventList*)obs.events();
//std::cout << obs << std::endl;
//std::cout << *ptr << std::endl;

    // Exit test
    return;
 
}


/***************************************************************************
 *  Test: LAT binned observation handling                                *
 ***************************************************************************/
void test_lat_binned(void)
{
	// Write header
    std::cout << "Test LAT binned data handling: ";

	// Declare observations
	GData           data;
	GLATObservation obs;

    // Load LAT binned observation
    try {
        obs.load_binned("data/lat/cntmap.fits.gz", "", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to load binned LAT observation." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Reload LAT binned observation from srcmap
    try {
        obs.load_binned("data/lat/srcmap.fits.gz", "", "");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to reload binned LAT observation." 
                  << std::endl;
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
        int sum = 0;
        for (GData::iterator event = data.begin(); event != data.end(); ++event) {
            num++;
            sum += event->counts();
        }
        if (sum != 2718 || num != 400000) {
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

    // Loop over all events using iterator
    try {
        int num = 0;
        int sum = 0;
        GLATEventCube *ptr = (GLATEventCube*)obs.events();
        for (GLATEventCube::iterator event = ptr->begin(); event != ptr->end(); ++event) {
//            std::cout << *((GLATEventBin*)&(*event));
            num++;
            sum += event->counts();
        }
        if (sum != 1359 || num != 200000) {
            std::cout << std::endl << 
                      "TEST ERROR: Wrong number of iterations in GLATEventCube::iterator."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to iterate GLATEventCube." << std::endl;
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
    test_lat_unbinned();
    test_lat_binned();

    // Return
    return 0;
}
