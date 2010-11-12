/***************************************************************************
 *           test_GApplication.cpp  -  test GApplication classes           *
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
 * @file test_GApplication.cpp
 * @brief Testing of application classes.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GammaLib.hpp"
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 * Test: GLog                                                              *
 ***************************************************************************/
void test_GLog(void)
{
    // Dump header
    std::cout << "Test GLog: ";

    // Test constructor
    try {
        GLog log1;
        GLog log2 = log1;
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct GLog object."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test stream logger
    try {
        GLog log;
        log.date(true);
        log.name("Test");
        log.open("test_GApplication.log", true);
        log << "1. This is a C++ string" << std::endl;
        log << "2. This is an integer: " << int(1) << std::endl;
        log << "3. This is a single precision: " << float(1.23456789) << std::endl;
        log << "4. This is a double precision: " << double(1.23456789) << std::endl;
        log << "5. This is a character: " << 'a' << std::endl;
        log << "6. This is a C string: " << "a" << std::endl;
        log << "7. This is a Boolean: " << true << std::endl;
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to stream into GLog object."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test C logger
    try {
        GLog log;
        log.date(true);
        log.name("Test");
        log.open("test_GApplication.log");
        log("%s", "8. This is a C++ string");
        log("%s %d", "9. This is an integer:", int(1));
        log("%s %f", "10. This is a single precision:", float(1.23456789));
        log("%s %f", "11. This is a double precision:", double(1.23456789));
    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Unable to write into GLog object."
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
    std::cout << "********************************" << std::endl;
    std::cout << "* GApplication classes testing *" << std::endl;
    std::cout << "********************************" << std::endl;

    // Execute tests
    test_GLog();

    // Return
    return 0;
}
