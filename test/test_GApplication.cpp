/***************************************************************************
 *           test_GApplication.cpp  -  test GApplication classes           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jurgen Knodlseder                                *
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
 * @file test_GApplication.cpp
 * @brief Testing of application classes.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GammaLib.hpp"
#include "test_GApplication.hpp"
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

/* __ Globals ____________________________________________________________ */

void TestGApplication::set(void){
    // Test name
    name("GApplication");

    //add tests
    add_test(static_cast<pfunction>(&TestGApplication::test_constructor),"GLog constructor");
    add_test(static_cast<pfunction>(&TestGApplication::test_stream_logger),"stream logger");
    add_test(static_cast<pfunction>(&TestGApplication::test_C_logger),"C logger");

    return;
}


// Test constructor
void TestGApplication::test_constructor(void)
{
    GLog log1;
    GLog log2 = log1;

    return; 
}

// Test stream logger
void TestGApplication::test_stream_logger(void)
{
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

    return;
}

// Test C logger
void TestGApplication::test_C_logger(void)
{
    GLog log;
    log.date(true);
    log.name("Test");
    log.open("test_GApplication.log");
    log("%s", "8. This is a C++ string");
    log("%s %d", "9. This is an integer:", int(1));
    log("%s %f", "10. This is a single precision:", float(1.23456789));
    log("%s %f", "11. This is a double precision:", double(1.23456789));

    return; 
}

/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuites("GApplication");

    bool was_successful=true;

    //Create a test suite
    TestGApplication test;

    //Append to the container
    testsuites.append(test);

    //Run
    was_successful=testsuites.run();

    //save xml report
    testsuites.save("reports/GApplication.xml");

    // Return
    return was_successful ? 0:1;
}
