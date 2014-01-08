/***************************************************************************
 *            test_GApplication.cpp - test GApplication classes            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @brief Testing of application classes
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "test_GApplication.hpp"

/* __ Globals ____________________________________________________________ */


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGApplication::set(void)
{
    // Test name
    name("GApplication");

    // Append tests
    append(static_cast<pfunction>(&TestGApplication::test_constructor), "Test GLog constructor");
    append(static_cast<pfunction>(&TestGApplication::test_stream_logger), "Test stream logger");
    append(static_cast<pfunction>(&TestGApplication::test_C_logger), "Test C logger");
    append(static_cast<pfunction>(&TestGApplication::test_GPar), "Test GPar class");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGApplication* TestGApplication::clone(void) const
{
    // Clone test suite
    return new TestGApplication(*this);
}


/***********************************************************************//**
 * @brief Test GLog constructor
 **************************************************************************/
void TestGApplication::test_constructor(void)
{
    // Void constructor
    GLog log1;

    // Copy constructor
    GLog log2 = log1;

    // Return
    return; 
}


/***********************************************************************//**
 * @brief Test stream logger
 **************************************************************************/
void TestGApplication::test_stream_logger(void)
{
    // Test
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test C logger
 **************************************************************************/
void TestGApplication::test_C_logger(void)
{
    // Test
    GLog log;
    log.date(true);
    log.name("Test");
    log.open("test_GApplication.log");
    log("%s", "8. This is a C++ string");
    log("%s %d", "9. This is an integer:", int(1));
    log("%s %f", "10. This is a single precision:", float(1.23456789));
    log("%s %f", "11. This is a double precision:", double(1.23456789));

    // Return
    return; 
}


/***********************************************************************//**
 * @brief Test GPar class
 **************************************************************************/
void TestGApplication::test_GPar(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GPar par;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test copy constructor
    test_try("Copy constructor");
    try {
        GPar par;
        GPar par2(par);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test parameter constructor
    test_try("Parameter constructor");
    try {
        GPar par("name", "r", "a", "1.0", "0.0", "2.0", "Parameter name");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test integer parameter exceptions
    GPar par;
    par = GPar("name", "i", "a", "INDEF", "0.0", "2.0", "Parameter name");
    test_assert(par.is_undefined(), "Check integer parameter INDEF.",
                par.value()+" found instead of undefined value.");
    par = GPar("name", "i", "a", "NONE", "0.0", "2.0", "Parameter name");
    test_assert(par.is_undefined(), "Check integer parameter NONE.",
                par.value()+" found instead of undefined value.");
    par = GPar("name", "i", "a", "UNDEF", "0.0", "2.0", "Parameter name");
    test_assert(par.is_undefined(), "Check integer parameter UNDEF.",
                par.value()+" found instead of undefined value.");
    par = GPar("name", "i", "a", "UNDEFINED", "0.0", "2.0", "Parameter name");
    test_assert(par.is_undefined(), "Check integer parameter UNDEFINED.",
                par.value()+" found instead of undefined value.");

    // Test floating point parameter exceptions
    par = GPar("name", "r", "a", "INDEF", "0.0", "2.0", "Parameter name");
    test_assert(par.is_undefined(), "Check floating point parameter INDEF.",
                par.value()+" found instead of undefined value.");
    par = GPar("name", "r", "a", "NONE", "0.0", "2.0", "Parameter name");
    test_assert(par.is_undefined(), "Check floating point parameter NONE.",
                par.value()+" found instead of undefined value.");
    par = GPar("name", "r", "a", "UNDEF", "0.0", "2.0", "Parameter name");
    test_assert(par.is_undefined(), "Check floating point parameter UNDEF.",
                par.value()+" found instead of undefined value.");
    par = GPar("name", "r", "a", "UNDEFINED", "0.0", "2.0", "Parameter name");
    test_assert(par.is_undefined(), "Check floating point parameter UNDEFINED.",
                par.value()+" found instead of undefined value.");
    par = GPar("name", "r", "a", "INF", "0.0", "2.0", "Parameter name");
    test_assert(par.is_notanumber(), "Check floating point parameter INF.",
                par.value()+" found instead of infinite value.");
    par = GPar("name", "r", "a", "INFINITY", "0.0", "2.0", "Parameter name");
    test_assert(par.is_notanumber(), "Check floating point parameter INFINITY.",
                par.value()+" found instead of infinite value.");
    par = GPar("name", "r", "a", "NAN", "0.0", "2.0", "Parameter name");
    test_assert(par.is_notanumber(), "Check floating point parameter NAN.",
                par.value()+" found instead of not a number.");

    // Return
    return; 
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("GApplication");

    // Initially assume that we pass all tests
    bool success = true;

    // Create a test suite
    TestGApplication test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GApplication.xml");

    // Return success status
    return (success ? 0 : 1);
}
