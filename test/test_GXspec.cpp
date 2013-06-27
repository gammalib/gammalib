/***************************************************************************
 *                  test_GXspec.hpp - Test Xspec module                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file test_GXspec.cpp
 * @brief Implementation of unit tests for Xspec module
 * @author Juergen Knoedlseder 
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "test_GXspec.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGXspec::set(void)
{
    // Test name
    name("GXspec");

    // Append tests
    append(static_cast<pfunction>(&TestGXspec::test_GPha), "Test GPha class");
    append(static_cast<pfunction>(&TestGXspec::test_GArf), "Test GArf class");

    // Return
    return; 
}


/***********************************************************************//**
 * @brief Test GPha class
 **************************************************************************/
void TestGXspec::test_GPha(void)
{
    // Test void constructor
    test_try("GPha void constructor");
    try {
        GPha pha;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test energy boundary constructor
    test_try("GPha energy boundary constructor");
    try {
        GEbounds ebds(10, GEnergy(0.1, "TeV"), GEnergy(10.0, "TeV"));
        GPha     pha(ebds);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test filling and accessing
    GEbounds ebds(9, GEnergy(1.0, "TeV"), GEnergy(10.0, "TeV"), false);
    GPha     pha(ebds);
    for (int i = 0; i < 12; i += 2) {
        GEnergy energy(double(i), "TeV");
        pha.fill(energy);
    }
    test_value(pha.counts(),    4.0);
    test_value(pha.underflow(), 1.0);
    test_value(pha.overflow(),  1.0);
    test_value(pha.outflow(),   0.0);
    for (int i = 0; i < 9; i += 2) {
        test_value(pha.at(i), 0.0);
        test_value(pha[i], 0.0);
    }
    for (int i = 1; i < 9; i += 2) {
        test_value(pha.at(i), 1.0);
        test_value(pha[i], 1.0);
    }
    pha[0] = 5.0;
    pha[1] = 3.7;
    test_value(pha[0], 5.0);
    test_value(pha[1], 3.7);

    // Test saving and loading
    pha.save("pha.fits", true);
    pha.load("pha.fits");
    test_value(pha[0], 5.0, 1.0e-6);
    test_value(pha[1], 3.7, 1.0e-6);
    test_value(pha.counts(),   11.7, 1.0e-6);
    test_value(pha.underflow(), 0.0, 1.0e-6);
    test_value(pha.overflow(),  0.0, 1.0e-6);
    test_value(pha.outflow(),   0.0, 1.0e-6);
    for (int i = 2; i < 9; i += 2) {
        test_value(pha.at(i), 0.0);
        test_value(pha[i], 0.0);
    }
    for (int i = 3; i < 9; i += 2) {
        test_value(pha.at(i), 1.0);
        test_value(pha[i], 1.0);
    }

    // Test constructing
    GPha pha2("pha.fits");
    test_value(pha2[0], 5.0, 1.0e-6);
    test_value(pha2[1], 3.7, 1.0e-6);
    test_value(pha2.counts(),   11.7, 1.0e-6);
    test_value(pha2.underflow(), 0.0, 1.0e-6);
    test_value(pha2.overflow(),  0.0, 1.0e-6);
    test_value(pha2.outflow(),   0.0, 1.0e-6);
    for (int i = 2; i < 9; i += 2) {
        test_value(pha2.at(i), 0.0);
        test_value(pha2[i], 0.0);
    }
    for (int i = 3; i < 9; i += 2) {
        test_value(pha2.at(i), 1.0);
        test_value(pha2[i], 1.0);
    }
    //std::cout << pha << std::endl;
    //std::cout << pha2 << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GArf class
 **************************************************************************/
void TestGXspec::test_GArf(void)
{
    // Test void constructor
    test_try("GArf void constructor");
    try {
        GArf arf;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test energy boundary constructor
    test_try("GArf energy boundary constructor");
    try {
        GEbounds ebds(10, GEnergy(0.1, "TeV"), GEnergy(10.0, "TeV"));
        GArf     arf(ebds);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test filling and accessing
    GEbounds ebds(9, GEnergy(1.0, "TeV"), GEnergy(10.0, "TeV"));
    GArf     arf(ebds);
    for (int i = 0; i < 9; i += 2) {
        arf[i] = 1.0;
    }
    for (int i = 0; i < 9; i += 2) {
        test_value(arf.at(i), 1.0);
        test_value(arf[i], 1.0);
    }
    for (int i = 1; i < 9; i += 2) {
        test_value(arf.at(i), 0.0);
        test_value(arf[i], 0.0);
    }
    arf[0] = 5.0;
    arf[1] = 3.7;
    test_value(arf[0], 5.0);
    test_value(arf[1], 3.7);

    // Test saving and loading
    arf.save("arf.fits", true);
    arf.load("arf.fits");
    test_value(arf[0], 5.0, 1.0e-6);
    test_value(arf[1], 3.7, 1.0e-6);
    for (int i = 2; i < 9; i += 2) {
        test_value(arf.at(i), 1.0);
        test_value(arf[i], 1.0);
    }
    for (int i = 3; i < 9; i += 2) {
        test_value(arf.at(i), 0.0);
        test_value(arf[i], 0.0);
    }

    // Test constructing
    GArf arf2("arf.fits");
    test_value(arf2[0], 5.0, 1.0e-6);
    test_value(arf2[1], 3.7, 1.0e-6);
    for (int i = 2; i < 9; i += 2) {
        test_value(arf2.at(i), 1.0);
        test_value(arf2[i], 1.0);
    }
    for (int i = 3; i < 9; i += 2) {
        test_value(arf2.at(i), 0.0);
        test_value(arf2[i], 0.0);
    }
    std::cout << arf << std::endl;
    std::cout << arf2 << std::endl;

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("Xspec module");

    // Initially assume that we pass all tests
    bool success = true;

    // Create a test suite
    TestGXspec test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GXspec.xml");

    // Return success status
    return (success ? 0 : 1);
}
