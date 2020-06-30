/***************************************************************************
 *                  test_GXspec.hpp - Test Xspec module                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2020 by Juergen Knoedlseder                         *
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
#include <cstdlib>     // getenv
#include "test_GXspec.hpp"
#include "GTools.hpp"


/* __ Constants __________________________________________________________ */
const std::string datadir = std::getenv("TEST_DATA");
const std::string rmfname = datadir + "/rmf.fits";


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
    append(static_cast<pfunction>(&TestGXspec::test_GRmf), "Test GRmf class");

    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGXspec* TestGXspec::clone(void) const
{
    // Clone test suite
    return new TestGXspec(*this);
}


/***********************************************************************//**
 * @brief Test GPha class
 **************************************************************************/
void TestGXspec::test_GPha(void)
{
    // Test void constructor
    GPha pha1;
    test_value(pha1.size(), 0, "Check if empty PHA has zero size");
    test_value(pha1.columns(), 0, "Check if empty PHA has zero columns");
    test_value(pha1.counts(), 0.0, "Check if empty PHA has zero counts");
    test_value(pha1.underflow(), 0.0, "Check if empty PHA has zero underflow");
    test_value(pha1.overflow(), 0.0, "Check if empty PHA has zero overflow");
    test_value(pha1.outflow(), 0.0, "Check if empty PHA has zero outflow");
    test_value(pha1.exposure(), 0.0, "Check if empty PHA has zero exposure");
    test_value(pha1.backfile(), "", "Check if empty PHA has empty background file");
    test_value(pha1.corrfile(), "", "Check if empty PHA has empty correction file");
    test_value(pha1.respfile(), "", "Check if empty PHA has empty RMF");
    test_value(pha1.ancrfile(), "", "Check if empty PHA has empty ARF");

    // Test energy boundary constructor, filling and accessing
    GEbounds ebds2(9, GEnergy(1.0, "TeV"), GEnergy(10.0, "TeV"), "LIN");
    GPha     pha2(ebds2);
    pha2.exposure(1234.5);
    pha2.backfile("background");
    pha2.corrfile("correction");
    pha2.respfile("rmf");
    pha2.ancrfile("arf");
    test_value(pha2.size(), 9, "Check if non-empty PHA has expected size");
    test_value(pha2.exposure(), 1234.5, "Check if non-empty PHA has expected exposure");
    test_value(pha2.backfile(), "background", "Check if non-empty PHA has expected background file");
    test_value(pha2.corrfile(), "correction", "Check if non-empty PHA has expected correction file");
    test_value(pha2.respfile(), "rmf", "Check if non-empty PHA has expected RMF");
    test_value(pha2.ancrfile(), "arf", "Check if non-empty PHA has expected ARF");
    for (int i = 0; i < 12; i += 2) {
        pha2.fill(GEnergy(double(i), "TeV"), 1.0);
    }
    test_value(pha2.counts(),    4.0);
    test_value(pha2.underflow(), 1.0);
    test_value(pha2.overflow(),  1.0);
    test_value(pha2.outflow(),   0.0);
    for (int i = 0; i < 9; i += 2) {
        test_value(pha2.at(i), 0.0);
        test_value(pha2[i], 0.0);
    }
    for (int i = 1; i < 9; i += 2) {
        test_value(pha2.at(i), 1.0);
        test_value(pha2[i], 1.0);
    }
    pha2[0] = 5.0;
    pha2[1] = 3.7;
    test_value(pha2[0], 5.0);
    test_value(pha2[1], 3.7);

    // Test content after saving and loading
    test_value(pha2.filename().url(), "", "Unexpected filename \""+pha2.filename().url()+"\".");
    pha2.save("test_pha.fits", true);
    test_value(pha2.filename().url(), "test_pha.fits",
               "Unexpected filename \""+pha2.filename().url()+"\".");
    pha2.clear();
    test_value(pha2.size(), 0, "Check if cleared PHA has zero size");
    test_value(pha2.counts(), 0.0, "Check if clear PHA has zero counts");
    test_value(pha2.underflow(), 0.0, "Check if empty PHA has zero underflow");
    test_value(pha2.overflow(), 0.0, "Check if empty PHA has zero overflow");
    test_value(pha2.outflow(), 0.0, "Check if empty PHA has zero outflow");
    test_value(pha2.exposure(), 0.0, "Check if empty PHA has zero exposure");
    test_value(pha2.backfile(), "", "Check if empty PHA has empty background file");
    test_value(pha2.corrfile(), "", "Check if empty PHA has empty correction file");
    test_value(pha2.respfile(), "", "Check if empty PHA has empty RMF");
    test_value(pha2.ancrfile(), "", "Check if empty PHA has empty ARF");
    pha2.load("test_pha.fits");
    test_value(pha2.filename().url(), "test_pha.fits",
               "Unexpected filename \""+pha2.filename().url()+"\".");
    test_value(pha2.size(), 9, "Check if loaded PHA has expected size");
    test_value(pha2.exposure(), 1234.5, "Check if loaded PHA has expected exposure");
    test_value(pha2.backfile(), "background", "Check if loaded PHA has expected background file");
    test_value(pha2.corrfile(), "correction", "Check if loaded PHA has expected correction file");
    test_value(pha2.respfile(), "rmf", "Check if loaded PHA has expected RMF");
    test_value(pha2.ancrfile(), "arf", "Check if loaded PHA has expected ARF");
    test_value(pha2[0], 5.0, 1.0e-6);
    test_value(pha2[1], 3.7, 1.0e-6);
    test_value(pha2.counts(),   11.7, 1.0e-6);
    test_value(pha2.underflow(), 1.0, 1.0e-6);
    test_value(pha2.overflow(),  1.0, 1.0e-6);
    test_value(pha2.outflow(),   0.0, 1.0e-6);
    for (int i = 2; i < 9; i += 2) {
        test_value(pha2.at(i), 0.0);
        test_value(pha2[i], 0.0);
    }
    for (int i = 3; i < 9; i += 2) {
        test_value(pha2.at(i), 1.0);
        test_value(pha2[i], 1.0);
    }

    // Test constructing from file
    GPha pha4("test_pha.fits");
    test_value(pha4.filename().url(), "test_pha.fits",
               "Unexpected filename \""+pha4.filename().url()+"\".");
    test_value(pha4.size(), 9, "Check if loaded PHA has expected size");
    test_value(pha4.exposure(), 1234.5, "Check if constructed PHA has expected exposure");
    test_value(pha4.backfile(), "background", "Check if constructed PHA has expected background file");
    test_value(pha4.corrfile(), "correction", "Check if constructed PHA has expected correction file");
    test_value(pha4.respfile(), "rmf", "Check if constructed PHA has expected RMF");
    test_value(pha4.ancrfile(), "arf", "Check if constructed PHA has expected ARF");
    test_value(pha4[0], 5.0, 1.0e-6);
    test_value(pha4[1], 3.7, 1.0e-6);
    test_value(pha4.counts(),   11.7, 1.0e-6);
    test_value(pha4.underflow(), 1.0, 1.0e-6);
    test_value(pha4.overflow(),  1.0, 1.0e-6);
    test_value(pha4.outflow(),   0.0, 1.0e-6);
    for (int i = 2; i < 9; i += 2) {
        test_value(pha4.at(i), 0.0);
        test_value(pha4[i], 0.0);
    }
    for (int i = 3; i < 9; i += 2) {
        test_value(pha4.at(i), 1.0);
        test_value(pha4[i], 1.0);
    }

    // Test addition operator
    GPha pha5 = pha4 + pha4;
    test_value(pha5[0], 10.0, 1.0e-6);
    test_value(pha5[1],  7.4, 1.0e-6);
    test_value(pha5.counts(),   23.4, 1.0e-6);
    test_value(pha5.underflow(), 2.0, 1.0e-6);
    test_value(pha5.overflow(),  2.0, 1.0e-6);
    test_value(pha5.outflow(),   0.0, 1.0e-6);

    // Test subtraction operator
    GPha pha6 = pha4 - pha4;
    test_value(pha6[0], 0.0, 1.0e-6);
    test_value(pha6[1], 0.0, 1.0e-6);
    test_value(pha6.counts(),    0.0, 1.0e-6);
    test_value(pha6.underflow(), 0.0, 1.0e-6);
    test_value(pha6.overflow(),  0.0, 1.0e-6);
    test_value(pha6.outflow(),   0.0, 1.0e-6);

    // Test scale operator
    GPha pha7 = pha4 * 2.0;
    test_value(pha7[0], 10.0, 1.0e-6);
    test_value(pha7[1],  7.4, 1.0e-6);
    test_value(pha7.counts(),   23.4, 1.0e-6);
    test_value(pha7.underflow(), 2.0, 1.0e-6);
    test_value(pha7.overflow(),  2.0, 1.0e-6);
    test_value(pha7.outflow(),   0.0, 1.0e-6);

    // Test scale operator
    GPha pha8 = 2.0 * pha4;
    test_value(pha8[0], 10.0, 1.0e-6);
    test_value(pha8[1],  7.4, 1.0e-6);
    test_value(pha8.counts(),   23.4, 1.0e-6);
    test_value(pha8.underflow(), 2.0, 1.0e-6);
    test_value(pha8.overflow(),  2.0, 1.0e-6);
    test_value(pha8.outflow(),   0.0, 1.0e-6);

    // Test division operator
    GPha pha9 = pha4 / 2.0;
    test_value(pha9[0],  2.5, 1.0e-6);
    test_value(pha9[1],  1.85, 1.0e-6);
    test_value(pha9.counts(),   5.85, 1.0e-6);
    test_value(pha9.underflow(), 0.5, 1.0e-6);
    test_value(pha9.overflow(),  0.5, 1.0e-6);
    test_value(pha9.outflow(),   0.0, 1.0e-6);

    // Test unary addition operator
    GPha pha10 = pha4;
    pha10     += pha4;
    test_value(pha10[0], 10.0, 1.0e-6);
    test_value(pha10[1],  7.4, 1.0e-6);
    test_value(pha10.counts(),   23.4, 1.0e-6);
    test_value(pha10.underflow(), 2.0, 1.0e-6);
    test_value(pha10.overflow(),  2.0, 1.0e-6);
    test_value(pha10.outflow(),   0.0, 1.0e-6);

    // Test unary subtraction operator
    GPha pha11 = pha4;
    pha11     -= pha4;
    test_value(pha11[0], 0.0, 1.0e-6);
    test_value(pha11[1], 0.0, 1.0e-6);
    test_value(pha11.counts(),    0.0, 1.0e-6);
    test_value(pha11.underflow(), 0.0, 1.0e-6);
    test_value(pha11.overflow(),  0.0, 1.0e-6);
    test_value(pha11.outflow(),   0.0, 1.0e-6);

    // Test unary scale operator
    GPha pha12 = pha4;
    pha12     *= 2.0;
    test_value(pha12[0], 10.0, 1.0e-6);
    test_value(pha12[1],  7.4, 1.0e-6);
    test_value(pha12.counts(),   23.4, 1.0e-6);
    test_value(pha12.underflow(), 2.0, 1.0e-6);
    test_value(pha12.overflow(),  2.0, 1.0e-6);
    test_value(pha12.outflow(),   0.0, 1.0e-6);

    // Test unary division operator
    GPha pha13 = pha4;
    pha13     /= 2.0;
    test_value(pha13[0],  2.5, 1.0e-6);
    test_value(pha13[1],  1.85, 1.0e-6);
    test_value(pha13.counts(),   5.85, 1.0e-6);
    test_value(pha13.underflow(), 0.5, 1.0e-6);
    test_value(pha13.overflow(),  0.5, 1.0e-6);
    test_value(pha13.outflow(),   0.0, 1.0e-6);

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
    test_assert(arf.filename().url() == "arf.fits",
                "Unexpected filename \""+arf.filename().url()+"\".");
    arf.load("arf.fits");
    test_assert(arf.filename().url() == "arf.fits",
                "Unexpected filename \""+arf.filename().url()+"\".");
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
    test_assert(arf2.filename().url() == "arf.fits",
                "Unexpected filename \""+arf2.filename().url()+"\".");
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

    // Test addition operator
    GArf arf3 = arf2 + arf2;
    test_value(arf3[0], 10.0, 1.0e-6);
    test_value(arf3[1],  7.4, 1.0e-6);

    // Test subtraction operator
    GArf arf4 = arf2 - arf2;
    test_value(arf4[0], 0.0, 1.0e-6);
    test_value(arf4[1], 0.0, 1.0e-6);

    // Test scale operator
    GArf arf5 = arf2 * 2.0;
    test_value(arf5[0], 10.0, 1.0e-6);
    test_value(arf5[1],  7.4, 1.0e-6);

    // Test scale operator
    GArf arf6 = 2.0 * arf2;
    test_value(arf6[0], 10.0, 1.0e-6);
    test_value(arf6[1],  7.4, 1.0e-6);

    // Test division operator
    GArf arf7 = arf2 / 2.0;
    test_value(arf7[0], 2.5, 1.0e-6);
    test_value(arf7[1], 1.85, 1.0e-6);

    // Test unary addition operator
    GArf arf8 = arf2;
    arf8 += arf2;
    test_value(arf8[0], 10.0, 1.0e-6);
    test_value(arf8[1],  7.4, 1.0e-6);

    // Test unary subtraction operator
    GArf arf9 = arf2;
    arf9 -= arf2;
    test_value(arf9[0], 0.0, 1.0e-6);
    test_value(arf9[1], 0.0, 1.0e-6);

    // Test unary scale operator
    GArf arf10 = arf2;
    arf10 *= 2.0;
    test_value(arf10[0], 10.0, 1.0e-6);
    test_value(arf10[1],  7.4, 1.0e-6);

    // Test unary divisison operator
    GArf arf11 = arf2;
    arf11 /= 2.0;
    test_value(arf11[0], 2.5, 1.0e-6);
    test_value(arf11[1], 1.85, 1.0e-6);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GRmf class
 **************************************************************************/
void TestGXspec::test_GRmf(void)
{
    // Test void constructor
    test_try("GRmf void constructor");
    try {
        GRmf rmf;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test filename constructor
    test_try("GRmf filename constructor");
    try {
        GRmf rmf(rmfname);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test energy boundary constructor
    test_try("GRmf energy boundary constructor");
    try {
        GEbounds ebds(10, GEnergy(0.1, "TeV"), GEnergy(10.0, "TeV"));
        GRmf     rmf(ebds, ebds);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test filling and accessing
    GEbounds ebds(9, GEnergy(1.0, "TeV"), GEnergy(10.0, "TeV"));
    GRmf     rmf(ebds, ebds);
    for (int i = 0; i < 9; ++i) {
        for (int k = i; k < i+3 && k < 9; ++k) {
            rmf(i, k) = 1.0;
        }
    }
    for (int i = 0; i < 9; ++i) {
        int k = 0;
        for (; k < i; ++k) {
            test_value(rmf.at(i,k), 0.0);
            test_value(rmf(i,k),    0.0);
        }
        for (; k < i+3 && k < 9; ++k) {
            test_value(rmf.at(i,k), 1.0);
            test_value(rmf(i,k),    1.0);
        }
        for (; k < 9; ++k) {
            test_value(rmf.at(i,k), 0.0);
            test_value(rmf(i,k),    0.0);
        }
    }

    // Test saving and loading
    rmf.save("rmf.fits", true);
    test_assert(rmf.filename().url() == "rmf.fits",
                "Unexpected filename \""+rmf.filename().url()+"\".");
    rmf.load("rmf.fits");
    test_assert(rmf.filename().url() == "rmf.fits",
                "Unexpected filename \""+rmf.filename().url()+"\".");
    for (int i = 0; i < 9; ++i) {
        int k = 0;
        for (; k < i; ++k) {
            test_value(rmf.at(i,k), 0.0);
            test_value(rmf(i,k),    0.0);
        }
        for (; k < i+3 && k < 9; ++k) {
            test_value(rmf.at(i,k), 1.0);
            test_value(rmf(i,k),    1.0);
        }
        for (; k < 9; ++k) {
            test_value(rmf.at(i,k), 0.0);
            test_value(rmf(i,k),    0.0);
        }
    }

    // Test constructing
    GRmf rmf2("rmf.fits");
    test_assert(rmf2.filename().url() == "rmf.fits",
                "Unexpected filename \""+rmf2.filename().url()+"\".");
    for (int i = 0; i < 9; ++i) {
        int k = 0;
        for (; k < i; ++k) {
            test_value(rmf.at(i,k), 0.0);
            test_value(rmf(i,k),    0.0);
        }
        for (; k < i+3 && k < 9; ++k) {
            test_value(rmf.at(i,k), 1.0);
            test_value(rmf(i,k),    1.0);
        }
        for (; k < 9; ++k) {
            test_value(rmf.at(i,k), 0.0);
            test_value(rmf(i,k),    0.0);
        }
    }

    // Test addition operator
    GRmf rmf3 = rmf2 + rmf2;
    test_value(rmf3(0,0), 2.0, 1.0e-6);
    test_value(rmf3(1,1), 2.0, 1.0e-6);

    // Test subtraction operator
    GRmf rmf4 = rmf2 - rmf2;
    test_value(rmf4(0,0), 0.0, 1.0e-6);
    test_value(rmf4(1,1), 0.0, 1.0e-6);

    // Test scale operator
    GRmf rmf5 = rmf2 * 2.0;
    test_value(rmf5(0,0), 2.0, 1.0e-6);
    test_value(rmf5(1,1), 2.0, 1.0e-6);

    // Test scale operator
    GRmf rmf6 = 2.0 * rmf2;
    test_value(rmf6(0,0), 2.0, 1.0e-6);
    test_value(rmf6(1,1), 2.0, 1.0e-6);

    // Test division operator
    GRmf rmf7 = rmf2 / 2.0;
    test_value(rmf7(0,0), 0.5, 1.0e-6);
    test_value(rmf7(1,1), 0.5, 1.0e-6);

    // Test unary addition operator
    GRmf rmf8 = rmf2;
    rmf8 += rmf2;
    test_value(rmf8(0,0), 2.0, 1.0e-6);
    test_value(rmf8(1,1), 2.0, 1.0e-6);

    // Test unary subtraction operator
    GRmf rmf9 = rmf2;
    rmf9 -= rmf2;
    test_value(rmf9(0,0), 0.0, 1.0e-6);
    test_value(rmf9(1,1), 0.0, 1.0e-6);

    // Test unary scale operator
    GRmf rmf10 = rmf2;
    rmf10 *= 2.0;
    test_value(rmf10(0,0), 2.0, 1.0e-6);
    test_value(rmf10(1,1), 2.0, 1.0e-6);

    // Test unary divisison operator
    GRmf rmf11 = rmf2;
    rmf11 /= 2.0;
    test_value(rmf11(0,0), 0.5, 1.0e-6);
    test_value(rmf11(1,1), 0.5, 1.0e-6);

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
