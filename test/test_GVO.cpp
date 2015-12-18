/***************************************************************************
 *                   test_GVO.cpp - Test VO module                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2015 by Juergen Knoedlseder                         *
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
 * @file test_GVO.hpp
 * @brief Implementation of unit tests for VO module
 * @author Juergen Knoedlseder 
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>       // for system() function
#include <unistd.h>      // for sleep() and usleep() function
#include "test_GVO.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGVO::set(void)
{
    // Test name
    name("GVO");

    // Append tests
    append(static_cast<pfunction>(&TestGVO::test_GVOHub), "Test GVOHub class");
    append(static_cast<pfunction>(&TestGVO::test_GVOClient), "Test GVOClient class");

    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGVO* TestGVO::clone(void) const
{
    // Clone test suite
    return new TestGVO(*this);
}


/***********************************************************************//**
 * @brief Test VO hub class
 **************************************************************************/
void TestGVO::test_GVOHub(void)
{
    // Kluge: wait until no Hub is online
    GRan ran;
    int  no_hub = 0;
    while (no_hub < 3) {
        GVOClient client;
        if (!client.ping_hub()) {
            no_hub++;
            usleep(ran.uniform()*1.0e6);
        }
        else {
            no_hub = 0;
        }
        usleep(ran.uniform()*1.0e6);
    }

   	// Create child process
    int pid = fork();

    // If we have a PID of 0 we are in the child process
    if (pid == 0) {

        // Create VO Hub
        GVOHub hub;

        // Start Hub
        hub.start();
            
        // Exit child process
        exit(0);
    
    }

    // ... otherwise if PID is >0 then we are in the parent process
    else if (pid > 0) {
    
        // Test Hub connection
        GVOClient client;
        for (int i = 0; i < 3; ++i) {
            sleep(1);
            client.connect();
            if (client.is_connected()) {
                break;
            }
        }

        // Test ping
        test_assert(client.ping_hub(), "Ping Hub.");

        // Disconnect client
        client.disconnect();

        // We now shutdown the Hub via the client
        client.shutdown_hub();

        // Sleep a bit to make sure that we exit without a Hub
        sleep(1);

    } // endelse: parent process

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test VO client class
 **************************************************************************/
void TestGVO::test_GVOClient(void)
{
    // Kluge: wait until no Hub is online
    GRan ran;
    int  no_hub = 0;
    while (no_hub < 3) {
        GVOClient client;
        if (!client.ping_hub()) {
            no_hub++;
            usleep(ran.uniform()*1.0e6);
        }
        else {
            no_hub = 0;
        }
        usleep(ran.uniform()*1.0e6);
    }

    // Test constructor
    test_try("GVOClient empty constructor");
    try {
        GVOClient client;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Create client
    GVOClient client;

    // Test ping
    test_assert(!client.ping_hub(), "Ping Hub (should not be alive).");

    // Connect client
    client.connect();
    test_assert(client.is_connected(), "Check for connection.");

    // Test ping
    test_assert(client.ping_hub(), "Ping Hub.");

    // Disconnect client
    client.disconnect();
    test_assert(!client.is_connected(), "Check for disconnection.");

    // Test ping
    test_assert(client.ping_hub(), "Ping Hub.");

    // Shutdown VO Hub
    client.shutdown_hub();

    // Sleep a bit to make sure that the Hub is dead
    sleep(1);

    // Test ping
    test_assert(!client.ping_hub(), "Ping Hub (should not be alive).");

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("VO module");

    // Unset HOME environment variable (so that the .samp file is written
    // into the test directory)
    unsetenv("HOME");

    // Initially assume that we pass all tests
    bool success = true;

    // Create a test suite
    TestGVO test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GVO.xml");

    // Return success status
    return (success ? 0 : 1);
}
