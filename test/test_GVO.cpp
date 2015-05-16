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
#include <unistd.h>      // for sleep() function
#include <pthread.h>
#include <signal.h>      // for pthread_kill() function
#include "test_GVO.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Thread start routine to test GVOHub class
 *
 * @param[in] ptr Start routine argument.
 **************************************************************************/
void* vo_thread(void* ptr)
{
    // Create VO Hub
    GVOHub hub;

    // Start Hub
    hub.start();

    // Return
    return NULL;
}


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
    // Declare thread
    pthread_t thread;
std::cout << "\n*** Before start of thread ***" << std::endl;
system("netstat -natp tcp");
    
    // Start Hub in a thread
    test_try("Start hub");
    try {
        int rc = pthread_create(&thread, NULL, vo_thread, NULL);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
std::cout << "\n*** After start of thread ***" << std::endl;
system("netstat -natp tcp");

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

    // Cancel thread. We do this here so that for the next test
    // there is no Hub anymore. We do some sleeping to assure that
    // the thread has been canceled.
    pthread_cancel(thread);
    sleep(1);
std::cout << "\n*** After thread canel ***" << std::endl;
system("netstat -natp tcp");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test VO client class
 **************************************************************************/
void TestGVO::test_GVOClient(void)
{
    // Test constructor
    test_try("GVOClient empty constructor");
    try {
        GVOClient client;
        //std::cout << client << std::endl;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

std::cout << "\n*** Before client connect ***" << std::endl;
system("netstat -natp tcp");
    // Create client
    GVOClient client;

    // Connect client
    client.connect();
    test_assert(client.is_connected(), "Check for connection.");
    //std::cout << client << std::endl;
std::cout << "\n*** After client connect ***" << std::endl;
system("netstat -natp tcp");

    // Test ping
    test_assert(client.ping_hub(), "Ping Hub.");

    // Disconnect client
    client.disconnect();
    test_assert(!client.is_connected(), "Check for disconnection.");
    //std::cout << client << std::endl;

    // Test ping
    test_assert(client.ping_hub(), "Ping Hub.");

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
