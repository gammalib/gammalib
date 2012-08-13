/***************************************************************************
 *    test_GSupport.cpp  -  test GSupport classes and GTools functions     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file test_GSupport.cpp
 * @brief Testing of support module
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception
#include "GammaLib.hpp"
#include "GTools.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/***********************************************************************//**
 * @brief Test expand_env function
 *
 * Test the environment variable expansion function.
 ***************************************************************************/
void test_expand_env(GTestSuite& testsuite)
{
    testsuite.test_try("Test Environment variable");
    try {

        // Declare strings that we use during the tests
        std::string s_in;
        std::string s_ref;
        std::string s_out;

        // Test 1: no environment variable
        testsuite.test_try("no environment variable");
        try {
            s_in  = "This string has no environment variable.";
            s_ref = s_in;
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // Get home environment variable
        std::string home(std::getenv("HOME"));

        // $ENV{HOME} environment variable within string
        testsuite.test_try("$ENV{HOME} environment variable within string");
        try {
            s_in  = "My $ENV{HOME} is my castle.";
            s_ref = "My "+home+" is my castle.";
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // ${HOME} environment variable within string
        testsuite.test_try("${HOME} environment variable within string");
        try {
            s_in  = "My ${HOME} is my castle.";
            s_ref = "My "+home+" is my castle.";
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // $HOME/path environment variable within string
        testsuite.test_try("$HOME/path environment variable within string");
        try {
            s_in  = "My $HOME/path is my castle.";
            s_ref = "My "+home+"/path is my castle.";
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // $HOME environment variable at end of string
        testsuite.test_try("$HOME environment variable at end of string");
        try {
        s_in  = "My $HOME";
        s_ref = "My "+home;
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
        }

        testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // Environment variable within single quotes
        testsuite.test_try("Environment variable within single quotes");
        try {
            s_in  = "My '$(HOME)' is my castle.";
            s_ref = s_in;
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // Environment variable within double quotes
        testsuite.test_try("Environment variable within double quotes");
        try {
           s_in  = "My \"$(HOME)\" is my castle.";
            s_ref = "My \""+home+"\" is my castle.";
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // Non existing environment variable within string
        testsuite.test_try("Non existing environment variable within string");
        try {
        s_in  = "My $(HOMEIXYZ) is my castle.";
        s_ref = "My  is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
        }

        testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // Empty environment variable within string
        testsuite.test_try("Empty environment variable within string");
        try {
            s_in  = "My $() is my castle.";
            s_ref = "My  is my castle.";
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }
    
            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // Two successive environment variables within string
        testsuite.test_try("Two successive environment variables within string");
        try {
            s_in  = "My $(HOME)$ENV{HOME} is my castle.";
            s_ref = "My "+home+home+" is my castle.";
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }
        
            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // $(HOME) at beginning of string
        testsuite.test_try("$(HOME) at beginning of string");
        try {
            s_in  = "$(HOME) is my castle.";
            s_ref = home+" is my castle.";
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }
            
            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // $ENV(HOME) at beginning of string
        testsuite.test_try("$ENV(HOME) at beginning of string");
        try {
            s_in  = "$ENV(HOME) is my castle.";
            s_ref = home+" is my castle.";
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // $(HOME) at end of string
        testsuite.test_try("$(HOME) at end of string");
        try {
            s_in  = "My castle is $(HOME)";
            s_ref = "My castle is "+home;
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // $ENV{HOME}${HOME}$ENV(HOME)$(HOME) only string
        testsuite.test_try("$ENV{HOME}${HOME}$ENV(HOME)$(HOME) only string");
        try {
            s_in  = "$ENV{HOME}${HOME}$ENV(HOME)$(HOME)";
            s_ref = home+home+home+home;
            s_out = expand_env(s_in);
            if (s_out != s_ref) {
                throw testsuite.exception_failure("Unexpected string '"+s_out+"' (expected '"+s_ref+"')");
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        // Debugging
        //std::cout << std::endl;
        //std::cout << s_in << std::endl;
        //std::cout << s_ref << std::endl;
        //std::cout << s_out << std::endl;

    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Error occured while testing "
                  << "environment variable expansion function"
                  << "." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    
    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GTools
 *
 * Test GTools support functions.
 ***************************************************************************/
void test_GTools(GTestSuite& testsuite)
{
    // Perform tests
    test_expand_env(testsuite);

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Main test entry point
 ***************************************************************************/
int main(void)
{
    // Create a test suite container
    GTestSuites testsuites("GammaLib support module testing");

    // Create a test suite
    GTestSuite testsuite("Test GTools");

    // Append testsuite
    testsuites.append(testsuite);

    // Execute tests
    test_GTools(testsuite);

    //End of the test suite
    testsuite.end_test();

    //save xml report
    testsuites.save("reports/GSupport.xml");

    // Return
    return 0;
}
