/***************************************************************************
 *    test_GSupport.cpp  -  test GSupport classes and GTools functions     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Jurgen Knodlseder                           *
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
#include "test_GSupport.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/***********************************************************************//**
 * @brief Set parameters and tests
 ***************************************************************************/
void TestGSupport::set(void){
    // Test name
    name("GSupport");


    //add tests
    add_test(static_cast<pfunction>(&TestGSupport::test_expand_env),"Test Environment variable");

    return;
}

/***********************************************************************//**
 * @brief Test expand_env function
 *
 * Test the environment variable expansion function.
 ***************************************************************************/
void TestGSupport::test_expand_env(void)
{

    // Declare strings that we use during the tests
    std::string s_in;
    std::string s_ref;
    std::string s_out;

    // Test 1: no environment variable
    s_in  = "This string has no environment variable.";
    s_ref = s_in;
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"no environment variable","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Get home environment variable
    std::string home(std::getenv("HOME"));

    // $ENV{HOME} environment variable within string
    s_in  = "My $ENV{HOME} is my castle.";
    s_ref = "My "+home+" is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"$ENV{HOME} environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // ${HOME} environment variable within string
    s_in  = "My ${HOME} is my castle.";
    s_ref = "My "+home+" is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"${HOME} environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $HOME/path environment variable within string
    s_in  = "My $HOME/path is my castle.";
    s_ref = "My "+home+"/path is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"$HOME/path environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $HOME environment variable at end of string
    s_in  = "My $HOME";
    s_ref = "My "+home;
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"$HOME environment variable at end of string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Environment variable within single quotes
    s_in  = "My '$(HOME)' is my castle.";
    s_ref = s_in;
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"Environment variable within single quotes","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Environment variable within double quotes
    s_in  = "My \"$(HOME)\" is my castle.";
    s_ref = "My \""+home+"\" is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"Environment variable within double quotes","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Non existing environment variable within string
    s_in  = "My $(HOMEIXYZ) is my castle.";
    s_ref = "My  is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"Non existing environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Empty environment variable within string
    s_in  = "My $() is my castle.";
    s_ref = "My  is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"Empty environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Two successive environment variables within string
    s_in  = "My $(HOME)$ENV{HOME} is my castle.";
    s_ref = "My "+home+home+" is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"Two successive environment variables within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $(HOME) at beginning of string
    s_in  = "$(HOME) is my castle.";
    s_ref = home+" is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"$(HOME) at beginning of string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $ENV(HOME) at beginning of string
    s_in  = "$ENV(HOME) is my castle.";
    s_ref = home+" is my castle.";
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"$ENV(HOME) at beginning of string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $(HOME) at end of string
    s_in  = "My castle is $(HOME)";
    s_ref = "My castle is "+home;
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"$(HOME) at end of string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $ENV{HOME}${HOME}$ENV(HOME)$(HOME) only string
    s_in  = "$ENV{HOME}${HOME}$ENV(HOME)$(HOME)";
    s_ref = home+home+home+home;
    s_out = expand_env(s_in);
    test_assert(s_out == s_ref,"$ENV{HOME}${HOME}$ENV(HOME)$(HOME) only string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Debugging
    //std::cout << std::endl;
    //std::cout << s_in << std::endl;
    //std::cout << s_ref << std::endl;
    //std::cout << s_out << std::endl;

    // Exit test
    return;
}

/***********************************************************************//**
 * @brief Main test entry point
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuites("GSupport");

    bool was_successful=true;

    //Create a test suite
    TestGSupport test;

    //Append to the container
    testsuites.append(test);

    //Run
    was_successful=testsuites.run();

    //save xml report
    testsuites.save("reports/GSupport.xml");

    // Return
    return was_successful ? 0:1;
}
