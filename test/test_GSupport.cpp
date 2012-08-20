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
void test_expand_env(void)
{
    // Main test protection
    try {
    
        // Declare strings that we use during the tests
        std::string s_in;
        std::string s_ref;
        std::string s_out;
        
        // Test 1: no environment variable
        s_in  = "This string has no environment variable.";
        s_ref = s_in;
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";
        
        // Get home environment variable
        std::string home(std::getenv("HOME"));

        // $ENV{HOME} environment variable within string
        s_in  = "My $ENV{HOME} is my castle.";
        s_ref = "My "+home+" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $ENV(HOME) environment variable within string
        s_in  = "My $ENV(HOME) is my castle.";
        s_ref = "My "+home+" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // ${HOME} environment variable within string
        s_in  = "My ${HOME} is my castle.";
        s_ref = "My "+home+" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $(HOME) environment variable within string
        s_in  = "My $(HOME) is my castle.";
        s_ref = "My "+home+" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $HOME environment variable within string
        s_in  = "My $HOME is my castle.";
        s_ref = "My "+home+" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $HOME/path environment variable within string
        s_in  = "My $HOME/path is my castle.";
        s_ref = "My "+home+"/path is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $HOME environment variable at end of string
        s_in  = "My $HOME";
        s_ref = "My "+home;
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // Environment variable within single quotes
        s_in  = "My '$(HOME)' is my castle.";
        s_ref = s_in;
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // Environment variable within double quotes
        s_in  = "My \"$(HOME)\" is my castle.";
        s_ref = "My \""+home+"\" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // Non existing environment variable within string
        s_in  = "My $(HOMEIXYZ) is my castle.";
        s_ref = "My  is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // Empty environment variable within string
        s_in  = "My $() is my castle.";
        s_ref = "My  is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // Two successive environment variables within string
        s_in  = "My $(HOME)$ENV{HOME} is my castle.";
        s_ref = "My "+home+home+" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $(HOME) at beginning of string
        s_in  = "$(HOME) is my castle.";
        s_ref = home+" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $ENV(HOME) at beginning of string
        s_in  = "$ENV(HOME) is my castle.";
        s_ref = home+" is my castle.";
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $(HOME) at end of string
        s_in  = "My castle is $(HOME)";
        s_ref = "My castle is "+home;
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

        // $ENV{HOME}${HOME}$ENV(HOME)$(HOME) only string
        s_in  = "$ENV{HOME}${HOME}$ENV(HOME)$(HOME)";
        s_ref = home+home+home+home;
        s_out = expand_env(s_in);
        if (s_out != s_ref) {
            std::cout << std::endl 
                      << "TEST ERROR: Unexpected string "
                      << "'"+s_out+"' "
                      << "(expected '"+s_ref+"')" << std::endl;
            throw;
        }
        std::cout << ".";

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
void test_GTools(void)
{
    // Dump header
    std::cout << "Test GTools: ";

    // Test GTools
    try {
    
        // Perform tests
        test_expand_env();

    }
    catch (std::exception &e) {
        std::cout << std::endl 
                  << "TEST ERROR: Error occured while testing "
                  << "GTools functions"
                  << "." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }

    // Signal final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Main test entry point
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "***********************************" << std::endl;
    std::cout << "* GammaLib support module testing *" << std::endl;
    std::cout << "***********************************" << std::endl;

    // Execute tests
    test_GTools();

    // Return
    return 0;
}
