/***************************************************************************
 *                 test_GSupport.cpp - test support module                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>   // getenv
#include <vector>
#include "GTools.hpp"
#include "test_GSupport.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/***********************************************************************//**
 * @brief Set parameters and tests
 ***************************************************************************/
void TestGSupport::set(void){

    // Set test name
    name("GSupport");

    // Add tests
    add_test(static_cast<pfunction>(&TestGSupport::test_expand_env), "Test Environment variable");
    add_test(static_cast<pfunction>(&TestGSupport::test_node_array), "Test GNodeArray");
    add_test(static_cast<pfunction>(&TestGSupport::test_url_file),   "Test GUrlFile");
    add_test(static_cast<pfunction>(&TestGSupport::test_url_string), "Test GUrlString");

    // Return
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
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"no environment variable","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Get home environment variable
    std::string home(std::getenv("HOME"));

    // $ENV{HOME} environment variable within string
    s_in  = "My $ENV{HOME} is my castle.";
    s_ref = "My "+home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$ENV{HOME} environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // ${HOME} environment variable within string
    s_in  = "My ${HOME} is my castle.";
    s_ref = "My "+home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"${HOME} environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $HOME/path environment variable within string
    s_in  = "My $HOME/path is my castle.";
    s_ref = "My "+home+"/path is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$HOME/path environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $HOME environment variable at end of string
    s_in  = "My $HOME";
    s_ref = "My "+home;
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$HOME environment variable at end of string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Environment variable within single quotes
    s_in  = "My '$(HOME)' is my castle.";
    s_ref = s_in;
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Environment variable within single quotes","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Environment variable within double quotes
    s_in  = "My \"$(HOME)\" is my castle.";
    s_ref = "My \""+home+"\" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Environment variable within double quotes","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Non existing environment variable within string
    s_in  = "My $(HOMEIXYZ) is my castle.";
    s_ref = "My  is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Non existing environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Empty environment variable within string
    s_in  = "My $() is my castle.";
    s_ref = "My  is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Empty environment variable within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // Two successive environment variables within string
    s_in  = "My $(HOME)$ENV{HOME} is my castle.";
    s_ref = "My "+home+home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Two successive environment variables within string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $(HOME) at beginning of string
    s_in  = "$(HOME) is my castle.";
    s_ref = home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$(HOME) at beginning of string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $ENV(HOME) at beginning of string
    s_in  = "$ENV(HOME) is my castle.";
    s_ref = home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$ENV(HOME) at beginning of string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $(HOME) at end of string
    s_in  = "My castle is $(HOME)";
    s_ref = "My castle is "+home;
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$(HOME) at end of string","Unexpected string '"+s_out+"' (expected '"+s_ref+"')");

    // $ENV{HOME}${HOME}$ENV(HOME)$(HOME) only string
    s_in  = "$ENV{HOME}${HOME}$ENV(HOME)$(HOME)";
    s_ref = home+home+home+home;
    s_out = gammalib::expand_env(s_in);
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
 * @brief Test GNodeArray class
 *
 * Test the GNodeArray class.
 ***************************************************************************/
void TestGSupport::test_node_array(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GNodeArray array;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test copy constructor
    test_try("Copy constructor");
    try {
        GNodeArray array;
        GNodeArray array2(array);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test array constructor
    test_try("Array constructor");
    try {
        double array[] = {0.0, 1.0, 2.0, 3.0, 5.0};
        GNodeArray array2(5, array);
        test_try_success();
        for (int i = 0; i < 5; ++i) {
            test_value(array2[i], array[i]);
        }
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GVector constructor
    test_try("GVector constructor");
    try {
        double array[] = {0.0, 1.0, 2.0, 3.0, 5.0};
        GVector vector(5);
        for (int i = 0; i < 5; ++i) {
            vector[i] = array[i];
        }
        GNodeArray array2(vector);
        test_try_success();
        for (int i = 0; i < 5; ++i) {
            test_value(array2[i], vector[i]);
        }
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test std::vector constructor
    test_try("std::vector constructor");
    try {
        double array[] = {0.0, 1.0, 2.0, 3.0, 5.0};
        std::vector<double> vector;
        for (int i = 0; i < 5; ++i) {
            vector.push_back(array[i]);
        }
        GNodeArray array2(vector);
        test_try_success();
        for (int i = 0; i < 5; ++i) {
            test_value(array2[i], vector[i]);
        }
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test linear interpolation
    double array_lin[] = {-1.0, 0.0, 1.0};
    test_node_array_interpolation(3, array_lin);

    // Test non-linear interpolation
    double array_nonlin[] = {-1.3, 0.0, 1.7};
    test_node_array_interpolation(3, array_nonlin);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GNodeArray class interpolation
 *
 * @param[in] num Number of nodes.
 * @param[in] nodes Nodes.
 *
 * Test the GNodeArray class interpolation method by comparing the
 * interpolation results for a linear function to the expected result.
 ***************************************************************************/
void TestGSupport::test_node_array_interpolation(const int&    num,
                                                 const double* nodes)
{
    // Setup function parameters
    const double slope  = 3.1;
    const double offset = 7.9;

    // Initialise node array
    GNodeArray array(num, nodes);

    // Setup node values
    std::vector<double> values;
    for (int i = 0; i < 3; ++i) {
        values.push_back(nodes[i] * slope + offset);
    }

    // Test values
    for (double value = -2.0; value <= +2.0; value += 0.2) {
        double expected = value * slope + offset;
        double result   = array.interpolate(value, values);
        test_value(result, expected);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GUrlFile class
 *
 * Test the GUrlFile class.
 ***************************************************************************/
void TestGSupport::test_url_file(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GUrlFile url;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test open constructor
    test_try("Open constructor");
    try {
        GUrlFile url("test_url.dat", "w");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test file writing
    GUrlFile url("test_url.dat", "w");
    test_value(url.write("abcd", 4), 4);
    url.putchar('e');
    url.printf("fghi%s%d%3.1fxyz", "jklm", 41, 9.9);
    url.close();

    // Test file reading using read() method
    char buffer[100];
    url.open("test_url.dat", "r");
    test_value(url.read(buffer, 99), 21);
    std::string result = std::string(buffer, 21);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\n");
    url.close();

    // Test file reading using getchar() method
    result.clear();
    url.open("test_url.dat", "r");
    int character = 0;
    do {
        character = url.getchar();
        if (character != EOF) {
            char c = (char)character;
            result.append(1, c);
        }
    } while (character != EOF);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\n");
    url.close();

    // Test file reading using scanf() method
    result.clear();
    url.open("test_url.dat", "r");
    url.scanf("%99s", buffer);
    result = std::string(buffer, 21);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\n");
    url.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GUrlString class
 *
 * Test the GUrlString class.
 ***************************************************************************/
void TestGSupport::test_url_string(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GUrlString url;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test open constructor
    test_try("Open constructor");
    try {
        GUrlString url("My nice text string");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test string writing
    GUrlString url;
    test_value(url.write("abcd", 4), 4);
    url.putchar('e');
    url.printf("fghi%s%d%3.1fxyz", "jklm", 41, 9.9);

    // Test string reading using read() method
    char buffer[100];
    url.rewind();
    test_value(url.read(buffer, 99), 21);
    std::string result = std::string(buffer, 21);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\n");

    // Test string reading using getchar() method
    result.clear();
    url.rewind();
    int character = 0;
    do {
        character = url.getchar();
        if (character != EOF) {
            char c = (char)character;
            result.append(1, c);
        }
    } while (character != EOF);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\n");

    // Test string reading using scanf() method
    char buffer2[100];
    result.clear();
    url.rewind();
    url.scanf("%99s", buffer2);
    result = std::string(buffer2, 21);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\n");
    url.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Main test entry point
 ***************************************************************************/
int main(void)
{
    // Allocate test suite container
    GTestSuites testsuites("Support module");

    // Initially assume that we pass all tests
    bool success = true;

    // Create and append test suite
    TestGSupport test;
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GSupport.xml");

    // Return success status
    return (success ? 0 : 1);
}
