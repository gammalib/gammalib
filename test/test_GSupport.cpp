/***************************************************************************
 *                 test_GSupport.cpp - test support module                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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
#include <cstdlib>     // getenv
#include <vector>
#include "GTools.hpp"
#include "test_GSupport.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir  = std::getenv("TEST_DATA");
const std::string csv_file = datadir + "/csv.dat";


/***********************************************************************//**
 * @brief Set parameters and tests
 ***************************************************************************/
void TestGSupport::set(void){

    // Set test name
    name("GSupport");

    // Append tests
    append(static_cast<pfunction>(&TestGSupport::test_tools),
           "Test GTools");
    append(static_cast<pfunction>(&TestGSupport::test_expand_env),
           "Test Environment variable");
    append(static_cast<pfunction>(&TestGSupport::test_node_array),
           "Test GNodeArray");
    append(static_cast<pfunction>(&TestGSupport::test_bilinear),
           "Test GBilinear");
    append(static_cast<pfunction>(&TestGSupport::test_url_file),
           "Test GUrlFile");
    append(static_cast<pfunction>(&TestGSupport::test_url_string),
           "Test GUrlString");
    append(static_cast<pfunction>(&TestGSupport::test_filename),
           "Test GFilename");
    append(static_cast<pfunction>(&TestGSupport::test_csv),
           "Test GCsv");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGSupport* TestGSupport::clone(void) const
{
    // Clone test suite
    return new TestGSupport(*this);
}


/***********************************************************************//**
 * @brief Test GTools
 *
 * Test the GTools functions
 ***************************************************************************/
void TestGSupport::test_tools(void)
{
    // Strip whitespace
    test_value(gammalib::strip_whitespace("  World  "), "World");
    test_value(gammalib::strip_whitespace("World  "), "World");
    test_value(gammalib::strip_whitespace("  World"), "World");

    // Strip characters
    test_value(gammalib::strip_chars("xxWorldyy", "xy"), "World");
    test_value(gammalib::strip_chars("Worldxy", "xy"), "World");
    test_value(gammalib::strip_chars("xxWorld", "x"), "World");

    // Right strip characters
    test_value(gammalib::rstrip_chars("xxWorldyy", "xy"), "xxWorld");
    test_value(gammalib::rstrip_chars("Worldxy", "xy"), "World");
    test_value(gammalib::rstrip_chars("xxWorld", "x"), "xxWorld");
    test_value(gammalib::rstrip_chars("0.0012300", "0"), "0.00123");

    // Test noun number conversion
    test_value(gammalib::number("World", 0), "Worlds");
    test_value(gammalib::number("World", 1), "World");
    test_value(gammalib::number("World", 2), "Worlds");

    // Test XML to string conversion
    std::string s_in  = "Hallo World, you \"are\" my 'nice' <planet> & place";
    std::string s_ref = "Hallo World, you &quot;are&quot; my &apos;nice&apos;"
                        " &lt;planet&gt; &amp; place";
    std::string s_out;
    s_out = gammalib::str2xml(s_in);
    test_value(s_out, s_ref);
    s_out = gammalib::xml2str(s_out);
    test_value(s_out, s_in);

    // Test power law flux computations
    test_value(gammalib::plaw_energy_flux(2.0, 3.0, 2.5, -1.0), 2.5);
    test_value(gammalib::plaw_energy_flux(2.0, 3.0, 2.5, -2.5), 2.56453824468);
    test_value(gammalib::plaw_photon_flux(2.0, 3.0, 2.5, -1.0), 1.01366277027);
    test_value(gammalib::plaw_photon_flux(2.0, 3.0, 2.5, -2.5), 1.06136118604);

    // Create vector of strings
    std::vector<std::string> strings;
    strings.push_back("This");
    strings.push_back("is");
    strings.push_back("a");
    strings.push_back("vector");
    test_assert(gammalib::contains(strings, "a"), "\"a\" in vector");
    test_assert(gammalib::contains(strings, "vector"), "\"vector\" in vector");
    test_assert(!gammalib::contains(strings, "b"), "\"b\" not in vector");
    test_assert(!gammalib::contains(strings, "isavector"), "\"isavector\" not in vector");

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
    test_assert(s_out == s_ref,"no environment variable",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // Get home environment variable
    std::string home(std::getenv("HOME"));

    // $ENV{HOME} environment variable within string
    s_in  = "My $ENV{HOME} is my castle.";
    s_ref = "My "+home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$ENV{HOME} environment variable within string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // ${HOME} environment variable within string
    s_in  = "My ${HOME} is my castle.";
    s_ref = "My "+home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"${HOME} environment variable within string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // $HOME/path environment variable within string
    s_in  = "My $HOME/path is my castle.";
    s_ref = "My "+home+"/path is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$HOME/path environment variable within string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // $HOME environment variable at end of string
    s_in  = "My $HOME";
    s_ref = "My "+home;
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$HOME environment variable at end of string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // Environment variable within single quotes
    s_in  = "My '$(HOME)' is my castle.";
    s_ref = s_in;
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Environment variable within single quotes",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // Environment variable within double quotes
    s_in  = "My \"$(HOME)\" is my castle.";
    s_ref = "My \""+home+"\" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Environment variable within double quotes",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // Non existing environment variable within string
    s_in  = "My $(HOMEIXYZ) is my castle.";
    s_ref = "My $(HOMEIXYZ) is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Non existing environment variable within string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // Empty environment variable within string
    s_in  = "My $() is my castle.";
    s_ref = "My $() is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Empty environment variable within string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // Two successive environment variables within string
    s_in  = "My $(HOME)$ENV{HOME} is my castle.";
    s_ref = "My "+home+home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"Two successive environment variables within string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // $(HOME) at beginning of string
    s_in  = "$(HOME) is my castle.";
    s_ref = home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$(HOME) at beginning of string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // $ENV(HOME) at beginning of string
    s_in  = "$ENV(HOME) is my castle.";
    s_ref = home+" is my castle.";
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$ENV(HOME) at beginning of string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // $(HOME) at end of string
    s_in  = "My castle is $(HOME)";
    s_ref = "My castle is "+home;
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$(HOME) at end of string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // $ENV{HOME}${HOME}$ENV(HOME)$(HOME) only string
    s_in  = "$ENV{HOME}${HOME}$ENV(HOME)$(HOME)";
    s_ref = home+home+home+home;
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"$ENV{HOME}${HOME}$ENV(HOME)$(HOME) only string",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // Test ~
    s_in  = "My castle is ~";
    s_ref = gammalib::expand_env("My castle is $(HOME)");
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"String with ~",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

    // Test ~+
    s_in  = "My castle is ~+";
    s_ref = gammalib::expand_env("My castle is $(PWD)");
    s_out = gammalib::expand_env(s_in);
    test_assert(s_out == s_ref,"String with ~+",
                "Unexpected string \""+s_out+"\" (expected \""+s_ref+"\")");

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

    // Test container manipulation
    test_try("Container manipulation");
    try {
        GNodeArray array;
        test_assert(array.is_empty(), "Node array should be empty");
        array.append(0.0);
        test_assert(!array.is_empty(), "Node array should not be empty");
        test_value(array.size(), 1);
        array.append(1.0);
        array.append(2.0);
        array.append(3.0);
        array.append(5.0);
        test_value(array.size(), 5);
        array.insert(4, 4.0);
        test_value(array.size(), 6);
        for (int i = 0; i < array.size(); ++i) {
            test_value(array[i], double(i));   // 0, 1, 2, 3, 4, 5
        }
        array.remove(1);
        array.remove(2);
        array.remove(3);
        test_value(array.size(), 3);
        for (int i = 0; i < array.size(); ++i) {
            test_value(array[i], double(2*i)); // 0, 2, 4
        }
        GNodeArray array2;
        array2.append(6.0);
        array2.append(8.0);
        test_value(array2.size(), 2);
        array.extend(array2);
        test_value(array.size(), 5);
        for (int i = 0; i < array.size(); ++i) {
            test_value(array[i], double(2*i)); // 0, 2, 4, 6, 8
        }
        test_try_success();
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

    // Test node array change using access operator[]
    GNodeArray nodes;
    nodes.append(1.0);
    nodes.append(2.0);
    nodes.append(3.0);
    nodes.set_value(2.5);
    test_value(nodes.inx_left(), 1, "Expected node 1");
    test_value(nodes.inx_right(), 2, "Expected node 2");
    test_value(nodes.wgt_left(), 0.5, 1.0e-6, "Expected weight 0.5");
    test_value(nodes.wgt_right(), 0.5, 1.0e-6, "Expected weight 0.5");
    nodes[1] = 3.0;
    nodes[2] = 4.0;
    nodes.set_value(2.5);
    test_value(nodes.inx_left(), 0, "Expected node 0");
    test_value(nodes.inx_right(), 1, "Expected node 1");
    test_value(nodes.wgt_left(), 0.25, 1.0e-6, "Expected weight 0.25");
    test_value(nodes.wgt_right(), 0.75, 1.0e-6, "Expected weight 0.75");

    // Remove test file
    GFilename filename("test_nodes.fits");
    filename.remove();

    // Check saving
    nodes.save("test_nodes.fits");
    GNodeArray load1("test_nodes.fits");
    test_value(load1.size(), 3, "GNodeArray should have 3 nodes.");
    test_value(load1[0], 1.0, 1.0e-7, "First node should be 1.");
    test_value(load1[1], 3.0, 1.0e-7, "First node should be 3.");
    test_value(load1[2], 4.0, 1.0e-7, "First node should be 4.");

    // Check saving in a different extnsion
    nodes.append(10.0);
    nodes.append(11.0);
    nodes.save("test_nodes.fits[NODE ARRAY]", true);
    GNodeArray load2("test_nodes.fits[NODE ARRAY]");
    test_value(load2.size(), 5, "GNodeArray should have 5 nodes.");
    test_value(load2[0], 1.0, 1.0e-7, "First node should be 1.");
    test_value(load2[1], 3.0, 1.0e-7, "First node should be 3.");
    test_value(load2[2], 4.0, 1.0e-7, "First node should be 4.");
    test_value(load2[3], 10.0, 1.0e-7, "First node should be 10.");
    test_value(load2[4], 11.0, 1.0e-7, "First node should be 11.");

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
 * @brief Test GBilinear class
 *
 * Test the GBilinear class.
 ***************************************************************************/
void TestGSupport::test_bilinear(void)
{
    // Test constructor
    GBilinear interpolator;

    // Test value setting and retrieving
    interpolator.index1()  = 1;
    interpolator.index2()  = 2;
    interpolator.index3()  = 3;
    interpolator.index4()  = 4;
    interpolator.weight1() = 1.0;
    interpolator.weight2() = 2.0;
    interpolator.weight3() = 3.0;
    interpolator.weight4() = 4.0;
    test_value(interpolator.index1(), 1);
    test_value(interpolator.index2(), 2);
    test_value(interpolator.index3(), 3);
    test_value(interpolator.index4(), 4);
    test_value(interpolator.weight1(), 1.0);
    test_value(interpolator.weight2(), 2.0);
    test_value(interpolator.weight3(), 3.0);
    test_value(interpolator.weight4(), 4.0);

    // Test copy constructor
    GBilinear interpolator2(interpolator);
    test_value(interpolator2.index1(), 1);
    test_value(interpolator2.index2(), 2);
    test_value(interpolator2.index3(), 3);
    test_value(interpolator2.index4(), 4);
    test_value(interpolator2.weight1(), 1.0);
    test_value(interpolator2.weight2(), 2.0);
    test_value(interpolator2.weight3(), 3.0);
    test_value(interpolator2.weight4(), 4.0);

    // Test interpolator
    double array[5]  = {1.0, 2.0, 3.0, 4.0, 5.0};
    double reference = 1.0 * array[1] +
                       2.0 * array[2] +
                       3.0 * array[3] +
                       4.0 * array[4];
    test_value(interpolator(array), reference);
    test_value(interpolator2(array), reference);

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
    url.put_char('e');
    url.printf("fghi%s%d%3.1fxyz", "jklm", 41, 9.9);
    url.close();

    // Test file reading using read() method
    char buffer[100];
    url.open("test_url.dat", "r");
    test_value(url.read(buffer, 99), 21);
    std::string result = std::string(buffer, 21);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\"");
    url.close();

    // Test file reading using get_char() method
    result.clear();
    url.open("test_url.dat", "r");
    int character = 0;
    do {
        character = url.get_char();
        if (character != EOF) {
            char c = (char)character;
            result.append(1, c);
        }
    } while (character != EOF);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\"");
    url.close();

    // Test file reading using scanf() method
    result.clear();
    url.open("test_url.dat", "r");
    url.scanf("%99s", buffer);
    result = std::string(buffer, 21);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\"");
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
    url.put_char('e');
    url.printf("fghi%s%d%3.1fxyz", "jklm", 41, 9.9);

    // Test string reading using read() method
    char buffer[100];
    url.rewind();
    test_value(url.read(buffer, 99), 21);
    std::string result = std::string(buffer, 21);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\"");

    // Test string reading using get_char() method
    result.clear();
    url.rewind();
    int character = 0;
    do {
        character = url.get_char();
        if (character != EOF) {
            char c = (char)character;
            result.append(1, c);
        }
    } while (character != EOF);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\"");

    // Test string reading using scanf() method
    char buffer2[100];
    result.clear();
    url.rewind();
    url.scanf("%99s", buffer2);
    result = std::string(buffer2, 21);
    test_assert(result.compare("abcdefghijklm419.9xyz") == 0,
                "Expected \"abcdefghijklm419.9xyz\" in file, found \""+
                result+"\"");
    url.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GFilename class
 *
 * Test the GFilename class.
 ***************************************************************************/
void TestGSupport::test_filename(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GFilename filename;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Set filename without extension
    GFilename filename = "myfile.fits";
    test_value(filename.url(), "myfile.fits");
    test_assert(!filename.has_extname(),
                "Expected that extension name is not set.");

    // Set filename with extension name
    filename = "myfile.fits[EVENTS]";
    test_value(filename.url(), "myfile.fits");
    test_value(filename.extname(), "EVENTS");
    test_assert(filename.has_extname(),
                "Expected that extension name is set.");

    // Set filename with extension name and version
    filename = "myfile.fits[EVENTS,2]";
    test_value(filename.url(), "myfile.fits");
    test_value(filename.extname(), "EVENTS");
    test_assert(filename.has_extname(),
                "Expected that extension name is set.");
    test_value(filename.extver(), 2);
    test_assert(filename.has_extver(),
                "Expected that extension version is set.");

    // Set filename with lower case extension name
    filename = "myfile.fits[events]";
    test_value(filename.extname(), "EVENTS");

    // Set filename with number in extension name
    filename = "myfile.fits[EVENTS_2D]";
    test_value(filename.extname(), "EVENTS_2D");

    // Set filename with extension number and version
    filename = "myfile.fits[1,2]";
    test_value(filename.url(), "myfile.fits");
    test_value(filename.extno(), 1);
    test_assert(filename.has_extno(),
                "Expected that extension number is set.");
    test_value(filename.extver(), 2);
    test_assert(filename.has_extver(),
                "Expected that extension version is set.");

    // Set filename with extension name and expression
    filename = "myfile.fits[EVENTS][ENERGY>0.1]";
    test_value(filename.url(), "myfile.fits");
    test_value(filename.extname(), "EVENTS");
    test_value(filename.expression(), "ENERGY>0.1");

    // Test missing closing symbol
    test_try("Missing ] symbol");
    try {
        GFilename filename("myfile.fits[");
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }

    // Test character after closing symbol
    test_try("Character after ] symbol");
    try {
        GFilename filename("myfile.fits[EVENTS]a");
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }

    // Test character after closing symbol
    test_try("Character after ] symbol after expression");
    try {
        GFilename filename("myfile.fits[EVENTS][ENERGY>0.1]a");
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }

    // Test empty extension name
    test_try("Empty extension name");
    try {
        GFilename filename("myfile.fits[]");
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }

    // Test invalid extension number
    test_try("Invalid extension number");
    try {
        GFilename filename("myfile.fits[-1]");
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }

    // Test missing extension version
    test_try("Missing extension version");
    try {
        GFilename filename("myfile.fits[EVENTS,]");
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }

    // Test invalid extension version
    test_try("Invalid extension version");
    try {
        GFilename filename("myfile.fits[EVENTS,-1]");
        test_try_failure();
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }

    // Test size operators
    filename = "myfile.fits";
    test_assert(!filename.is_empty(), "Non empty file name expected.");
    test_value(filename.length(), 11);

    // Test default extension name
    filename = "myfile.fits";
    test_value(filename.extname("EVENTS"), "EVENTS");
    filename = "myfile.fits[LUCKY]";
    test_value(filename.extname("EVENTS"), "LUCKY");

    // Test default extension number
    filename = "myfile.fits";
    test_value(filename.extno(1), 1);
    filename = "myfile.fits[3]";
    test_value(filename.extno(1), 3);

    // Test default extension version
    filename = "myfile.fits";
    test_value(filename.extver(1), 1);
    filename = "myfile.fits[3,3]";
    test_value(filename.extver(1), 3);

    // Test type method
    filename = "myfile.fits";
    test_value(filename.type(), "fits");
    filename = "myfile.fit";
    test_value(filename.type(), "fits");
    filename = "myfile.fits.gz";
    test_value(filename.type(), "fits");
    filename = "myfile.fits.Z";
    test_value(filename.type(), "fits");
    filename = "myfile.fits.z";
    test_value(filename.type(), "fits");
    filename = "myfile.fits.zip";
    test_value(filename.type(), "fits");
    filename = "myfile.csv";
    test_value(filename.type(), "");
    filename = "myfile.csv.gz";
    test_value(filename.type(), "");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GCsv class
 *
 * Test the GCsv class.
 ***************************************************************************/
void TestGSupport::test_csv(void)
{
    // Test void constructor
    GCsv csv1;
    test_value(csv1.size(), 0);
    test_value(csv1.ncols(), 0);
    test_value(csv1.nrows(), 0);

    // Test rows and columns constructor
    GCsv csv2(3,4);
    test_value(csv2.size(), 12);
    test_value(csv2.ncols(), 4);
    test_value(csv2.nrows(), 3);

    // Test filename constructor
    GCsv csv3(csv_file, ",");
    test_value(csv3.size(), 12);
    test_value(csv3.ncols(), 3);
    test_value(csv3.nrows(), 4);

    // Test copy constructor
    GCsv csv4(csv3);
    test_value(csv3.size(), 12);
    test_value(csv3.ncols(), 3);
    test_value(csv3.nrows(), 4);

    // Test clear method
    csv4.clear();
    test_value(csv4.size(), 0);
    test_value(csv4.ncols(), 0);
    test_value(csv4.nrows(), 0);

    // Test clone method
    GCsv* csv5 = csv3.clone();
    test_value(csv5->size(), 12);
    test_value(csv5->ncols(), 3);
    test_value(csv5->nrows(), 4);

    // Test append method
    std::vector<std::string> row1;
    std::vector<std::string> row2;
    row1.push_back("ra");
    row1.push_back("dec");
    row1.push_back("durations");
    row2.push_back("10.0");
    row2.push_back("-10.0");
    row2.push_back("1000.0");
    csv4.append(row1);
    csv4.append(row2);
    test_value(csv4.size(), 6);
    test_value(csv4.ncols(), 3);
    test_value(csv4.nrows(), 2);

    // Test operator access method
    test_value(csv4(0,0), "ra");
    test_value(csv4(0,1), "dec");
    test_value(csv4(0,2), "durations");
    test_value(csv4(1,0), "10.0");
    test_value(csv4(1,1), "-10.0");
    test_value(csv4(1,2), "1000.0");

    // Test access methods
    test_value(csv4.string(1,0), "10.0");
    test_value(csv4.real(1,0), 10.0);
    test_value(csv4.integer(1,0), 10);

    // Test set methods
    csv4.string(1,0,"11.0");
    test_value(csv4.string(1,0), "11.0");
    csv4.real(1,0, 12.0);
    test_value(csv4.real(1,0), 12.0);
    csv4.integer(1,0, 13);
    test_value(csv4.integer(1,0), 13);

    // Test save method (make sure that former result is overwritten)
    csv4.save("test_csv.dat", ",", true);
 
    // Test load method
    GCsv csv6;
    csv6.load("test_csv.dat", ",");
    test_value(csv6.size(), 6);
    test_value(csv6.ncols(), 3);
    test_value(csv6.nrows(), 2);
    test_value(csv6(0,0), "ra");
    test_value(csv6(0,1), "dec");
    test_value(csv6(0,2), "durations");
    test_value(csv6(1,0), "13");
    test_value(csv6(1,1), "-10.0");
    test_value(csv6(1,2), "1000.0");

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

    // Append test to the container
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GSupport.xml");

    // Return success status
    return (success ? 0 : 1);
}
