/***************************************************************************
 *              GTestSuites.cpp - Test suite container class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Jean-Baptiste Cayrou                        *
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
 * @file GTestSuites.cpp
 * @brief Test suite container class implementation
 * @author Jean-Baptiste Cayrou
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <ctime>
#include "GTestSuites.hpp"
#include "GTestCase.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                          "GTestSuites::at(int&)"
#define G_SET                           "GTestSuites::set(int&, GTestSuite&)"
#define G_INSERT         "GTestSuite* GTestSuites::insert(int&, GTestSuite&)"
#define G_REMOVE                                  "GTestSuites::remove(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GTestSuites::GTestSuites(void)
{
    // Initialise members
    init_members();

    //Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] suites Test suite container.
 ***************************************************************************/
GTestSuites::GTestSuites(const GTestSuites& suites)
{
    // Initialise members
    init_members();
    
    // Copy members
    copy_members(suites);
    
    //Return
    return;
}


/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Name of the test suites.
 ***************************************************************************/
GTestSuites::GTestSuites(const std::string& name)
{
    // Initialise members
    init_members();
    
    // Set name
    m_name = name;

    //Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GTestSuites::~GTestSuites()
{
    //Free members
    free_members();

    //Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] suites Test suite container.
 * @return Test suite container.
 ***************************************************************************/
GTestSuites& GTestSuites::operator=(const GTestSuites& suites)
{
    // Execute only if object is not identical
    if (this != &suites) {
        
        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(suites);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear test suites
 ***************************************************************************/
void GTestSuites::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suites
 *
 * @return Pointer to deep copy of test suites.
 ***************************************************************************/
GTestSuites* GTestSuites::clone(void) const
{
    // Clone object
    return new GTestSuites(*this);
}


/***********************************************************************//**
 * @brief Returns pointer to test suite
 *
 * @param[in] index Test suite index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Test suite index is out of range.
 ***************************************************************************/
GTestSuite* GTestSuites::at(const int& index)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Test suite index", index, size());
    }
    #endif

    // Return test suite
    return (m_testsuites[index]);
}


/***********************************************************************//**
 * @brief Returns pointer to test suite (const version)
 *
 * @param[in] index Test Suite index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Test suite index is out of range.
 ***************************************************************************/
const GTestSuite* GTestSuites::at(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Test suite index", index, size());
    }
    #endif

    // Return test suite
    return (m_testsuites[index]);
}


/***********************************************************************//**
 * @brief Set test suite in container
 *
 * @param[in] index Test suite index [0,...,size()-1].
 * @param[in] suite Test suite.
 * @return Pointer to deep copy of test suite.
 *
 * @exception GException::out_of_range
 *            Test suite index is out of range.
 *
 * Set test suite in the container. A deep copy of the test suite will be
 * made.
 ***************************************************************************/
GTestSuite* GTestSuites::set(const int& index, const GTestSuite& suite)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET, "Test suite index", index, size());
    }
    #endif

    // Delete any existing test suite
    if (m_testsuites[index] != NULL) delete m_testsuites[index];

    // Assign new test suite by cloning
    m_testsuites[index] = suite.clone();

    // Return pointer to test suite
    return m_testsuites[index];
}


/***********************************************************************//**
 * @brief Append test suite to container
 *
 * @param[in] suite Test suite.
 * @return Pointer to deep copy of test suite.
 *
 * Appends test suite to the container by making a deep copy of the test
 * suite and storing its pointer.
 ***************************************************************************/
GTestSuite* GTestSuites::append(const GTestSuite& suite)
{
    // Create deep copy of test suite
    GTestSuite* ptr = suite.clone();

    // Append deep copy to container
    m_testsuites.push_back(ptr);

    // Return pointer to test suite
    return ptr;
}


/***********************************************************************//**
 * @brief Insert test suite into container
 *
 * @param[in] index Test suite index [0,...,size()-1].
 * @param[in] suite Test suite.
 * @return Pointer to deep copy of test suite.
 *
 * @exception GException::out_of_range
 *            Test suite index is out of range.
 *
 * Inserts a test @p suite into the container before the test suite with the
 * specified @p index.
 ***************************************************************************/
GTestSuite* GTestSuites::insert(const int& index, const GTestSuite& suite)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (isempty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "Test suite index", index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "Test suite index", index, size());
        }
    }
    #endif

    // Create deep copy of test suite
    GTestSuite* ptr = suite.clone();

    // Append deep copy to container
    m_testsuites.insert(m_testsuites.begin()+index, ptr);

    // Return pointer to test suite
    return ptr;
}


/***********************************************************************//**
 * @brief Remove test suite from container
 *
 * @param[in] index Test suite index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Test suite index is out of range.
 *
 * Remove test suite of specified @p index from container.
 ***************************************************************************/
void GTestSuites::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Test suite index", index, size());
    }
    #endif

    // Erase test suite from container
    m_testsuites.erase(m_testsuites.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append test suite container
 *
 * @param[in] suites Test suite container.
 *
 * Append test suite container to the container.
 ***************************************************************************/
void GTestSuites::extend(const GTestSuites& suites)
{
    // Do nothing if test suite container is empty
    if (!suites.isempty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = suites.size();

        // Reserve enough space
        reserve(size() + num);

        // Append deep copies of test suites to container 
        for (int i = 0; i < num; ++i) {
            m_testsuites.push_back(suites[i]->clone());
        }

    } // endif: test suite container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return the total number of errors in all test suites
 ***************************************************************************/
int GTestSuites::errors(void) const
{
    // Initialise number of errors
    int errors = 0;
    
    // Add up the errors from all test suites
    for (int i = 0; i < m_testsuites.size(); ++i) {
        errors += m_testsuites[i]->errors(); 
    }

    // Return errors
    return errors;
}


/***********************************************************************//**
 * @brief Return the total number of failures in all test suites
 ***************************************************************************/
int GTestSuites::failures(void) const
{
    // Initialise number of failures
    int failures=0;
    
    // Add up the failures from all test suites
    for (int i = 0; i < m_testsuites.size(); ++i) {
        failures += m_testsuites[i]->failures(); 
    }

    // Return failures
    return failures;
}


/***********************************************************************//**
 * @brief Return the total number of tests they are in all test suites
 ***************************************************************************/
int GTestSuites::tests(void) const
{
    // Initialise number of tests
    int tests = 0;
    
    // Add up the number of tests from all test suites
    for (int i = 0; i < m_testsuites.size(); ++i) {
        tests += m_testsuites[i]->size(); 
    }

    // Return number of tests
    return tests; 
}


/***********************************************************************//**
 * @brief Run all tests
 *
 * Runs all tests that are in the test suite container. The method returns
 * true if all test suites succeeded. If one test suite failed, false is
 * returned.
 ***************************************************************************/
bool GTestSuites::run(void)
{
    // Set header
    std::string text  = "* " + name() + " *";
    std::string frame = gammalib::fill("*", text.length());
    
    // Log name of test suites
    std::cout << std::endl;
    std::cout << frame << std::endl;
    std::cout << text  << std::endl;
    std::cout << frame << std::endl;

    // Initialise success flag
    bool success = true;
    
    // Run all test suites
    for (int i = 0; i < m_testsuites.size(); ++i) {
        if (!m_testsuites[i]->run()) {
            success = false;
        }
    }

    // Return success flag
    return success;
}


/***********************************************************************//**
 * @brief Save test report into XML file.
 *
 * @param[in] filename Name of XML file.
 *
 * Saves the test results in a JUnit compliant format into an XML file.
 ***************************************************************************/
void GTestSuites::save(const std::string& filename) const
{
    // Declare empty XML document
    GXml xml;

    // Write observations into XML file
    write(xml);

    // Save XML document
    xml.save(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print test suites information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing test suites information.
 ***************************************************************************/
std::string GTestSuites::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GTestSuites ===");

        // Append information
        result.append("\n"+gammalib::parformat("Name")+m_name);
        result.append("\n"+gammalib::parformat("Number of test suites"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of tests"));
        result.append(gammalib::str(tests()));
        result.append("\n"+gammalib::parformat("Number of errors"));
        result.append(gammalib::str(errors()));
        result.append("\n"+gammalib::parformat("Number of failures"));
        result.append(gammalib::str(failures()));

        // Append test suites
        for (int i = 0; i < size(); ++i) {
            result.append("\n");
            result.append((*this)[i]->print(chatter));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GTestSuites::init_members(void)
{
    // Initialise members
    m_name      = "Unnamed Test Suites";
    m_testsuites.clear();
    m_timestamp = time(NULL);
    m_log.clear();

    // Set logger parameters
    m_log.max_size(1);
    m_log.cout(true);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] suites Test suites container.
 *
 * This method just clone the container not the test suite.
 ***************************************************************************/
void GTestSuites::copy_members(const GTestSuites& suites)
{
    // Copy members
    m_name       = suites.m_name;
    m_testsuites = suites.m_testsuites;
    m_timestamp  = suites.m_timestamp;
    m_log        = suites.m_log;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTestSuites::free_members(void)
{
    // Close logger
    m_log.close();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Test Suites into XML document
 *
 * @param[in] xml XML document.
 *
 * Write test suites into the XML document.
 ***************************************************************************/
void GTestSuites::write(GXml& xml) const
{
    // Append XML element for Test Suites
    GXmlElement* testsuites = xml.append("testsuites");

    // Set attributes
    testsuites->attribute("name", "GammaLib");
    /*
    testsuites->attribute("test",str(tests()));
    testsuites->attribute("errors",str(errors()));
    testsuites->attribute("failures",str(failures()));
    testsuites->attribute("time","0"); // not used
    testsuites->attribute("timestamp",str(timestamp()));
    */

    // Loop over all test suites in the container
    for (int i = 0; i < m_testsuites.size(); ++i) {

        // Get a pointer on the current test suite
        GTestSuite* testsuite = m_testsuites[i];

        // Append XML element for this test suite
        GXmlElement* element_testsuite = testsuites->append("testsuite");

        // Set attributes
        element_testsuite->attribute("disabled","");  // not used
        element_testsuite->attribute("errors",gammalib::str(testsuite->errors()));
        element_testsuite->attribute("failures",gammalib::str(testsuite->failures()));
        element_testsuite->attribute("hostname","");  // not used
        element_testsuite->attribute("id",gammalib::str(i));
        element_testsuite->attribute("name",testsuite->name()); 
        element_testsuite->attribute("package","");  // not used
        element_testsuite->attribute("skipped","");  // not used
        element_testsuite->attribute("tests",gammalib::str(testsuite->size()));
        element_testsuite->attribute("time",gammalib::str(testsuite->duration()));
        element_testsuite->attribute("timestamp",gammalib::str(testsuite->timestamp()));

        // Loop over all test cases in the test suite
        for (int j = 0; j < testsuite->size(); ++j) {

            // Reference to the current test case
            GTestCase& testcase = (*testsuite)[j];

            // Append XML element for this test case
            GXmlElement* element_testcase = element_testsuite->append("testcase");

            // Set attributes
            element_testcase->attribute("assertions",""); // not used
            element_testcase->attribute("classname",name());
            element_testcase->attribute("name",testcase.name());
            element_testcase->attribute("status","");  // not used
            element_testcase->attribute("time",gammalib::str(testcase.duration()));

            // If a failure or error occured then append the message to the
            // XML element.
            if (!testcase.passed()) {

                // Append XML element for the test case problem
                GXmlElement* element_testcase_problem = element_testcase->append("error");

                // Set attributes
                element_testcase_problem->attribute("message",testcase.message());
                element_testcase_problem->attribute("type",testcase.type());

                // Set tag name dependent on the type of problem
                if (testcase.kind() == GTestCase::ERROR_TEST) {
                    element_testcase_problem->name("error");
                }
                else {
                    element_testcase_problem->name("failure");
                }

            } // endif: failure or error occured

        } // endfor: looped over all test cases

    } // endfor: looped over all test suites

    // Return
    return;
}
