/***************************************************************************
 *             GTestSuite.cpp - Abstract test suite base class             *
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
 * @file GTestSuite.cpp
 * @brief Abstract test suite base class implementation
 * @author Jean-Baptiste Cayrou
 */
 
/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo>
#include "GTestSuite.hpp"
#include "GTools.hpp"
#include "GLog.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#ifdef HAVE_OPENMP_DARWIN_KLUGE
#include <pthread.h>
pthread_attr_t gomp_thread_attr;
#endif
#endif

/* __ Method name definitions ____________________________________________ */
#define G_OP_ACCESS                            "GTestSuite::operator[](int&)"
#define G_TRY_SUCCESS                        "GTestSuite::test_try_success()"
#define G_TRY_FAILURE1           "GTestSuite::test_try_failure(std::string&,"\
                                                             " std::string&)"
#define G_TRY_FAILURE2        "GTestSuite::test_try_failure(std::exception&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GTestSuite::GTestSuite(void)
{
    // Initialise members
    init_members();
    
    //Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] suite Test Suite.
 ***************************************************************************/
GTestSuite::GTestSuite(const GTestSuite& suite)
{
    // Initialise members
    init_members();
    
    // Copy members
    copy_members(suite);
    
    //Return
    return;
}


/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Test suite name.
 ***************************************************************************/
GTestSuite::GTestSuite(const std::string& name)
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
GTestSuite::~GTestSuite(void)
{
    // Free members
    free_members();

    // Return
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
 * @param[in] suite Test suite.
 * @return Test suite.
 ***************************************************************************/
GTestSuite& GTestSuite::operator=(const GTestSuite& suite)
{
    // Execute only if object is not identical
    if (this != &suite) {
        
        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(suite);

    } // endif: object was not identical

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Returns reference to test case
 *
 * @param[in] index Test case index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Test case index is out of range.
 ***************************************************************************/
GTestCase& GTestSuite::operator[](const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OP_ACCESS, "Test case index", index, size());
    }
    #endif

    // Return reference
    return *(m_tests[index]);
}

/***********************************************************************//**
 * @brief Returns reference to test case
 *
 * @param[in] index Test case index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Test case index is out of range.
 ***************************************************************************/
const GTestCase& GTestSuite::operator[](const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OP_ACCESS, "Test case index", index, size());
    }
    #endif

    // Return reference
    return *(m_tests[index]);
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear test suite
 ***************************************************************************/
void GTestSuite::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append test functions to test suite
 *
 * @param[in] function Test function pointer.
 * @param[in] name Test name.
 *
 * This method adds test functions to the test suite. The test functions will
 * be executed when the run method is called.
 ***************************************************************************/
void GTestSuite::append(const pfunction function, const std::string& name)
{
    // Add test function pointer and name to suite
    m_functions.push_back(function);
    m_names.push_back(name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Run all tests in test suite
 *
 * @return True is all tests were successful, false otherwise.
 *
 * Executes all test functions that have been appended to the test suite.
 * For each test function a test case is added to the test suite.
 ***************************************************************************/
bool GTestSuite::run(void)
{
    // Setup the test functions. This is a pure virtual function that needs
    // to be implemented in the derived class. It sets the function
    // pointers and function names for all test functions that should be
    // executed.
    set();

    // Initialise success flag
    bool success = true;

    // Loop over all functions in suite
    for (m_index = 0; m_index < m_functions.size(); ++m_index) {

        // Continue only if function is valid
        if (m_functions[m_index] != NULL) {

            // Save the number of errors and failures before test
            // execution. We use this after the test to see if
            // any failures occured.
            int old_errors   = errors();
            int old_failures = failures();

            // Log the name of the test
            std::cout << m_names[m_index] << ": ";

            // Create a test of error type for function testing
            GTestCase* test = new GTestCase(GTestCase::ERROR_TEST, m_names[m_index]);

            // Set start time
            #ifdef _OPENMP
            double t_start = omp_get_wtime();
            #else
            clock_t t_start = clock();
            #endif

            // Execute test function
            try {
                (this->*(m_functions[m_index]))();
            }
            catch (std::exception& e) {

                // Signal that test did not succeed
                test->passed(false);

                // Set test message to exception message
                test->message(e.what());

                // Set type as class name
                test->type(typeid(e).name());

            }
            catch (...)
            {
                // For other exceptions
                test->passed(false);
            }

            // Compute elapsed time
            #ifdef _OPENMP
            double t_elapse = omp_get_wtime()-t_start;
            #else
            double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
            #endif

            // Set test duration
            test->duration(t_elapse);

            // Increment number of errors if the test did not pass
            if (!test->passed()) {
                m_errors++;
            }

            // Log the result (".","F" or, "E")
            std::cout << test->print();

            // Log if there are errors or failures
            if ((m_errors == old_errors && m_failures == old_failures)) {
                std::cout << " ok" << std::endl;
            }
            else {
                std::cout << " NOK" << std::endl;
                success = false;
            }

            // Add test case to test suite
            m_tests.push_back(test);

        } // endif: test case has a function pointer

    } // endwhile: there were tests on the stack

    // Reset index
    m_index = 0;

    // Return success flag
    return success;
}


/***********************************************************************//**
 * @brief Test an assert
 *
 * @param[in] assert Assert (true/false).
 * @param[in] name Test case name (defaults to "").
 * @param[in] message Test case name (defaults to "").
 *
 * Tests if a condition is true or false. This method adds a test case of
 * type "failure" to the test suite.
 *
 * Examples:
 *   test_assert(x>3, "Test if x > 3");
 *   test_assert(x>3 && x<10, "Test if  3 < x < 10 ");
 ***************************************************************************/
void GTestSuite::test_assert(const bool&        assert,
                             const std::string& name,
                             const std::string& message)
{
    // Create a test case of failure type
    GTestCase* testcase = new GTestCase(GTestCase::FAIL_TEST, format_name(name));

    // If assert is false then signal that the test is not passed and
    // increement the number of failures in this test suite 
    if (!assert) {
        testcase->passed(false);
        m_failures++;
    }

    // Set message
    testcase->message(message);

    // Log the result (".","F" or, "E")
    std::cout << testcase->print();

    // Add test case to test suite
    m_tests.push_back(testcase);

    // Return
    return;

}

/***********************************************************************//**
 * @brief Test a value
 *
 * @param[in] value Value to test.
 * @param[in] expected Expected value.
 * @param[in] name Test case name (defaults to "").
 * @param[in] message Test case message (defaults to "").
 *
 * Test if the value is comprised in the interval
 * [expected-eps, expected+eps].
 ***************************************************************************/
void GTestSuite::test_value(const int&         value,
                            const int&         expected,
                            const std::string& name,
                            const std::string& message)
{
    // Set test case name. If no name is specify then build the name from
    // the actual test parameters.
    std::string formated_name;
    if (name != "") {
        formated_name = format_name(name);
    }
    else {
        formated_name = format_name("Test if " + gammalib::str(value) + " is " + 
                                    gammalib::str(expected));
    }

    // Create a test case of failure type
    GTestCase* testcase = new GTestCase(GTestCase::FAIL_TEST, formated_name);

    // If value is not the expected one then signal test as failed and
    // increment the number of failures
    if (value != expected) {
        testcase->passed(false);
        m_failures++;
    }

    // If no message is specified then build message from test result
    std::string formated_message;
    if (message != "") {
        formated_message = message;
    }
    else {
        formated_message = "Value " + gammalib::str(value) + " equals not the " +
                           "expected value of " + gammalib::str(expected) + ".";
    }

    // Set message
    testcase->message(formated_message);

    // Log the result (".","F" or, "E")
    std::cout << testcase->print();

    // Add test case to test suite
    m_tests.push_back(testcase);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test a value
 *
 * @param[in] value Value to test.
 * @param[in] expected Expected value.
 * @param[in] eps Precision of the test (default 1e-30).
 * @param[in] name Test case name (defaults to "").
 * @param[in] message Test case message (defaults to "").
 *
 * Test if the value is comprised in the interval
 * [expected-eps, expected+eps].
 ***************************************************************************/
void GTestSuite::test_value(const double&      value,
                            const double&      expected,
                            const double&      eps,
                            const std::string& name,
                            const std::string& message)
{
    // Set test case name. If no name is specify then build the name from
    // the actual test parameters.
    std::string formated_name;
    if (name != "") {
        formated_name = format_name(name);
    }
    else {
        formated_name = format_name("Test if " + gammalib::str(value) + " is within " + 
                                    gammalib::str(expected) + " +/- " + gammalib::str(eps));
    }

    // Create a test case of failure type
    GTestCase* testcase = new GTestCase(GTestCase::FAIL_TEST, formated_name);

    // If value is not between in interval [expected-eps, expected+eps]
    // then signal test as failed and increment the number of failures
    if (value > expected + eps || value < expected - eps) {
        testcase->passed(false);
        m_failures++;
    }

    // If no message is specified then build message from test result
    std::string formated_message;
    if (message != "") {
        formated_message = message;
    }
    else {
        formated_message = "Value " + gammalib::str(value) + " not within " +
                           gammalib::str(expected) + " +/- " + gammalib::str(eps) +
                           " (value-expected = " + gammalib::str(value-expected) + ").";
    }

    // Set message
    testcase->message(formated_message);

    // Log the result (".","F" or, "E")
    std::cout << testcase->print();

    // Add test case to test suite
    m_tests.push_back(testcase);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test an try block
 *
 * @param[in] name Test case name (defaults to "").
 *
 * @see test_try_sucess() 
 * @see test_try_failure(const std::string& message,const std::string& type)
 * @see test_try_failure(const std::exception& e)
 *
 * Call before testing a try/catch block.
 *
 * Example: 
 *       test_try("Test a try block");
 *       try {
 *          ... //someting to test
 *          test_try_success();
 *       }
 *       catch(...) {
 *          test_try_failure();
 *       }
 ***************************************************************************/
void GTestSuite::test_try(const std::string& name)
{
    // Create a test case of error type
    GTestCase* testcase = new GTestCase(GTestCase::ERROR_TEST, format_name(name));

    // Add test case to try stack of test suite
    m_stack_try.push_back(testcase);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Notice when a try block succeeded
 *
 * @exception GException::test_nested_try_error
 *            Test case index is out of range.
 *
 * @see test_try(const std::string& name)
 * @see test_try_failure(const std::string& message, const std::string& type)
 * @see test_try_failure(const std::exception& e)
 *
 * Call this method at the last line of a try
 *
 * Example: 
 *       test_try("Test a try block");
 *       try {
 *          ... //someting to test
 *          test_try_success();
 *       }
 *       catch(...) {
 *          test_try_failure();
 *       }
 ***************************************************************************/
void GTestSuite::test_try_success(void)
{
    // If the stack is empty
    if (m_stack_try.empty()) {
        throw GException::test_nested_try_error(G_TRY_SUCCESS, 
              "Called \"G_TRY_SUCCESS\" without a previous call to test_try()");
    }

    // Add test case to test suite
    m_tests.push_back(m_stack_try.back());

    // Delete the test case from the try stack
    m_stack_try.pop_back();

    // Log the result (".","F" or, "E")
    std::cout << m_tests.back()->print();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Notice when a try block failed
 *
 * @param[in] message Message to explain why test failed (defaults to "").
 * @param[in] type Type of message (defaults to "").
 *
 * @exception GException::test_nested_try_error
 *            Test case index is out of range.
 *
 * @see test_try_sucess()
 * @see test_try(const std::string& name)
 * @see test_try_failure(const std::exception& e)
 *
 * Call this method in the catch block.
 *
 * Example: 
 *       test_try("Test a try block");
 *       try {
 *          ... //someting to test
 *          test_try_success();
 *       }
 *       catch(...) {
 *          test_try_failure();
 *       }
 ***************************************************************************/
void GTestSuite::test_try_failure(const std::string& message,
                                  const std::string& type)
{
    // If the stack is empty then create an eception test case
    if (m_stack_try.empty()) {
        GTestCase* testcase = new GTestCase(GTestCase::ERROR_TEST, "Exception test");
        m_stack_try.push_back(testcase);
    }

    // Signal that test is not ok
    m_stack_try.back()->passed(false);

    // Increment the number of errors
    m_errors++;

    // Set test type
    m_stack_try.back()->kind(GTestCase::ERROR_TEST);
    
    // Set message
    m_stack_try.back()->message(message);

    // Set type of message
    m_stack_try.back()->type(type);

    // Add test case to test suite
    m_tests.push_back(m_stack_try.back());
    
    // Delete the test case from the stack
    m_stack_try.pop_back();

    // Log the result ( ".","F" or, "E")
    std::cout << m_tests.back()->print();

    //Return
    return;
}


/***********************************************************************//**
 * @brief Notice when a try block failed
 *
 * @param[in] e Exception.
 *
 * @exception GException::test_nested_try_error
 *            Test case index is out of range.
 *
 * @see test_try_sucess()
 * @see test_try(const std::string& name)
 * @see test_try_failure(const std::string& message, const std::string& type)
 *
 * Call this method in a catch block.
 *
 * Example: 
 *       test_try("Test a try block");
 *       try {
 *          ... //someting to test
 *          test_try_success();
 *       }
 *       catch(exception& e) {
 *          test_try_failure(e);
 *       }
 ***************************************************************************/
void GTestSuite::test_try_failure(const std::exception& e)
{
    // Extract message of exception and class name
    test_try_failure(e.what(), typeid(e).name());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return a failure exception
 *
 * @param[in] message Message.
 *
 * @see test_try()
 * @see test_error(const std::string& message)
 *
 * It can be use in a try test
 *
 * Example: 
 *       test_try("Test a try block");
 *       try {
 *          throw exception_failure("a failure");
 *          test_try_success();
 *       }
 *       catch(exception& e) {
 *          test_try_failure(e);
 *       }
 ***************************************************************************/
GException::test_failure& GTestSuite::exception_failure(const std::string& message)
{
    // Return exception
    return *(new GException::test_failure(m_stack_try.back()->name(), message));
}


/***********************************************************************//**
 * @brief Return an error exception
 *
 * @param[in] message Message.
 *
 * @see test_try()
 * @see test_failure(const std::string& message)
 *
 * It can be use in a try test
 *
 * Example: 
 *       test_try("Test a try block");
 *       try {
 *          throw exception_error("an error");
 *          test_try_success();
 *       }
 *       catch(exception& e) {
 *          test_try_failure(e);
 *       }
 ***************************************************************************/
GException::test_error& GTestSuite::exception_error(const std::string& message)
{
    // Return exception
    return *(new GException::test_error(m_stack_try.back()->name(),message));
}


/***********************************************************************//**
 * @brief Return the number of successful tests
 ***************************************************************************/
int GTestSuite::success(void) const
{
    // Return successes
    return size()-(m_errors+m_failures);
}


/***********************************************************************//**
 * @brief Return the total duration of all tests
 *
 * This method sums up all test durations and returns the result.
 ***************************************************************************/
double GTestSuite::duration(void) const
{
    // Initialise duration
    double duration = 0.0;

    // Add up the durations of all tests
    for (int i = 0; i < m_tests.size(); ++i) {
        duration += m_tests[i]->duration();
    }

    // Return duration
    return duration;
}


/***********************************************************************//**
 * @brief Print test suite information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing test suite information.
 ***************************************************************************/
std::string GTestSuite::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GTestSuite ===");

        // Append information
        result.append("\n"+gammalib::parformat("Name")+m_name);
        result.append("\n"+gammalib::parformat("Number of functions"));
        result.append(gammalib::str(m_names.size()));
        result.append("\n"+gammalib::parformat("Number of executed tests"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of errors"));
        result.append(gammalib::str(errors()));
        result.append("\n"+gammalib::parformat("Number of failures"));
        result.append(gammalib::str(failures()));

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
void GTestSuite::init_members(void)
{
    // Initialise members
    m_name      = "Unnamed Test Suite";
    m_names.clear();
    m_functions.clear();
    m_tests.clear();
    m_stack_try.clear();
    m_index     = 0;
    m_failures  = 0;
    m_errors    = 0;
    m_log.clear();
    m_timestamp = time(NULL);

    // Set logger parameters
    cout(true);
    m_log.max_size(1);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] suite Test suite.
 *
 * This method just clone the container not the test case.
 ***************************************************************************/
void GTestSuite::copy_members(const GTestSuite& suite)
{
    // Copy members
    m_name       = suite.m_name;
    m_names      = suite.m_names;
    m_functions  = suite.m_functions;
    m_index      = suite.m_index;
    m_failures   = suite.m_failures;
    m_errors     = suite.m_errors;
    m_log        = suite.m_log;
    m_timestamp  = suite.m_timestamp;

    // Clone test cases
    for (int i = 0; i < suite.m_tests.size(); ++i) {
        m_tests[i] = suite.m_tests[i]->clone(); 
    }

    // Clone try stack
    for (int i = 0; i < suite.m_stack_try.size(); ++i) {
        m_stack_try[i] = suite.m_stack_try[i]->clone(); 
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTestSuite::free_members(void)
{
    // Delete test cases
    for (int i = 0; i < m_tests.size(); ++i) {
        delete m_tests[i];
        m_tests[i] = NULL;
    }

    // Delete try stack
    for (int i = 0; i < m_stack_try.size(); ++i) {
        delete m_stack_try[i];
        m_stack_try[i] = NULL;
    }

    // Close logger
    m_log.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Format Name
 *
 * Return a string with the format
 * "TestFunctionName:TestTryname1:TestTryName2: name"
 ***************************************************************************/
std::string GTestSuite::format_name(const std::string& name)
{
    // Initialise format name
    std::string format_name;

    // Set name of the try blocks
    if (!m_stack_try.empty()) {
        format_name = m_stack_try.back()->name();
    }
    else {
        // Set name of the test
        format_name = m_names[m_index];
    }

    // Append test suite name
    format_name += ": " + name;

    // Return format
    return format_name;
}
