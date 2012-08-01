/***************************************************************************
 *         GTestSuite.cpp  - Test Suite class for GammaLib                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jean-Baptiste Cayrou                             *
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
 * @brief Test Suite class implementation
 * @author Jean-Baptiste Cayrou
 */
 
/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTestSuite.hpp"
#include "GLog.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_ACCESS                         "GTestSuite::operator[](int&)"

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
 * @param[in] testsuite Test Suite.
 ***************************************************************************/
GTestSuite::GTestSuite(const GTestSuite& testsuite)
{
    // Initialise members
    init_members();
    
    // Copy members
    copy_members(testsuite);
    
    //Return
    return;
}

/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Name of the Test Suite
 ***************************************************************************/
GTestSuite::GTestSuite(const std::string& name)
{
    // Initialise members
    init_members();
    
    //Set name
    m_name=name;
    
    //Return
    return;
}

/***********************************************************************//**
 * @brief Destructor
 * Exit(1) if they are failures or errors for compatibility with make check
 ***************************************************************************/
GTestSuite::~GTestSuite(void)
{
    //For make check
    if(m_errors+m_failures>0){
        exit(1);
    }
    
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
 * @param[in] testsuite Test Suite.
 ***************************************************************************/
 GTestSuite&  GTestSuite::operator= (const GTestSuite& testsuite)
{
    // Execute only if object is not identical
    if (this != &testsuite) {
        
        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(testsuite);

    } // endif: object was not identical

    // Return
    return *this;
}

/***********************************************************************//**
 * @brief Returns reference to test case
 *
 * @param[in] index Test case index [0,...,tests()-1].
 *
 * @exception GException::out_of_range
 *            Test case index is out of range.
 ***************************************************************************/
GTestCase& GTestSuite::operator[](const int& index)
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= tests()) {
        throw GException::out_of_range(G_OP_ACCESS, index, 0, tests()-1);
    }

    // Return reference
    return *(m_tests[index]);
}

/***********************************************************************//**
 * @brief Returns reference to test case
 *
 * @param[in] index Test case index [0,...,tests()-1].
 *
 * @exception GException::out_of_range
 *            Test case index is out of range.
 ***************************************************************************/
const GTestCase& GTestSuite::operator[](const int& index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= tests()) {
        throw GException::out_of_range(G_OP_ACCESS, index, 0, tests()-1);
    }

    // Return reference
    return *(m_tests[index]);
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return Test Suite name
 ***************************************************************************/
std::string GTestSuite::name(void) const
{
    return m_name;
}

/***********************************************************************//**
 * @brief Set Test Suite name
 * @param[in] name Parameter name.
 ***************************************************************************/
void GTestSuite::name(const std::string& name)
{
    // Set name
    m_name=name;
    
    // Return
    return;
}

/***********************************************************************//**
 * @brief Enables/disables logging into standard output stream
 *
 * @param[in] flag Enable/disable logging (true/false).
 *
 * Enables or disables logging into the standard output stream.
 ***************************************************************************/
void GTestSuite::cout(bool cout)
{
    // Enables or disables logging into the standard output stream
    m_log.cout(cout);
    
    //Return
    return;
}

/***********************************************************************//**
 * @brief Test an assert
 * @param[in] assert Assert (true/false).
 * @param[in] name Test case name. (defaults to "")
 * Examples :
 *          GTestSuite testsuite;
 *          testsuite.test_assert(x>3,"Test if x > 3");
 *          testsuite.test_assert(x>3&&x<10,"Test if  3 < x < 10 ");
 ***************************************************************************/
void GTestSuite::test_assert(bool assert, const std::string& name,const std::string& message)
{
    // If it is the first test, print test suite name
    if(tests()==0){
        m_log<<this->name()<<": ";
    }
    
    // Create a test case of failure type
    GTestCase* testcase = new GTestCase(GTestCase::FAIL_TEST,name);
    
    // If assert is false 
    if(!assert)
    {
        // test case is not ok
        testcase->m_passed=false;
        
        //increment number of failures
        m_failures++;
    }
    
    //Set message
    testcase->message(message);
    
    // Show the result ( ".","F" or, "E")
    m_log<<testcase->print_result();
    
    // Add the test in the container
    m_tests.push_back(testcase);
    
    // Return
    return;
    
}

/***********************************************************************//**
 * @brief Test an try block
 * @param[in] name Test case name. (defaults to "")
 * @see test_try_sucess() 
 * @see test_try_failure(const std::string& message)
 * @see test_try_failure(const std::exception& e)
 * Call before testing a try/catch block.
 * Example: 
 *       GTestSuite testsuite;
 *       testsuite.test_try("Test a try block");
 *       try{
 *          ... //someting to test
            testsuite.test_try_success();
 *       }
 *      catch(...)
 *      {
 *          testsuite.test_try_failure();
 *      }
 *
 ***************************************************************************/
void GTestSuite::test_try(const std::string& name)
{ 
    // If it is the first test, print test suite name
    if(tests()==0){
        m_log<<this->name()<<": ";
    }
    
    // Create a test case of error type
    GTestCase* testcase = new GTestCase(GTestCase::ERROR_TEST,name);
    
    //Add test case to container
    m_tests.push_back(testcase);
    
    // Return
    return;
}

/***********************************************************************//**
 * @brief Notice when a try block successed
 * @see test_try(const std::string& name)
 * @see test_try_failure(const std::string& message)
 * @see test_try_failure(const std::exception& e)
 *
 * Call this method at the last line of a try
 * Example: 
 *       GTestSuite testsuite;
 *       testsuite.test_try("Test a try block");
 *       try{
 *          ... //someting to test
            testsuite.test_try_success();
 *       }
 *      catch(...)
 *      {
 *          testsuite.test_try_failure();
 *      }
 *
 ***************************************************************************/
void GTestSuite::test_try_success(){
    // The try block test is ok
    m_log<<m_tests.back()->print_result();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Notice when a try block failed
 * @param[in] message Message to explain why it failed (defaults to "")
 * @param[in] type Type of message (defaults to "")
 * @see test_try_sucess()
 * @see test_try(const std::string& name)
 * @see test_try_failure(const std::exception& e)
 *
 * Call this method in the catch block.
 *
 * Example: 
 *       GTestSuite testsuite;
 *       testsuite.test_try("Test a try block");
 *       try{
 *          ... //someting to test
            testsuite.test_try_success();
 *       }
 *      catch(...)
 *      {
 *          testsuite.test_try_failure("Problem with ...");
 *      }
 *
 ***************************************************************************/
void GTestSuite::test_try_failure(const std::string& message, const std::string& type){
    // Test is not ok
    m_tests.back()->m_passed=false;
    
    // Increment the number of errors
    m_errors++;
    
    // Set a message
    m_tests.back()->message(message);
 
    // Set a type of message
    m_tests.back()->message_type(type);
    
    // Show the result ( ".","F" or, "E")
    m_log<<m_tests.back()->print_result();
    
    //Return
    return;
}

/***********************************************************************//**
 * @brief Notice when a try block failed
 * @param[in] exception Exception caught
 * @see test_try_sucess()
 * @see test_try(const std::string& name)
 * @see test_try_failure(const std::string& message)
 *
 * Call this method in a catch block.
 *
 * Example: 
 *       GTestSuite testsuite;
 *       testsuite.test_try("Test a try block");
 *       try{
 *          ... //someting to test
            testsuite.test_try_success();
 *       }
 *      catch(exception& e)
 *      {
 *          testsuite.test_try_failure(e);
 *      }
 *
 ***************************************************************************/
void GTestSuite::test_try_failure(const std::exception& e){
    
    // Extract message of exception and class name
    test_try_failure(e.what(),typeid(e).name());
    
    //Return
    return;
    
}

/***********************************************************************//**
 * @brief Test a function
 * @param[in] function Pointer to the function to test ( (void)(*)(void) )
 * @param[in] name Test case name. (defaults to "")
 *
 * The function to test should have only one throw instruction else the different tests could not be separate.
 * 
 * Example :
 *
 *      void foo()
 *      {
 *          if (...)
 *              throw GException::bad_alloc("Test");
 *      }
 *
 *     in main :
 *          GTestSuite testsuite;
 *          testsuite.test_function(&foo,"Test function foo");
 *
 * 
 *
 ***************************************************************************/
void GTestSuite::test_function(const pfunction function,const std::string& name)
{
     // If it is the first test, print test suite name
    if(tests()==0){
        m_log<<this->name()<<": ";
    }

    // Create a test case of error type (by default with GTestCase)
    GTestCase* testcase = new GTestCase(function,name);

    // Add test case to container
    m_tests.push_back(testcase);

    // Call the function.
    testcase->run();

    // Imcrement number of errors if the test is not passed
    if(!testcase->is_passed()){
        m_errors++;
    }

    // Show the result ( ".","F" or, "E")
    m_log<<testcase->print_result();

    // Return
    return;
}

/***********************************************************************//**
 * @brief Notice they are not other test to do.
 ***************************************************************************/
void GTestSuite::end_test()
{
    if((m_errors+m_failures)==0){
        m_log<<" ok";
    }
    else{
        m_log<<" NOK\n";
    }
}

/***********************************************************************//**
 * @brief Return the number of tests
 ***************************************************************************/
int GTestSuite::tests() const
{
    return m_tests.size(); 
}

/***********************************************************************//**
 * @brief Return the number of errors
 ***************************************************************************/
int GTestSuite::errors() const
{
    return m_errors; 
}

/***********************************************************************//**
 * @brief Return the number of failures
 ***************************************************************************/
int GTestSuite::failures() const
{
    return m_failures; 
}

/***********************************************************************//**
 * @brief Return the number of success
 ***************************************************************************/
int GTestSuite::success() const
{
    return tests()-(m_errors+m_failures);
}

/***********************************************************************//**
 * @brief Return the timestamp. Set at the creation of the object
 ***************************************************************************/
std::time_t GTestSuite::timestamp() const
{
    return m_timestamp;
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
    m_name="Unamed Test Suite";;
    m_log.clear();
    cout(true);
    m_log.max_size(1);
    m_tests.clear();
    m_errors=0;
    m_failures=0;
    m_timestamp=time(NULL);
    
    // Return
    return;
}

/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] testsuite Test suite.
 *
 * This method just clone the container not the test case.
 ***************************************************************************/
void GTestSuite::copy_members(const GTestSuite& testsuite)
{
    m_name=testsuite.m_name;
    m_tests=testsuite.m_tests;
    m_log=testsuite.m_log;
    m_errors=testsuite.m_errors;
    m_failures=testsuite.m_failures;
    
    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTestSuite::free_members(void)
{
    //Add a new line
    m_log<<"\n";
    
    // Delete test cases
    for(int i=0;i<m_tests.size();++i)
    {
        delete m_tests[i]; 
    }
    // Return
    return;
}
