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
#define G_TRY_SUCCESS                       "GTestSuite::test_try_success()"
#define G_TRY_FAILURE                       "GTestSuite::test_try_failure(std::string&,std::string&)"

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
 * @brief Run tests
 ***************************************************************************/
bool GTestSuite::run(void)
{
    // Set
    set();

    //Show the TestSuite name
    m_log<<this->name()<<": ";

    bool was_successful = true;

    while(m_stack_test.size()>0)
    {
        if(m_stack_test.front()->ptr_function()!=NULL){

            // Add the test in the m_tests container
            m_tests.push_back(m_stack_test.front());

            // Run the test
            m_stack_test.front()->run();

            // Imcrement number of errors if the test is not passed
            if(!m_stack_test.front()->is_passed()){
                m_errors++;
            }

            // Show the result ( ".","F" or, "E")
            m_log<<m_stack_test.front()->print_result();

            // Erase the test case
            m_stack_test.erase(m_stack_test.begin());
        }
    }

    // If they are errors or failures
    if((m_errors+m_failures)==0){
        m_log<<" ok\n";
        was_successful=true;
    }
    else{
        m_log<<" NOK\n";
        was_successful=false;
    }

    return was_successful;
}

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
 *          test_assert(x>3,"Test if x > 3");
 *          test_assert(x>3&&x<10,"Test if  3 < x < 10 ");
 ***************************************************************************/
void GTestSuite::test_assert(bool assert, const std::string& name,const std::string& message)
{

    // Create a test case of failure type
    GTestCase* testcase = new GTestCase(GTestCase::FAIL_TEST,format_name(name));

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
 * @see test_try_failure(const std::string& message,const std::string& type)
 * @see test_try_failure(const std::exception& e)
 * Call before testing a try/catch block.
 * Example: 
 *       test_try("Test a try block");
 *       try{
 *          ... //someting to test
            test_try_success();
 *       }
 *      catch(...)
 *      {
 *          test_try_failure();
 *      }
 *
 ***************************************************************************/
void GTestSuite::test_try(const std::string& name)
{
    // Create a test case of error type
    GTestCase* testcase = new GTestCase(GTestCase::ERROR_TEST,format_name(name));

    //Add test case to container
    m_stack_try.push_back(testcase);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Notice when a try block successed
 * @see test_try(const std::string& name)
 * @see test_try_failure(const std::string& message, const std::string& message)
 * @see test_try_failure(const std::exception& e)
 *
 * Call this method at the last line of a try
 * Example: 
 *       test_try("Test a try block");
 *       try{
 *          ... //someting to test
            test_try_success();
 *       }
 *      catch(...)
 *      {
 *          test_try_failure();
 *      }
 *
 ***************************************************************************/
void GTestSuite::test_try_success(void){

    // If the stack is empty
    if(m_stack_try.size()==0)
        throw GException::test_nested_try_error(G_TRY_SUCCESS, "Error : test_try_success() without test_try()");

    // Add the test in the container
    m_tests.push_back(m_stack_try.back());

    // Delete the test case of the stack
    m_stack_try.pop_back();

    // The try block test is ok
    m_log<<m_tests.back()->print_result();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Notice when a try block failed
 * @param[in] message Message to explain why it failed (defaults to "")
 * @param[in] message_type Type of message (defaults to "")
 * @param[in] type Type of case (GTestCase::ERROR_TEST or GTestCase::FAIL_TEST)
 * @see test_try_sucess()
 * @see test_try(const std::string& name)
 * @see test_try_failure(const std::exception& e)
 *
 * Call this method in the catch block.
 *
 * Example: 
 *       test_try("Test a try block");
 *       try{
 *          ... //someting to test
            test_try_success();
 *       }
 *      catch(...)
 *      {
 *          test_try_failure("Problem with ...");
 *      }
 *
 ***************************************************************************/
void GTestSuite::test_try_failure(const std::string& message, const std::string& message_type){

    // If the stack is empty
    if(m_stack_try.size()==0)
        throw GException::test_nested_try_error(G_TRY_FAILURE,"test_try_success() call without test_try()");

    // Test is not ok
    m_stack_try.back()->m_passed=false;

    // Increment the number of errors
    m_errors++;

    //Set type
    m_stack_try.back()->type(GTestCase::ERROR_TEST);
    
    // Set a message
    m_stack_try.back()->message(message);

    // Set a type of message
    m_stack_try.back()->message_type(message_type);

    // Add the test in the container
    m_tests.push_back(m_stack_try.back());
    
    // Delete the test case of the stack
    m_stack_try.pop_back();

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
 * @brief Return a failure exception
 * @param[in] message Message of the exception
 * @see test_try()
 * @see test_error(const std::string& message)
 *
 * It can be use in a try test
 *
 * Example: 
 *       GTestSuite testsuite;
 *       testsuite.test_try("Test a try block");
 *       try{
 *          throw testsuite.exception_failure("a failure");
            testsuite.test_try_success();
 *       }
 *      catch(exception& e)
 *      {
 *          testsuite.test_try_failure(e);
 *      }
 *
 ***************************************************************************/

GException::test_failure& GTestSuite::exception_failure(const std::string& message)
{
    return *(new GException::test_failure(m_stack_try.back()->name(),message));
}

/***********************************************************************//**
 * @brief Return an error exception
 * @param[in] message Message of the exception
 * @see test_try()
 * @see test_failure(const std::string& message)
 *
 * It can be use in a try test
 *
 * Example: 
 *       testsuite;
 *       ttest_try("Test a try block");
 *       try{
 *          throw exception_error("an error");
            test_try_success();
 *       }
 *      catch(exception& e)
 *      {
 *          test_try_failure(e);
 *      }
 *
 ***************************************************************************/
 
GException::test_error& GTestSuite::exception_error(const std::string& message)
{
    return *(new GException::test_error(m_stack_try.back()->name(),message));
}

void GTestSuite::add_test(const pfunction function,const std::string& name)
{
    GTestCase * testcase = new GTestCase(function,name,this);

    // Add test case to container
    m_stack_test.push_back(testcase);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Return the number of tests
 ***************************************************************************/
int GTestSuite::tests(void) const
{
    return m_tests.size(); 
}

/***********************************************************************//**
 * @brief Return the number of errors
 ***************************************************************************/
int GTestSuite::errors(void) const
{
    return m_errors; 
}

/***********************************************************************//**
 * @brief Return the number of failures
 ***************************************************************************/
int GTestSuite::failures(void) const
{
    return m_failures; 
}

/***********************************************************************//**
 * @brief Return the number of success
 ***************************************************************************/
int GTestSuite::success(void) const
{
    return tests()-(m_errors+m_failures);
}

/***********************************************************************//**
 * @brief Return the timestamp. Set at the creation of the object
 ***************************************************************************/
time_t GTestSuite::timestamp(void) const
{
    return m_timestamp;
}

/***********************************************************************//**
 * @brief Return the total duration of the tests
 ***************************************************************************/
double GTestSuite::duration(void) const
{
    double time=0;
    for(int i=0;i<m_tests.size();i++){
        time+=m_tests[i]->time();
    }

    return time;
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
    // Delete test cases
    for(int i=0;i<m_tests.size();++i)
    {
        delete m_tests[i]; 
    }
    // Return
    return;
}

/***********************************************************************//**
 * @brief Format Name
 * Return a string with the format: "TestFunctionName:TestTryname1:TestTryName2: name"
 ***************************************************************************/
std::string GTestSuite::format_name(const std::string& name)
{
    std::string format_name;

    // Add the names of the try blocks
    if(m_stack_try.size()>0){
        format_name+=m_stack_try.back()->name();
    }
    else
    {
         // Add the name of the test before
        if(m_stack_test.size()>0){
            format_name=(m_stack_test.front())->name();
        }

    }

    format_name+=": "+name;

    return format_name;
}
