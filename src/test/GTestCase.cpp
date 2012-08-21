/***************************************************************************
 *         GTestCase.cpp  - Test case class for GammaLib                 *
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
 * @file GTestCase.cpp
 * @brief Test Case class implementation
 * @author Jean-Baptiste Cayrou
 */
 
/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTestCase.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#ifdef __APPLE__ & __MACH__
#include <pthread.h>
pthread_attr_t gomp_thread_attr;
#endif
#endif

/* __ Method name definitions ____________________________________________ */

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
 GTestCase::GTestCase(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}

/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] testsuite Test Suite.
 ***************************************************************************/
GTestCase::GTestCase(const GTestCase& testcase)
{
    // Initialise members
    init_members();
    
    // Copy members
    copy_members(testcase);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Type constructor
 *
 * @param[in] type Type of the Test case ( GTestCase::FAIL_TEST or GTestCase::ERROR_TEST )
 ***************************************************************************/
GTestCase::GTestCase(ErrorType type)
{
    //Initialise members
    init_members();

    // Set type
    m_type = type;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Function and name constructor
 *
 * @param[in] ptr_function Pointer to function (void)(*)(void)
 * @param[in] name Name of the Test case
 ***************************************************************************/
GTestCase::GTestCase(const pfunction ptr_function, const std::string& name, GTestSuite* testsuite)
{
    // Initialise members
    init_members();
    
    // Set type as ERROR (test a function)
    m_type = ERROR_TEST;
    
    // Set the pointer function
    m_ptr_function=ptr_function;
    
    // Set name
    m_name=name;

    // Set testsuite
    m_testsuite = testsuite;
    // Return
    return;
}

/***********************************************************************//**
 * @brief Type and name constructor
 *
 * @param[in] type Type of the Test case (GTestCase::FAIL_TEST or GTestCase::ERROR_TEST)
 * @param[in] name Name of the Test case
 ***************************************************************************/
GTestCase::GTestCase(ErrorType type, const std::string& name)
{
    //Initialise members
    init_members();
    
    // Set type
    m_type = type;
    
    // Set name
    m_name=name;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GTestCase::~GTestCase(void)
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
 * @param[in] testcase Test Case.
 ***************************************************************************/
 GTestCase&  GTestCase::operator= (const GTestCase& testcase)
{
    // Execute only if object is not identical
    if (this != &testcase) {
        
        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(testcase);

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
 * @brief Return Test Case name
 ***************************************************************************/
std::string GTestCase::name(void) const{
    return m_name;
}

/***********************************************************************//**
 * @brief Set Test Case name
 * @param[in] name Parameter name.
 ***************************************************************************/
void GTestCase::name(const std::string& name){
    m_name=name;
    return;
}

/***********************************************************************//**
 * @brief Return Test Case message
 ***************************************************************************/
std::string GTestCase::message(void) const{
    return m_message;
}

/***********************************************************************//**
 * @brief Set Test Case mesage
 * @param[in] message Parameter message.
 ***************************************************************************/
void GTestCase::message( const std::string& message){
    m_message=message;
    return;
}

/***********************************************************************//**
 * @brief Return Test Case type of message
 ***************************************************************************/
std::string GTestCase::message_type(void) const{
    return m_message_type;
}

/***********************************************************************//**
 * @brief Set Test Case type of message
 * @param[in] message_type Parameter type of message.
 ***************************************************************************/
void GTestCase::message_type( const std::string& message_type){
    m_message_type=message_type;
    return;
}

/***********************************************************************//**
 * @brief Set Test Case pointer to function
 * @param[in] function Pointer to function (void)(*)(void)
 ***************************************************************************/
void GTestCase::ptr_function(const pfunction function){
    m_ptr_function=function;
    return;
}

/***********************************************************************//**
 * @brief Return Test Case pointer to function
 ***************************************************************************/
GTestSuite& GTestCase::testsuite(void) const{
    return *m_testsuite;
}

/***********************************************************************//**
 * @brief Set Test Case pointer to function
 * @param[in] function Pointer to function (void)(*)(void)
 ***************************************************************************/
void GTestCase::testsuite(GTestSuite* testsuite){
    m_testsuite = testsuite;
    return;
}

/***********************************************************************//**
 * @brief Return Test Case pointer to function
 ***************************************************************************/
pfunction GTestCase::ptr_function(void) const{
    return m_ptr_function;
}

/***********************************************************************//**
 * @brief Return Test Case type
 ***************************************************************************/
GTestCase::ErrorType GTestCase::type(void) const{
    return m_type; 
}

/***********************************************************************//**
 * @brief Set Test Case type
 * @param[in] type Type of the test case
 ***************************************************************************/
void GTestCase::type(ErrorType type){
    m_type=type;
    return;
}

/***********************************************************************//**
 * @brief Return True if the test passed
 ***************************************************************************/
bool GTestCase::is_passed(void) const
{
    return m_passed;
}
/***********************************************************************//**
 * @brief Return the duration of the test case
 ***************************************************************************/
double GTestCase::time(void) const
{
    return m_time;
}
/***********************************************************************//**
 * @brief Print test case result
 * Return a string : ".", "F" or "E"
 ***************************************************************************/
std::string GTestCase::print_result(void) const{
    std::string result;
    if(m_passed){
        result.append(".");
    }
    else if (m_type==ERROR_TEST){
        result.append("E");
    }
    else{
        result.append("F");
    }

    return result;
}

/***********************************************************************//**
 * @brief Call the function m_ptr_function and catch exception.
 ***************************************************************************/
void GTestCase::run(void)
{
    #ifdef _OPENMP
        double t_start = omp_get_wtime();
    #else
        clock_t t_start = clock();
    #endif

    try{
        // If they are a function
        if(m_ptr_function!=NULL){
            // Call it
            (*m_testsuite.*(m_ptr_function))();
        }
    }
    catch(std::exception& e)
    {
        // Test not success
        m_passed=false;

        //Set message
        m_message=e.what();

        //Set type as class name
        m_message_type=typeid(e).name();
    }
    catch(...)
    {
        //For other exceptions
        m_passed=false;
    }

    #ifdef _OPENMP
        double t_elapse = omp_get_wtime()-t_start;
    #else
        double t_elapse = (double)(clock() - t_start) / (double)CLOCKS_PER_SEC;
    #endif

    //Set duration
    m_time=t_elapse;

    return;
}

/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GTestCase::init_members(void)
{
    m_name="Unamed Test Case";
    m_type=ERROR_TEST;
    m_ptr_function=NULL;
    m_passed=true;
    m_time=0;

    return;
}

/***********************************************************************//**
 * @brief Copy class members
 * @param[in] testcase Test case to copy.
 ***************************************************************************/
void GTestCase::copy_members(const GTestCase& testcase)
{
    m_name = testcase.name();
    m_ptr_function = testcase.ptr_function();
    m_message = testcase.m_message;
    m_message_type = testcase.m_message_type;

    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTestCase::free_members(void)
{
    return;
}
