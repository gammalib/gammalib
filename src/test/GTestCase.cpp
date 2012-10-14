/***************************************************************************
 *                    GTestCase.cpp  - Test case class                     *
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
 * @param[in] test Test case.
 ***************************************************************************/
GTestCase::GTestCase(const GTestCase& test)
{
    // Initialise members
    init_members();
    
    // Copy members
    copy_members(test);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Type and name constructor
 *
 * @param[in] kind Test kind (FAIL_TEST or ERROR_TEST)
 * @param[in] name Test case name (default: "")
 ***************************************************************************/
GTestCase::GTestCase(ErrorKind kind, const std::string& name)
{
    // Initialise members
    init_members();
    
    // Set kind
    m_kind = kind;
    
    // Set name
    m_name = name;

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
 * @param[in] test Test Case.
 ***************************************************************************/
GTestCase&  GTestCase::operator= (const GTestCase& test)
{
    // Execute only if object is not identical
    if (this != &test) {
        
        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(test);

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
 * @brief Clear instance
 ***************************************************************************/
void GTestCase::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GTestCase* GTestCase::clone(void) const
{
    // Clone this image
    return new GTestCase(*this);
}


/***********************************************************************//**
 * @brief Return test case name
 ***************************************************************************/
std::string GTestCase::name(void) const
{
    // Return name
    return m_name;
}


/***********************************************************************//**
 * @brief Set test case name
 *
 * @param[in] name Test case name.
 ***************************************************************************/
void GTestCase::name(const std::string& name)
{
    // Set name
    m_name = name;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return test case message
 ***************************************************************************/
std::string GTestCase::message(void) const
{
    // Return message
    return m_message;
}


/***********************************************************************//**
 * @brief Set test case message
 *
 * @param[in] message Test case message.
 ***************************************************************************/
void GTestCase::message(const std::string& message)
{
    // Set message
    m_message = message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return test case type
 ***************************************************************************/
std::string GTestCase::type(void) const
{
    // Return test case type
    return m_type;
}


/***********************************************************************//**
 * @brief Set type of test case
 *
 * @param[in] type Type of test case.
 ***************************************************************************/
void GTestCase::type(const std::string& type)
{
    // Set type of test case
    m_type = type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return kind of test case
 *
 * Returns whether this test case is for failure testing (FAIL_TEST) or
 * error testing (ERROR_TEST).
 ***************************************************************************/
GTestCase::ErrorKind GTestCase::kind(void) const
{
    // Return kind
    return m_kind; 
}


/***********************************************************************//**
 * @brief Set kind of test case
 *
 * @param[in] kind Kind of test case (FAIL_TEST or ERROR_TEST)
 *
 * Specifies whether this test case is for failure testing (FAIL_TEST) or
 * error testing (ERROR_TEST).
 ***************************************************************************/
void GTestCase::kind(ErrorKind kind)
{
    // Set kind
    m_kind = kind;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return whether the test passed
 *
 * This method returns true if the test has passed, false otherwise.
 ***************************************************************************/
bool GTestCase::passed(void) const
{
    // Return
    return m_passed;
}


/***********************************************************************//**
 * @brief Set if test passed
 *
 * @param[in] passed Test passed (true or false)
 ***************************************************************************/
void GTestCase::passed(const bool& passed)
{
    // Set passed flag
    m_passed = passed;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return test case duration
 ***************************************************************************/
double GTestCase::duration(void) const
{
    // Return duration
    return m_duration;
}


/***********************************************************************//**
 * @brief Set test duration
 *
 * @param[in] duration Test duration.
 ***************************************************************************/
void GTestCase::duration(const double& duration)
{
    // Set duration
    m_duration = duration;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print test case result
 *
 * Returns either ".", "F" or "E", dependent on the test result.
 ***************************************************************************/
std::string GTestCase::print(void) const
{
    // Initialize string
    std::string result;

    // Set string dependent on test result
    if (m_passed) {
        result.append(".");
    }
    else if (m_kind == ERROR_TEST) {
        result.append("E");
    }
    else {
        result.append("F");
    }

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
void GTestCase::init_members(void)
{
    // Initialise members
    m_name     = "Unnamed Error Test Case";
    m_message  = "";
    m_type     = "";
    m_passed   = true;
    m_kind     = ERROR_TEST;
    m_duration = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] test Test case.
 *
 * Copies all members from the test case class.
 ***************************************************************************/
void GTestCase::copy_members(const GTestCase& test)
{
    // Copy members
    m_name     = test.m_name;
    m_message  = test.m_message;
    m_type     = test.m_type;
    m_passed   = test.m_passed;
    m_kind     = test.m_kind;
    m_duration = test.m_duration;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTestCase::free_members(void)
{
    // Return
    return;
}
