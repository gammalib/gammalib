/***************************************************************************
 *                    GTestCase.cpp - Test case class                      *
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
 * @brief Type and name constructor
 *
 * @param[in] kind Test kind (FAIL_TEST or ERROR_TEST)
 * @param[in] name Test case name (default: "")
 ***************************************************************************/
GTestCase::GTestCase(const ErrorKind& kind, const std::string& name)
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
 * @param[in] test Test case.
 * @return Test case.
 ***************************************************************************/
GTestCase&  GTestCase::operator=(const GTestCase& test)
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
 * @brief Clear test case
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
 * @brief Clone test case
 *
 * @return Pointer to deep copy of test case.
 ***************************************************************************/
GTestCase* GTestCase::clone(void) const
{
    // Clone this image
    return new GTestCase(*this);
}


/***********************************************************************//**
 * @brief Print test case result
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing test case information.
 *
 * Returns either ".", "F" or "E", dependent on the test result.
 ***************************************************************************/
std::string GTestCase::print(const GChatter& chatter) const
{
    // Initialize string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

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
