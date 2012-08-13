/***************************************************************************
 *               GException_test.cpp  -  Unit Test exception handlers      *
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
 * @file GException_test.cpp
 * @brief Implement exceptions for unit test module
 * @author Jean-Baptiste Cayrou
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Test nested try error
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::test_nested_try_error::test_nested_try_error(std::string origin, std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Nested try error ("+message+")";

    // Return
    return;
}

/***********************************************************************//**
 * @brief Failure test
 * @param[in] message Optional error message.
 * Use in unit tests to notice a fail in a try block
 ***************************************************************************/
GException::test_failure::test_failure(std::string origin, std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Failure Test ("+message+")";

    // Return
    return;
}

/***********************************************************************//**
 * @brief Failure test
 * @param[in] message Optional error message.
 * Use in unit tests to notice a fail in a try block
 ***************************************************************************/
GException::test_error::test_error(std::string origin, std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Error Test ("+message+")";

    // Return
    return;
}
