/***************************************************************************
 *                   GException.cpp  -  exception handlers                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2013 by Juergen Knoedlseder                         *
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
 * @file GException.cpp
 * @brief Exception handler interface implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Exception message.
 ***************************************************************************/
const char* GExceptionHandler::what() const throw()
{
    // Set error message
    std::string message = "*** ERROR in " + m_origin + ": " + m_message;

    // Return message as C character array
    return (gammalib::tochar(message));
}


/***********************************************************************//**
 * @brief Invalid value
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Error message.
 ***************************************************************************/
GException::invalid_value::invalid_value(const std::string& origin,
                                         const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message string
    m_message = "Invalid value.";
    if (message.length() > 0) {
        m_message += ("\n" + message);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid argument
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Error message.
 ***************************************************************************/
GException::invalid_argument::invalid_argument(const std::string& origin,
                                               const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message string
    m_message = "Invalid argument.";
    if (message.length() > 0) {
        m_message += ("\n" + message);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Index is out of range [0,elements-1]
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] what Describes what is out of range.
 * @param[in] index Index.
 * @param[in] elements Number of elements.
 * @param[in] message Optional error message.
 *
 * The @p what string specifies the index type that is out of range. For
 * example what="Vector index" will lead to "Vector index <index> is ...".
 * The first letter in @p what is expected to be a capital letter.
 ***************************************************************************/
GException::out_of_range::out_of_range(const std::string& origin,
                                       const std::string& what,
                                       const int&         index,
                                       const int&         elements,
                                       const std::string& message)
{
    // Set origin
    m_origin = origin;

    // Set message
    if (elements > 0) {
        m_message = gammalib::strip_whitespace(what) + " " +
                    gammalib::str(index) + " is outside the"
                    " valid range [0," + gammalib::str(elements-1) + "].";
    }
    else {
        m_message = "Invalid access to empty object with " +
                    gammalib::tolower(gammalib::strip_whitespace(what)) +
                    " " + gammalib::str(index) + ".";
    }
    if (message.length() > 0) {
        m_message += ("\n" + message);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Feature not implement
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional message.
 *
 * This exception signals features that are not yet implemented. It may be
 * thrown by modules that are still under development.
 ***************************************************************************/
GException::feature_not_implemented::feature_not_implemented(std::string origin,
                                                             std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message string
    if (message.length() > 0) {
        m_message = message;
    }
    else {
        m_message = "Feature not implemented.";
    }
    m_message += " In case that you need this feature for your application"
                 " please submit a feature request on"
                 " https://sourceforge.net/projects/gammalib/,"
                 " join this error message and provide a detailed"
                 " description of your needs.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid argument
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] argument Argument name.
 * @param[in] message Optional message.
 *
 * This exception signals that a specified argument was not valid.
 ***************************************************************************/
GException::invalid_argument::invalid_argument(const std::string& origin,
                                               const std::string& argument,
                                               const std::string& message)
{
    // Set origin
    m_origin = origin;

    // Set message string
    m_message = "Invalid argument \""+argument+"\" value.";

    // Add optional error message
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid type conversion
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::bad_type::bad_type(std::string origin, std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message string
    m_message = "Invalid type conversion.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Environment variable not found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] envname Name of environment variable.
 * @param[in] message Optional message.
 *
 * This exception signals that are required environment variable was not
 * found.
 ***************************************************************************/
GException::env_not_found::env_not_found(std::string origin,
                                         std::string envname,
                                         std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message string
    m_message = "Environment variable \""+envname+"\" not found.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Memory allocation exception.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] num Requested amount of memory.
 ***************************************************************************/
GException::mem_alloc::mem_alloc(std::string origin, unsigned num)
{
    // Set origin
    m_origin = origin;

    // Set message string
    m_message = "Memory allocation error (" + gammalib::str((int)num) + " elements)";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Not enough nodes in node array
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] num Number of nodes in node array.
 ***************************************************************************/
GException::not_enough_nodes::not_enough_nodes(std::string origin, int num)
{
    // Set origin
    m_origin = origin;

    // Set message string
    m_message = "Not enough nodes in node array (" + gammalib::str(num) + " nodes).";

    // Return
    return;
}


/***********************************************************************//**
 * @brief File not found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Filename.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::file_not_found::file_not_found(std::string origin,
                                           std::string filename,
                                           std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message string
    m_message = "File \"" + filename +"\" not found.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief File open error
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Filename.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::file_open_error::file_open_error(std::string origin,
                                             std::string filename,
                                             std::string message)
{
    // Set origin
    m_origin  = origin;

    // Set message string
    m_message = "Unable to open file \"" + filename +"\"";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Directory not found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] dirname Directory name.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::directory_not_found::directory_not_found(std::string origin,
                                                     std::string dirname,
                                                     std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message string
    m_message = "Directory \"" + dirname +"\" not found.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Directory not accessible
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] dirname Directory name.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::directory_not_accessible::directory_not_accessible(std::string origin,
                                                               std::string dirname,
                                                               std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message string
    m_message = "Directory \"" + dirname +"\" not accessible.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bad CSV columns encountered
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Filename.
 * @param[in] row CSV file row.
 * @param[in] cols Expected number of columns.
 * @param[in] elements Found number of columns.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::csv_bad_columns::csv_bad_columns(std::string origin,
                                             std::string filename,
                                             int         row,
                                             int         cols,
                                             int         elements,
                                             std::string message)
{
    m_origin  = origin;
    m_message = "Inconsistent number of columns found in row "+gammalib::str(row)+
                " of file \""+filename+"\" ("+gammalib::str(cols)+" expected, "+
                gammalib::str(elements)+" encountered). "+message;
    return;
}
