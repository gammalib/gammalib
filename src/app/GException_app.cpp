/***************************************************************************
 *           GException_app.cpp  -  Application exception handlers         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GException_app.cpp
 * @brief Application exception handler interface definition
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Generic application error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::app_error::app_error(const std::string& origin,
                                 const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter file not found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename Filename that was not found.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::par_file_not_found::par_file_not_found(const std::string& origin,
                                                   const std::string& filename,
                                                   const std::string& message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Parameter file \""+filename+"\" not found. "
                "Make sure that PFILES environment variable is set "
                "correctly.";

    // Add optional error message
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unable to open parameter file
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename Filename that could not be opened.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::par_file_open_error::par_file_open_error(const std::string& origin,
                                                     const std::string& filename,
                                                     const std::string& message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Unable to open parameter file \""+filename+"\".";

    // Add optional error message
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unable to determine users home directory
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::home_not_found::home_not_found(const std::string& origin,
                                           const std::string& message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Unable to determine users home directory.";

    // Add optional error message
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unable to create pfiles directory in users home directory
 *
 * @param[in] origin Method that throws the error.
 * @param[in] home Users home directory.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::could_not_create_pfiles::could_not_create_pfiles(const std::string& origin,
                                                             const std::string& home,
                                                             const std::string& message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Unable to create \""+home+"\".";

    // Add optional error message
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unable to change access rights for pfiles directory
 *
 * @param[in] origin Method that throws the error.
 * @param[in] home Users home directory.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::pfiles_not_accessible::pfiles_not_accessible(const std::string& origin,
                                                         const std::string& home,
                                                         const std::string& message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Could not make \""+home+"\" write accessible for "
                "writing of the applications parameter file.";

    // Add optional error message
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Syntax error encountered in parameter file line
 *
 * @param[in] origin Method that throws the error.
 * @param[in] line Parameter file line in which syntax error occured.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::par_file_syntax_error::par_file_syntax_error(const std::string& origin,
                                                         const std::string& line,
                                                         const std::string& message)
{
    // Set origin and message
    m_origin  = origin;
    if (message.length() > 0) {
        m_message = "Syntax error occured in the following line of the "
                    "parameter file ("+message+"): "+line;
    }
    else {
        m_message = "Syntax error occured in the following line of the "
                    "parameter file: "+line;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Error encountered in parameter definition
 *
 * @param[in] origin Method that throws the error.
 * @param[in] name Parameter name.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::par_error::par_error(const std::string& origin,
                                 const std::string& name,
                                 const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Parameter \""+name+"\"";

    // Add optional error message
    if (message.length() > 0) {
        m_message += ": " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid command line parameter
 *
 * @param[in] origin Method that throws the error.
 * @param[in] arg Command line argument.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::bad_cmdline_argument::bad_cmdline_argument(const std::string& origin,
                                                       const std::string& arg,
                                                       const std::string& message)
{
    // Set origin and message
    m_origin  = origin;
    if (message.length() > 0) {
        m_message = "Invalid command line parameter encountered ("+message+
                    "): "+arg;
    }
    else {
        m_message = "Invalid command line parameter encountered: "+arg;
    }

    // Return
    return;
}

