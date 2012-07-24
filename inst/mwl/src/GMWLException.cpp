/***************************************************************************
 *          GMWLException.cpp  - Multi-wavelength exception handler        *
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
 * @file GMWLException.cpp
 * @brief MWL exception handler interface implementation.
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GMWLException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Exception message
 ***************************************************************************/
const char* GMWLExceptionHandler::what() const throw()
{
    // Set error message
    std::string message = "*** ERROR in " + m_origin + ": " + m_message;

    // Return message as C character array
    return tochar(message);
}


/***********************************************************************//**
 * @brief Error while opening file
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Name of file that could not be opened.
 * @param[in] message Optional error message.
 ***************************************************************************/
GMWLException::file_open_error::file_open_error(std::string origin,
                                                std::string filename,
                                                std::string message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Unable to open file \""+filename+"\".";

    // Append optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief File has invalid format
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GMWLException::bad_file_format::bad_file_format(std::string origin,
                                                std::string message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "File has invalid format.";

    // Append optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid or unsupported unit encountered
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] unit Unit string.
 * @param[in] message Optional error message.
 ***************************************************************************/
GMWLException::invalid_unit::invalid_unit(std::string origin,
                                          std::string unit,
                                          std::string message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Invalid or unsupported unit \""+unit+"\" encountered.";

    // Append optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response is not a MWL response
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GMWLException::bad_response_type::bad_response_type(std::string origin,
                                                    std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Response is not of type GMWLResponse.";

    // Append optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}
