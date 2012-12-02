/***************************************************************************
 *              GCOMException.cpp  - COMPTEL exception handler             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCOMException.cpp
 * @brief COMPTEL exception handler interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GCOMException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Exception message
 ***************************************************************************/
const char* GCOMExceptionHandler::what() const throw()
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
GCOMException::file_open_error::file_open_error(const std::string& origin,
                                                const std::string& filename,
                                                const std::string& message)
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
 * @brief No sky pixels found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCOMException::no_sky::no_sky(const std::string& origin,
                              const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "No sky pixels have been found.";

    // Append optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief No energy boundaries found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCOMException::no_ebds::no_ebds(const std::string& origin,
                                const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "No energy boundaries have been found.";

    // Append optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief No Good Time Intervals found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCOMException::no_gti::no_gti(const std::string& origin,
                              const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "No Good Time Intervals (GTIs) have been found.";

    // Append optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief No sky directions found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCOMException::no_dirs::no_dirs(const std::string& origin,
                                const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Sky direction vector has not been setup."
                " Cannot access event information.";

    // Append optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observation is not a COMPTEL observation
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCOMException::bad_observation_type::bad_observation_type(const std::string& origin,
                                                          const std::string& message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Observation is not of type GCOMObservation.";
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Event is not a COMPTEL event
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCOMException::bad_event_type::bad_event_type(const std::string& origin,
                                              const std::string& message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Event is not of type GCOMEventAtom or GCOMEventBin.";
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}
