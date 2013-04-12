/***************************************************************************
 *                 GCTAException.cpp  - CTA exception handler              *
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
 * @file GCTAException.cpp
 * @brief CTA exception handler interface implementation.
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GCTAException.hpp"
#include "GTools.hpp"
#include "GCTAModelRadialRegistry.hpp"


/***********************************************************************//**
 * @brief Exception message.
 ***************************************************************************/
const char* GCTAExceptionHandler::what() const throw()
{
    // Set error message
    std::string message = "*** ERROR in " + m_origin + ": " + m_message;

    // Return message as C character array
    return gammalib::tochar(message);
}


/***********************************************************************//**
 * @brief Error while opening file.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Name of file that could not be opened.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::file_open_error::file_open_error(std::string origin,
                                                std::string filename,
                                                std::string message)
{
    m_origin  = origin;
    m_message = "Unable to open file '"+filename+"'. "+message;
    return;
}


/***********************************************************************//**
 * @brief Member not set.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Error message.
 ***************************************************************************/
GCTAException::no_member::no_member(const std::string& origin,
                                    const std::string& message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = message;

    // Return
    return;
}




/***********************************************************************//**
 * @brief Pointing not set.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::no_pointing::no_pointing(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No pointing found in CTA observation. "+message;
    return;
}


/***********************************************************************//**
 * @brief Instrument response not set.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::no_response::no_response(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "Attempt to use CTA response, but no response function"
                " has been found. Use GCTAObservation::response() method"
                " to set response function. "+message;
    return;
}


/***********************************************************************//**
 * @brief No sky pixels found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::no_sky::no_sky(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No sky pixels have been found. "+message;
    return;
}


/***********************************************************************//**
 * @brief No energy boundary information found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::no_ebds::no_ebds(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No energy boundaries have been found. "+message;
    return;
}


/***********************************************************************//**
 * @brief No Good Time Intervals found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::no_gti::no_gti(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No Good Time Intervals (GTIs) have been found. "+message;
    return;
}


/***********************************************************************//**
 * @brief No energies set
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::no_energies::no_energies(std::string origin, std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Energy vector has not been setup."
                " Cannot access event information.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief No sky directions set
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::no_dirs::no_dirs(std::string origin, std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Sky direction vector has not been setup."
                " Cannot access event information.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observation is not a CTA observation
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_observation_type::bad_observation_type(std::string origin,
                                                          std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Observation is not of type GCTAObservation.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Event is not a CTA event
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_event_type::bad_event_type(std::string origin, std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Event is not of type GCTAEventAtom or GCTAEventBin.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief ROI is not a CTA ROI
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_roi_type::bad_roi_type(std::string origin, std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "ROI is not of type GCTARoi.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Instrument direction is not a CTA instrument direction
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_instdir_type::bad_instdir_type(std::string origin,
                                                  std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Instrument direction is not of type GCTAInstDir.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Pointing is not a CTA pointing
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_pointing_type::bad_pointing_type(std::string origin,
                                                    std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Pointing is not of type GCTAPointing.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response is not a CTA response
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_response_type::bad_response_type(std::string origin,
                                                    std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Response is not of type GCTAResponse.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid radial model type
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Radial model type that has been encountered.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::model_invalid_radial::model_invalid_radial(std::string origin,
                                                          std::string type,
                                                          std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid radial CTA model type \""+type+"\" encountered. " +
                message;

    // Add list of valid spatial models
    GCTAModelRadialRegistry registry;
    if (registry.size() > 0) {
        m_message += "The following models are registered: ";
        for (int i = 0; i < registry.size(); ++i) {
            if (i > 0) {
                m_message += ", ";
            }
            m_message += "\"" + registry.name(i) + "\"";
        }
        m_message += ".";
    }
    else {
        m_message += "No models are registered.";
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response is not a CTA response
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_model_type::bad_model_type(const std::string& origin,
                                              const std::string& message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Model is not of expected type.";

    // Optional message
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bad response table dimension
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] found Found table dimension.
 * @param[in] expected Expected table dimension.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_rsp_table_dim::bad_rsp_table_dim(std::string origin,
                                                    int         found,
                                                    int         expected,
                                                    std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Response table dimension "+
                gammalib::str(found)+
                " is smaller than the expected dimension "+
                gammalib::str(expected)+
                ".";
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bad response table format
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::bad_rsp_table_format::bad_rsp_table_format(std::string origin,
                                                          std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Response table has invalid format.";
    if (message.length() > 0) {
        m_message += " "+message;
    }

    // Return
    return;
}
