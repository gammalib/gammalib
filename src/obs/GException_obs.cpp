/***************************************************************************
 *          GException_obs.cpp  -  observations exception handlers         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Juergen Knoedlseder                         *
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
 * @file GException_obs.cpp
 * @brief Observation exception handler interface implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief No valid instrument response found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_response::no_response(std::string origin, std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "No valid instrument response found in observation." \
                " Set instrument response before using this method.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief No valid event container found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_events::no_events(std::string origin, std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "No valid event container found in observation." \
                " Set event container before using this method.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief No valid event list found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_list::no_list(std::string origin, std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "No valid event list found in observation." \
                " Set event list before using this method.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief No valid event cube found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_cube::no_cube(std::string origin, std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "No valid event cube found in observation." \
                " Set event cube before using this method.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief No valid region of interest found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::no_roi::no_roi(std::string origin, std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "No valid region of interest (ROI) found in observation." \
                " Set region of interest before using this method.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Mismatch between gradient vector and number of model parameters
 *
 * @param[in] origin Method that throws the error.
 * @param[in] nsize Gradient vector size.
 * @param[in] npars Number of parameters in model.
 ***************************************************************************/
GException::gradient_par_mismatch::gradient_par_mismatch(std::string origin,
                                                         int         nsize,
                                                         int         npars)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "Gradient vector size "+str(nsize)+
                " mismatches number "+str(npars)+" of model parameters.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Calibration database directory not found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] caldb Name of calibration database.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::caldb_not_found::caldb_not_found(std::string origin,
                                             std::string caldb,
                                             std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "Calibration database \""+caldb+"\" not found.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response error: invalid response type specified
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Specified response type.
 ***************************************************************************/
GException::rsp_invalid_type::rsp_invalid_type(std::string origin,
                                               std::string type)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "Invalid response type \""+type+"\" specified.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid Good Time Interval found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] tstart Good Time Interval start time.
 * @param[in] tstop Good Time Interval stop time.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::gti_invalid::gti_invalid(std::string origin, GTime tstart,
                                     GTime tstop, std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "Invalid Good Time Interval (" + tstart.print() + "-" +
                tstop.print() + ") specified.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid energy range found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] emin Minimum energy in MeV.
 * @param[in] emax Maximum energy in MeV.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::erange_invalid::erange_invalid(std::string origin, double emin, 
                                           double emax, std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "Invalid energy range (Emin="+str(emin)+" MeV, Emax="+
                str(emax)+" MeV specified.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid optimization statistics
 *
 * @param[in] origin Method that throws the error.
 * @param[in] statistics Encountered statistics.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::invalid_statistics::invalid_statistics(std::string origin,
                                                   std::string statistics,
                                                   std::string message)
{
    // Set method name
    m_origin = origin;

    // Set error message
    m_message = "Invalid optimization statistics \""+statistics+
                "\" specified.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid instrument encountered
 *
 * @param[in] origin Method that throws the error.
 * @param[in] instrument Encountered instrument.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::invalid_instrument::invalid_instrument(std::string origin,
                                                   std::string instrument,
                                                   std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid instrument \""+instrument+"\" encountered.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Add list of valid instruments
    GObservationRegistry registry;
    if (registry.size() > 0) {
        m_message += "The following instruments are registered: ";
        for (int i = 0; i < registry.size(); ++i) {
            if (i > 0) {
                m_message += ", ";
            }
            m_message += "\"" + registry.name(i) + "\"";
        }
        m_message += ".";
    }
    else {
        m_message += "No instruments are registered.";
    }

    // Return
    return;
}
