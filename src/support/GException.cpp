/***************************************************************************
 *                   GException.cpp  -  exception handlers                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2021 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsCfitsio.hpp"
#include "GEnergy.hpp"


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
        m_message += (" " + message);
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
        m_message += (" " + message);
    }

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
        m_message += (" " + message);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid return value
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Error message.
 ***************************************************************************/
GException::invalid_return_value::invalid_return_value(const std::string& origin,
                                                       const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message string
    m_message = "Invalid return value.";
    if (message.length() > 0) {
        m_message += (" " + message);
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
        m_message += (" " + message);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Runtime error
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Error message.
 ***************************************************************************/
GException::runtime_error::runtime_error(const std::string& origin,
                                         const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message string
    m_message = "Runtime error.";
    if (message.length() > 0) {
        m_message += (" " + message);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief General FITS error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] status cfitsio status.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::fits_error::fits_error(const std::string& origin,
                                   const int&         status,
                                   const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set FITS message
    char err_text[31];
    __ffgerr(status, err_text);
    m_message  = std::string(err_text);
    m_message += " (status=" + gammalib::str(status) + ").";

    // Add optional error message
    if (message.length() > 0) {
        m_message += (" " + message);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief File error
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Error message.
 ***************************************************************************/
GException::file_error::file_error(const std::string& origin,
                                   const std::string& message)
{
    // Set origin
    m_origin  = origin;

    // Set message string
    m_message = "File error.";
    if (message.length() > 0) {
        m_message += (" " + message);
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
GException::feature_not_implemented::feature_not_implemented(const std::string& origin,
                                                             const std::string& message)
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
                 " https://cta-redmine.irap.omp.eu/projects/gammalib/,"
                 " join this error message and provide a detailed"
                 " description of your needs.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks energy interval
 *
 * @param[in] origin Method performing the check.
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * @exception GException::invalid_argument
 *            Minimum energy is equal or larger than maximum energy.
 *
 * Checks that the minimum energy is smaller than the maximum energy.
 ***************************************************************************/
void gammalib::check_energy_interval(const std::string& origin,
                                     const GEnergy&     emin,
                                     const GEnergy&     emax)
{
    // Throw an exception if energy interval is invalid
    if (emin == emax) {
        std::string msg = "Minimum energy "+emin.print()+" is equal to "
                          "maximum energy "+emax.print()+". Please specify "
                          "a minimum energy that is smaller than the maximum "
                          "energy.";
        throw GException::invalid_argument(origin, msg);
    }
    else if (emin > emax) {
        std::string msg = "Minimum energy "+emin.print()+" is larger than "
                          "maximum energy "+emax.print()+". Please specify a "
                          "minimum energy that is smaller than the maximum "
                          "energy.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return
    return;
}
