/***************************************************************************
 *                   GException.cpp  -  exception handlers                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2011 by Jurgen Knodlseder                           *
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
 * @brief Exception handler interface implementation.
 * @author J. Knodlseder
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
    return tochar(message);
}


/***********************************************************************//**
 * @brief Feature not implement
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::feature_not_implemented::feature_not_implemented(std::string origin,
                                                             std::string message)
{
    m_origin   = origin;
    if (message.length() > 0)
        m_message = message;
    else
        m_message = "Feature not implemented.";
    m_message += " In case that you need this feature for your application"
                 " please submit a feature request on"
                 " https://sourceforge.net/projects/gammalib/,"
                 " join this error message and provide a detailed"
                 " description of your needs.";
    return;
}


/***********************************************************************//**
 * @brief Invalid argument
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::invalid_argument::invalid_argument(std::string origin,
                                               std::string message)
{
    // Set origin
    m_origin   = origin;

    // Set message string
    m_message = "Invalid argument encountered.";
    if (message.length() > 0)
        m_message += " " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bad type
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional message.
 ***************************************************************************/
GException::bad_type::bad_type(std::string origin, std::string message)
{
    m_origin   = origin;
    m_message = "Invalid type conversion. " + message;
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
    m_origin  = origin;
    m_message = "Memory allocation error (" + str((int)num) + " elements)";
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
    m_origin  = origin;
    m_message = "Not enough nodes in node array (" + str(num) + " nodes).";
    return;
}


/***********************************************************************//**
 * @brief File not found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Filename.
 ***************************************************************************/
GException::file_not_found::file_not_found(std::string origin, std::string filename)
{
    m_origin  = origin;
    m_message = "File \"" + filename +"\" not found.";
    return;
}


/***********************************************************************//**
 * @brief File open error
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Filename.
 ***************************************************************************/
GException::file_open_error::file_open_error(std::string origin, std::string filename)
{
    m_origin  = origin;
    m_message = "Unable to open file \"" + filename +"\"";
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
    m_message = "Inconsistent number of columns found in row "+str(row)+
                " of file \""+filename+"\" ("+str(cols)+" expected, "+
                str(elements)+" encountered). "+message;
    return;
}
