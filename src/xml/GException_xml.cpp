/***************************************************************************
 *               GException_xml.cpp  -  XML exception handlers             *
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
 * @file GException_xml.cpp
 * @brief Implement exceptions for the XML module
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief XML syntax error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] segment XML segment for which syntax error was encountered.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::xml_syntax_error::xml_syntax_error(std::string origin,
                                               std::string segment,
                                               std::string message)
{
    // Set origin and message
    m_origin  = origin;
    
    // Set message
    m_message = "XML syntax error ("+message+")";

    // Set segment
    if (segment.length() > 0) {
        m_message += " occured in " + segment;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid XML attribute value
 *
 * @param[in] origin Method that throws the error.
 * @param[in] value XML attribute value.
 ***************************************************************************/
GException::xml_attribute_value::xml_attribute_value(std::string origin,
                                                     std::string value)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid XML attribute value: "+value;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid XML node type
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type XML node type.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::xml_bad_node_type::xml_bad_node_type(std::string origin,
                                                 std::string type,
                                                 std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid XML node type ("+type+").";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML element name not found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] name XML element name.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::xml_name_not_found::xml_name_not_found(std::string origin,
                                                   std::string name,
                                                   std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "XML element name \""+name+"\" not found.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid number of parameters in XML element
 *
 * @param[in] origin Method that throws the error.
 * @param[in] xml XML element.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::xml_invalid_parnum::xml_invalid_parnum(std::string origin,
                                                   GXmlElement xml,
                                                   std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid number of parameters found in XML element.";
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid parameter names in XML element
 *
 * @param[in] origin Method that throws the error.
 * @param[in] xml XML element.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::xml_invalid_parnames::xml_invalid_parnames(std::string origin,
                                                       GXmlElement xml,
                                                       std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid parameter names found in XML element." ;
    if (message.length() > 0) {
        m_message += " " + message;
    }

    // Return
    return;
}


