/***************************************************************************
 *         GXmlAttribute.cpp - XML attribute class implementation          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GXmlAttribute.cpp
 * @brief XML attribute class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include "GException.hpp"
#include "GTools.hpp"
#include "GXmlAttribute.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GXmlAttribute::GXmlAttribute(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] attr Element attribute.
 ***************************************************************************/
GXmlAttribute::GXmlAttribute(const GXmlAttribute& attr)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(attr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Name-Value pair constructor
 *
 * @param[in] name Attribute name.
 * @param[in] value Attribute value.
 *
 * Construct attribute form a @p name and a @p value. Predefined entities
 * (e.g. &quot;) in attribute values are automatically converted into normal
 * characters. The constructor strips any existing " or ' leading and
 * trailing hyphens from the @p value string.
 ***************************************************************************/
GXmlAttribute::GXmlAttribute(const std::string& name, const std::string& value)
{
    // Initialise members
    init_members();

    // Create working copy of attribute value
    std::string v(value);

    // Strip any pair of leading and trailing hyphens
    int n = v.length();
    if (n >= 2) {
        if (((v[0] ==  '"') && (v[n-1] ==  '"')) ||
            ((v[0] == '\'') && (v[n-1] == '\''))) {
            if (n > 2) {
                v = v.substr(1, n-2);
            }
            else {
                v = "";
            }
        }
    }

    // Set attribute
    this->name(name);
    this->value(gammalib::xml2str(v));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlAttribute::~GXmlAttribute(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] attr Element attribute.
 * @return Element attribute.
 ***************************************************************************/
GXmlAttribute& GXmlAttribute::operator=(const GXmlAttribute& attr)
{
    // Execute only if object is not identical
    if (this != &attr) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(attr);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

 /***********************************************************************//**
 * @brief Clear element attribute
 *
 * Resets element attribute to a clean initial state.
 ***************************************************************************/
void GXmlAttribute::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone element attribute
 *
 * @return Pointer to deep copy of an element attribute.
 ***************************************************************************/
GXmlAttribute* GXmlAttribute::clone(void) const
{
    // Clone XML attribute
    return new GXmlAttribute(*this);
}


/***********************************************************************//**
 * @brief Write attribute into URL
 *
 * @param[in] url Unified Resource Locator.
 *
 * Writes the element attribute into the @p url. Special characters are
 * automatically transformed into predefined entities (e.g. &quot;).
 ***************************************************************************/
void GXmlAttribute::write(GUrl& url) const
{
    // Convert value into XML format and add hyphens
    std::string value = "\""+gammalib::str2xml(m_value)+"\"";

    // Write attribute into URL
    url.printf(" %s=%s", m_name.c_str(), value.c_str());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print element attribute
 *
 * @param[in] chatter Chattiness.
 * @return String containing element attribute.
 ***************************************************************************/
std::string GXmlAttribute::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;
    
    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append attribute to string
        result.append(" "+m_name+"=\""+m_value+"\"");

    } // endif: chatter was not silent

    // Return
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GXmlAttribute::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_value.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] attr Element attribute.
 ***************************************************************************/
void GXmlAttribute::copy_members(const GXmlAttribute& attr)
{
    // Copy attributes
    m_name  = attr.m_name;
    m_value = attr.m_value;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXmlAttribute::free_members(void)
{
    // Return
    return;
}
