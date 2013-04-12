/***************************************************************************
 *        GXmlDocument.cpp - XML document node class implementation        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GXmlDocument.cpp
 * @brief XML document node class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GXmlDocument.hpp"
#include "GTools.hpp"

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
GXmlDocument::GXmlDocument(void) : GXmlNode()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] node XML document.
 ***************************************************************************/
GXmlDocument::GXmlDocument(const GXmlDocument& node) : GXmlNode(node)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(node);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlDocument::~GXmlDocument(void)
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
 * @param[in] node XML document.
 * @return XML document.
 ***************************************************************************/
GXmlDocument& GXmlDocument::operator=(const GXmlDocument& node)
{
    // Execute only if object is not identical
    if (this != &node) {

        // Copy base class members
        this->GXmlNode::operator=(node);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(node);

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
 * @brief Clear XML document
 *
 * Resets the XML document to a clean initial state.
 ***************************************************************************/
void GXmlDocument::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GXmlNode::free_members();

    // Initialise members
    this->GXmlNode::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone XML document
 *
 * @return Pointer to deep copy of XML document.
 ***************************************************************************/
GXmlDocument* GXmlDocument::clone(void) const
{
    // Clone document
    return new GXmlDocument(*this);
}


/***********************************************************************//**
 * @brief Write XML document into URL
 *
 * @param[in] url Unified Resource Locator.
 * @param[in] indent Text indentation (default = 0).
 *
 * Writes the XML document into a @p url object.
 ***************************************************************************/
void GXmlDocument::write(GUrl& url, const int& indent) const
{
    // Write document header into URL
    url.printf("<?xml version=\"%s\" encoding=\"%s\" standalone=\"%s\"?>\n",
               version().c_str(),
               encoding().c_str(),
               standalone().c_str());

    // Write children into URL
    for (int i = 0; i < m_nodes.size(); ++i) {
        m_nodes[i]->write(url, indent);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print XML document
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @param[in] indent Text indentation (default to 0).
 * @return String containing XML document.
 ***************************************************************************/
std::string GXmlDocument::print(const GChatter& chatter,
                                const int&      indent) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Initialise result string
        result = gammalib::fill(" ", indent);

        // Append document to string
        result.append("GXmlDocument::");
        result.append("version=" + version());
        result.append(" encoding=" + encoding());
        result.append(" standalone=" + standalone());

        // Append children
        for (int i = 0; i < m_nodes.size(); ++i) {
            result.append("\n" + m_nodes[i]->print(chatter, indent));
        }

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
void GXmlDocument::init_members(void)
{
    // Initialise members
    m_version.clear();
    m_encoding.clear();
    m_standalone.clear();
    m_version.name("version");
    m_encoding.name("encoding");
    m_standalone.name("m_standalone");
    m_version.value("\"1.0\"");
    m_encoding.value("\"UTF-8\"");
    m_standalone.value("\"no\"");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node Object from which members which should be copied.
 ***************************************************************************/
void GXmlDocument::copy_members(const GXmlDocument& node)
{
    // Copy attributes
    m_version    = node.m_version;
    m_encoding   = node.m_encoding;
    m_standalone = node.m_standalone;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXmlDocument::free_members(void)
{
    // Return
    return;
}
