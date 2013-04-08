/***************************************************************************
 *                  GXmlText.cpp - XML text node class                     *
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
 * @file GXmlText.cpp
 * @brief XML text node class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GXmlText.hpp"
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
GXmlText::GXmlText(void) : GXmlNode()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] node XML text.
 ***************************************************************************/
GXmlText::GXmlText(const GXmlText& node) : GXmlNode(node)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(node);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Text constructor
 *
 * @param[in] text Text string.
 *
 * Construct object by attributing a text string.
 ***************************************************************************/
GXmlText::GXmlText(const std::string& text) : GXmlNode()
{
    // Initialise members
    init_members();

    // Set text
    m_text = text;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlText::~GXmlText(void)
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
 * @param[in] node XML text.
 * @return XML text.
 ***************************************************************************/
GXmlText& GXmlText::operator=(const GXmlText& node)
{
    // Execute only if object is not identical
    if (this != &node) {

        // Copy base class members
        this->GXmlNode::operator=(node);

        // Free members
        free_members();

        // Initialise members
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
 * @brief Clear XML text
 *
 * Resets the XML text to a clean initial state.
 ***************************************************************************/
void GXmlText::clear(void)
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
 * @brief Clone XML text
 *
 * @return Pointer to deep copy of XML text.
 ***************************************************************************/
GXmlText* GXmlText::clone(void) const
{
    // Clone XML text
    return new GXmlText(*this);
}


/***********************************************************************//**
 * @brief Write XML text into URL
 *
 * @param[in] url Unified Resource Locator.
 * @param[in] indent Text indentation (parameter ignored).
 *
 * Writes the text into the URL.
 ***************************************************************************/
void GXmlText::write(GUrl& url, const int& indent) const
{
    // Write text
    url.printf("%s", m_text.c_str());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print XML text
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @param[in] indent Text indentation (default to 0).
 * @return String containing XML text.
 ***************************************************************************/
std::string GXmlText::print(const GChatter& chatter,
                            const int&      indent) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Initialise result string
        result = fill(" ", indent);

        // Append text to string
        result.append("GXmlText::"+m_text);

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
void GXmlText::init_members(void)
{
    // Initialise members
    m_text.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node XML text.
 ***************************************************************************/
void GXmlText::copy_members(const GXmlText& node)
{
    // Copy attributes
    m_text = node.m_text;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXmlText::free_members(void)
{
    // Return
    return;
}

