/***************************************************************************
 *               GXmlPI.cpp - XML PI node class implementation             *
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
 * @file GXmlPI.cpp
 * @brief XML PI node class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GXmlPI.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PARSE                                 "GXmlPI::parse(std::string&)"

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
GXmlPI::GXmlPI(void) : GXmlNode()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] node XML Processing Instruction.
 ***************************************************************************/
GXmlPI::GXmlPI(const GXmlPI& node) : GXmlNode(node)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(node);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Segment constructor
 *
 * @param[in] segment Text for instance building.
 ***************************************************************************/
GXmlPI::GXmlPI(const std::string& segment) : GXmlNode()
{
    // Initialise private members for clean destruction
    init_members();

    // Parse segment
    parse(segment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlPI::~GXmlPI(void)
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
 * @param[in] node XML Processing Instruction.
 * @return XML Processing Instruction.
 ***************************************************************************/
GXmlPI& GXmlPI::operator=(const GXmlPI& node)
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
 * @brief Clear XML Processing Instruction
 *
 * Resets the XML Processing Instruction to an initial state.
 ***************************************************************************/
void GXmlPI::clear(void)
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
 * @brief Clone XML Processing Instruction
 *
 * @return Pointer to deep copy of XML Processing Instruction.
 ***************************************************************************/
GXmlPI* GXmlPI::clone(void) const
{
    // Clone XML PI
    return new GXmlPI(*this);
}


/***********************************************************************//**
 * @brief Write Processing Instruction into URL
 *
 * @param[in] url Unified Resource Locator.
 * @param[in] indent Text indentation (default = 0).
 ***************************************************************************/
void GXmlPI::write(GUrl& url, const int& indent) const
{
    // Prepend indentation
    for (int k = 0; k < indent; ++k) {
        url.printf(" ");
    }

    // Write Processing Instruction into URL
    url.printf("<?%s?>\n", m_pi.c_str());

    // Return
    return;
}



/***********************************************************************//**
 * @brief Print XML Processing Instruction
 *
 * @param[in] indent Text indentation (default = 0).
 * @return String containing XML Processing Instruction.
 ***************************************************************************/
std::string GXmlPI::print(const int& indent) const
{
    // Initialise result string
    std::string result = fill(" ", indent);

    // Append comment to string
    result.append("GXmlPI::"+m_pi);

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
void GXmlPI::init_members(void)
{
    // Initialise members
    m_pi.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node Object from which members which should be copied.
 ***************************************************************************/
void GXmlPI::copy_members(const GXmlPI& node)
{
    // Copy attributes
    m_pi = node.m_pi;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXmlPI::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Parse comment segment string
 *
 * @param[in] segment Segment string.
 *
 * @exception GException::xml_syntax_error
 *            XML syntax error.
 *
 * Parse the segment string.
 ***************************************************************************/
void GXmlPI::parse(const std::string& segment)
{
    // Get length of segment
    int n = segment.length();

    // Check on existence of brackets
    if (n < 4 || (segment.compare(0,2,"<?")   != 0) ||
                 (segment.compare(n-2,2,"?>") != 0)) {
        throw GException::xml_syntax_error(G_PARSE, segment,
              "missing or invalid Processing Instruction brackets");
    }

    // Set comment
    if (n > 4) {
        m_pi = segment.substr(2, n-4);
    }
    else {
        m_pi.clear();
    }

    // Return
    return;
}
