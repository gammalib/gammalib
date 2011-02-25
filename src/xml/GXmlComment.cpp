/***************************************************************************
 *         GXmlComment.cpp - XML comment node class implementation         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GXmlComment.cpp
 * @brief XML comment node class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GXmlComment.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PARSE                            "GXmlComment::parse(std::string&)"

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
GXmlComment::GXmlComment(void) : GXmlNode()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] node XML comment.
 ***************************************************************************/
GXmlComment::GXmlComment(const GXmlComment& node) : GXmlNode(node)
{
    // Initialise members
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
GXmlComment::GXmlComment(const std::string& segment) : GXmlNode()
{
    // Initialise members
    init_members();

    // Parse segment
    parse(segment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlComment::~GXmlComment(void)
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
 * @param[in] node XML comment.
 ***************************************************************************/
GXmlComment& GXmlComment::operator= (const GXmlComment& node)
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
 * @brief Clear object.
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GXmlComment::clear(void)
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
 * @brief Clone class
***************************************************************************/
GXmlComment* GXmlComment::clone(void) const
{
    return new GXmlComment(*this);
}


/***********************************************************************//**
 * @brief Write comment into file
 *
 * @param[in] fptr File pointer.
 * @param[in] indent Text indentation.
 ***************************************************************************/
void GXmlComment::write(FILE* fptr, int indent) const
{
    // Write comment into file
    for (int k = 0; k < indent; ++k)
        std::fprintf(fptr, " ");
    std::fprintf(fptr, "<!--%s-->\n", m_comment.c_str());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print comment in output stream
 *
 * @param[in] indent Text indentation.
 ***************************************************************************/
std::string GXmlComment::print(int indent) const
{
    // Initialise result string
    std::string result = fill(" ", indent);

    // Append comment to string
    result.append("GXmlComment::"+m_comment);

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
void GXmlComment::init_members(void)
{
    // Initialise members
    m_comment.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node XML comment.
 ***************************************************************************/
void GXmlComment::copy_members(const GXmlComment& node)
{
    // Copy attributes
    m_comment = node.m_comment;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXmlComment::free_members(void)
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
 * Parse the segment string and extract the comment.
 *
 * @todo Check validity of characters in comment string
 ***************************************************************************/
void GXmlComment::parse(const std::string& segment)
{
    // Initialise comment string
    m_comment.clear();

    // Get length of segment
    int n = segment.length();

    // Do nothing if string is empty
    if (n > 0) {

        // If string starts with brackets then check that the brackets are
        // valid comment brackets
        if (segment[0] == '<') {
            if (n < 7 || (segment.compare(0,4,"<!--") != 0) &&
                          (segment.compare(n-3,3,"-->") != 0))
                throw GException::xml_syntax_error(G_PARSE, segment,
                                  "invalid comment brackets");
            else
                m_comment = segment.substr(4, n-7);
        }
        else
            m_comment = segment;

        //@todo Check validity of characters comment string

    } // endif: string is not empty

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

