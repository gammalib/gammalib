/***************************************************************************
 *         GXmlComment.cpp - XML comment node class implementation         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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
//#include <stdio.h>            // fopen, fgets, fclose, etc...
#include <iostream>
#include <cstdio>             // fopen, fgets, fclose, etc...
#include "GException.hpp"
#include "GXmlComment.hpp"

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
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] node Object from which the instance should be built.
 ***************************************************************************/
GXmlComment::GXmlComment(const GXmlComment& node) : GXmlNode(node)
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
 * @param[in] text Text for instance building.
 ***************************************************************************/
GXmlComment::GXmlComment(const std::string& segment) : GXmlNode()
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
 * @param[in] node Object which should be assigned.
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
 * @brief Write comment into file
 *
 * @param[in] fptr File pointer.
 ***************************************************************************/
void GXmlComment::write(FILE* fptr, int indent) const
{
    // Write comment into file
    for (int k = 0; k < indent; ++k)
        fprintf(fptr, " ");
    fprintf(fptr, "<!--%s-->\n", m_comment.c_str());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print comment in output stream
 *
 * @param[in] os Output stream into which the node will be printed.
 ***************************************************************************/
void GXmlComment::print(std::ostream& os, int indent) const
{
    // Put comment into output stream
    for (int k = 0; k < indent; ++k)
        os << " ";
    os << "GXmlComment::" << m_comment;

    // Return
    return;
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
 * @param[in] node Object from which members which should be copied.
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
 * @brief Clone class
***************************************************************************/
GXmlComment* GXmlComment::clone(void) const
{
    return new GXmlComment(*this);
}


/***********************************************************************//**
 * @brief Parse comment segment string
 *
 * @param[in] segement Segment string.
 *
 * @exception GException::xml_syntax_error
 *            XML syntax error.
 *
 * Parse the segment string and extract the comment.
 *
 * @TODO Check validity of characters in comment string
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

        //TODO Check validity of characters comment string

    } // endif: string is not empty

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

