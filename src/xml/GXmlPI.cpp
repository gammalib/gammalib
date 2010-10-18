/***************************************************************************
 *               GXmlPI.cpp - XML PI node class implementation             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GXmlPI.cpp
 * @brief XML PI node class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GXmlPI.hpp"

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
 * @param[in] node Object from which the instance should be built.
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
 * @param[in] text Text for instance building.
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
 * @param[in] node Object which should be assigned.
 ***************************************************************************/
GXmlPI& GXmlPI::operator= (const GXmlPI& node)
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
 * @brief Print node in output stream
 *
 * @param[in] os Output stream into which the node will be printed.
 ***************************************************************************/
void GXmlPI::print(std::ostream& os, int indent) const
{
    // Put element name in output stream
    for (int k = 0; k < indent; ++k)
        os << " ";
    os << "GXmlPI::" << m_pi;

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
 * @brief Clone class
***************************************************************************/
GXmlPI* GXmlPI::clone(void) const
{
    return new GXmlPI(*this);
}


/***********************************************************************//**
 * @brief Parse comment segment string
 *
 * @param[in] segement Segment string.
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
    if (n < 4 || (segment.compare(0,2,"<?") != 0) &&
                 (segment.compare(n-2,2,"?>") != 0))
        throw GException::xml_syntax_error(G_PARSE, segment,
              "missing or invalid Processing Instruction brackets");

    // Set comment
    if (n > 4)
        m_pi = segment.substr(2, n-4);
    else
        m_pi.clear();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

