/***************************************************************************
 *          GXmlElement.cpp - XML element node class implementation        *
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
 * @file GXmlElement.cpp
 * @brief XML element node class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <iostream>
#include "GException.hpp"
#include "GXmlElement.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PARSE_START                "GXmlElement::parse_start(std::string&)"
#define G_PARSE_STOP                  "GXmlElement::parse_stop(std::string&)"

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
GXmlElement::GXmlElement(void) : GXmlNode()
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
GXmlElement::GXmlElement(const GXmlElement& node) : GXmlNode(node)
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
 * @param[in] segment XML segment from which instance is built.
 ***************************************************************************/
GXmlElement::GXmlElement(const std::string& segment) : GXmlNode()
{
    // Initialise private members for clean destruction
    init_members();

    // Parse start element
    parse_start(segment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlElement::~GXmlElement(void)
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
GXmlElement& GXmlElement::operator= (const GXmlElement& node)
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
void GXmlElement::clear(void)
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
void GXmlElement::print(std::ostream& os, int indent) const
{
    // Put element name in output stream
    for (int k = 0; k < indent; ++k)
        os << " ";
    os << "GXmlElement::" << m_name;

    // Put children in stream
    for (int i = 0; i < size(); ++i) {
        os << std::endl;
        m_nodes[i]->print(os, indent+4);
    }

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
void GXmlElement::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_parent = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node Object from which members which should be copied.
 ***************************************************************************/
void GXmlElement::copy_members(const GXmlElement& node)
{
    // Copy attributes
    m_name   = node.m_name;
    m_parent = node.m_parent;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXmlElement::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GXmlElement* GXmlElement::clone(void) const
{
    return new GXmlElement(*this);
}


/***********************************************************************//**
 * @brief Parse element start segment string
 *
 * @param[in] segement Segment string.
 *
 * @exception GException::xml_syntax_error
 *            XML syntax error.
 *
 * Parse the segment string and set class members based on the information
 * that is found. The method also performs syntax checking.
 ***************************************************************************/
void GXmlElement::parse_start(const std::string& segment)
{
    // Get length of segment
    int n = segment.length();

    // Check on existence of brackets
    if (n < 2 || (segment.compare(0,1,"<") != 0) ||
                 (segment.compare(n-1,1,">") != 0))
        throw GException::xml_syntax_error(G_PARSE_START, segment,
                                           "tag brackets missing");

    // Extract element name
    size_t pos = segment.find_first_of("\x20\x09\x0d\x0a>", 1);
    if (pos == 1)
        throw GException::xml_syntax_error(G_PARSE_START, segment,
                          "no whitespace allowed after '<'");
    if (pos == std::string::npos)
        throw GException::xml_syntax_error(G_PARSE_START, segment,
                          "element name not found");
    m_name = segment.substr(1, pos-1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Parse element stop segment string
 *
 * @param[in] segement Segment string.
 *
 * @exception GException::xml_syntax_error
 *            XML syntax error.
 *
 * Parse the stop segment string and verify the syntax.
 ***************************************************************************/
void GXmlElement::parse_stop(const std::string& segment)
{
    // Get length of segment
    int n = segment.length();

    // Check on existence of brackets
    if (n < 3 || (segment.compare(0,2,"</") != 0) ||
                 (segment.compare(n-1,1,">") != 0))
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "incorrect or missing tag brackets");

    // Extract and verify element name
    size_t pos = segment.find_first_of("\x20\x09\x0d\x0a>", 2);
    if (pos == 2)
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "no whitespace allowed after '</'");
    if (pos == std::string::npos)
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "element name not found");
    std::string name = segment.substr(2, pos-2);
    if (name != m_name)
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "invalid name in element stop tag"
                          " (found "+name+", expected "+m_name);

    // Verify that no further characters exist in element stop tag
    size_t pos2 = segment.find_first_of("\x20\x09\x0d\x0a>", pos);
    if (pos2 != n-1)
        throw GException::xml_syntax_error(G_PARSE_STOP, segment,
                          "invalid characters found after element name");

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

