/***************************************************************************
 *        GXmlDocument.cpp - XML document node class implementation        *
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
 * @file GXmlDocument.cpp
 * @brief XML document node class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>             // std::fprintf
#include <iostream>
#include "GXmlDocument.hpp"

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
 * @param[in] node Object from which the instance should be built.
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
 * @param[in] node Object which should be assigned.
 ***************************************************************************/
GXmlDocument& GXmlDocument::operator= (const GXmlDocument& node)
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
 * @brief Clone class
***************************************************************************/
GXmlDocument* GXmlDocument::clone(void) const
{
    return new GXmlDocument(*this);
}


/***********************************************************************//**
 * @brief Write node into file
 *
 * @param[in] fptr File pointer.
 * @param[in] indent Text indentation.
 ***************************************************************************/
void GXmlDocument::write(FILE* fptr, int indent) const
{
    // Write document header in file
    std::fprintf(fptr, "<?xml version=\"%s\" encoding=\"%s\" standalone=\"%s\"?>\n",
                 version().c_str(),
                 encoding().c_str(),
                 standalone().c_str());

    // Write children in file
    for (int i = 0; i < children(); ++i)
        m_nodes[i]->write(fptr, indent);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print node in output stream
 *
 * @param[in] os Output stream.
 * @param[in] indent Text indentation.
 ***************************************************************************/
void GXmlDocument::print(std::ostream& os, int indent) const
{
    // Put document header in output stream
    for (int k = 0; k < indent; ++k)
        os << " ";
    os << "GXmlDocument::";
    os << "version=" << version() << " ";
    os << "encoding=" << encoding() << " ";
    os << "standalone=" << standalone() << " ";

    // Put children in stream
    for (int i = 0; i < children(); ++i) {
        os << std::endl;
        m_nodes[i]->print(os, indent);
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


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
