/***************************************************************************
 *                     GXml.cpp - XML class implementation                 *
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
 * @file GXml.cpp
 * @brief XML class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include <string>
#include <iostream>
#include "GXml.hpp"
#include "GXmlNode.hpp"
#include "GXmlDocument.hpp"
#include "GXmlText.hpp"
#include "GXmlElement.hpp"
#include "GXmlComment.hpp"
#include "GXmlPI.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                                     "GXml::load(std::string&)"
#define G_SAVE                                     "GXml::load(std::string&)"
#define G_PARSE                                          "GXml::parse(FILE*)"
#define G_PROCESS              "GXml::process(GXmlNode*, const std::string&)"

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
GXml::GXml(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] xml XML object.
 ***************************************************************************/
GXml::GXml(const GXml& xml)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename XML file name.
 ***************************************************************************/
GXml::GXml(const std::string& filename)
{
    // Initialise members
    init_members();

    // Load XML file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXml::~GXml(void)
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
 * @param[in] xml XML object.
 ***************************************************************************/
GXml& GXml::operator= (const GXml& xml)
{
    // Execute only if object is not identical
    if (this != &xml) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(xml);

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
 * Reset object to a clean initial state.
 ***************************************************************************/
void GXml::clear(void)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append child node to XML document.
 *
 * @param[in] node Child node.
 ***************************************************************************/
void GXml::append(GXmlNode* node)
{
    // Append node
    m_root.append(node);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load XML file.
 *
 * @param[in] filename Name of file to be loaded.
 *
 * @exception GException::file_not_found
 *            XML file not found.
 * @exception GException::file_open_error
 *            Unable to open XML file (read access requested).
 *
 * Loads XML file by reading all lines from the XML file.
 ***************************************************************************/
void GXml::load(const std::string& filename)
{
    // Clear object
    clear();

    // Check if file exists
    if (!file_exists(filename))
        throw GException::file_not_found(G_LOAD, filename);

    // Open parameter file
    FILE* fptr = fopen(filename.c_str(), "r");
    if (fptr == NULL)
        throw GException::file_open_error(G_LOAD, filename);

    // Parse file
    parse(fptr);

    // Close file
    fclose(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save XML file.
 *
 * @param[in] filename Name of file to be saved.
 *
 * @exception GException::file_open_error
 *            Unable to open XML file (write access requested).
 *
 * Save XML document into file.
 ***************************************************************************/
void GXml::save(const std::string& filename)
{
    // Open parameter file
    FILE* fptr = fopen(filename.c_str(), "w");
    if (fptr == NULL)
        throw GException::file_open_error(G_SAVE, filename);

    // Write XML document
    m_root.write(fptr);

    // Close file
    fclose(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of children in document root
 ***************************************************************************/
int GXml::children(void) const
{
    // Return number of document root children
    return m_root.children();
}


/***********************************************************************//**
 * @brief Return pointer on document child node
 *
 * @param[in] index Index of node (0,1,2,...)
 ***************************************************************************/
GXmlNode* GXml::child(int index) const
{
    // Return pointer
    return m_root.child(index);
}


/***********************************************************************//**
 * @brief Return number of child elements in document root
 ***************************************************************************/
int GXml::elements(void) const
{
    // Return number
    return m_root.elements();
}


/***********************************************************************//**
 * @brief Return number of child elements of given name in document root
 *
 * @param[in] name Name of child elements.
 ***************************************************************************/
int GXml::elements(const std::string& name) const
{
    // Return number
    return m_root.elements(name);
}


/***********************************************************************//**
 * @brief Return pointer on child element
 *
 * @param[in] index Index of child element (0,1,2,...)
 ***************************************************************************/
GXmlElement* GXml::element(int index) const
{
    // Return pointer
    return (GXmlElement*)m_root.element(index);
}


/***********************************************************************//**
 * @brief Return pointer on child element of a given name
 *
 * @param[in] name Name of child elements.
 * @param[in] index Index of child element (0,1,2,...)
 ***************************************************************************/
GXmlElement* GXml::element(const std::string& name, int index) const
{
    // Return pointer
    return (GXmlElement*)m_root.element(name, index);
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GXml::init_members(void)
{
    // Initialise members
    m_root.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] xml Object from which members which should be copied.
 ***************************************************************************/
void GXml::copy_members(const GXml& xml)
{
    // Copy attributes
    m_root = xml.m_root;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXml::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Parse XML file
 *
 * @param[in] fptr Pointer to file to be parsed.
 *
 * @exception GException::xml_syntax_error
 *            XML syntax error.
 *
 * Parse the XML file and add nodes corresponding to the content to the
 * object. The XML file is split into segments, made either of text or of
 * tags.
 ***************************************************************************/
void GXml::parse(FILE* fptr)
{
    // Initialise parser
    int         c;
    bool        in_markup  = false;
    std::string segment;
    GXmlNode*   current = &m_root;

    // Main parsing loop
    while ((c = fgetc(fptr)) != EOF) {

        // Convert special characters into line feeds
        if (c == '\x85' || c == L'\x2028') {
            if (in_markup)
                 throw GException::xml_syntax_error(G_PARSE, segment,
                                   "invalid character encountered");
            else
                c = '\x0a';
        }

        // If we are not within a markup and if a markup is reached then add
        // the text segment to the nodes and switch to in_markup mode
        if (in_markup == false) {

            // Markup start reached?
            if (c == '<') {

                // Add text segment to nodes (ignores empty segments)
                process_text(&current, segment);

                // Prepare new segment and signal that we are within tag
                segment.clear();
                segment.append(1, (char)c);
                in_markup = true;
            }

            // Markup stop encountered?
            else if (c == '>') {
                 segment.append(1, (char)c);
                 throw GException::xml_syntax_error(G_PARSE, segment,
                       "unexpected closing bracket \">\" encountered");
            }

            // ... otherwise add character to segment
            else
                segment.append(1, (char)c);
        }

        // If we are within a markuo and if a markup end is reached then process
        // the markup and switch to not in_tag mode
        else {

            // Markup stop reached?
            if (c == '>') {

                // Process markup
                segment.append(1, (char)c);
                process_markup(&current, segment);

                // Prepare new segment and signal that we are not within markup
                segment.clear();
                in_markup = false;
            }

            // Markup start encountered?
            else if (c == '<') {
                 segment.append(1, (char)c);
                 throw GException::xml_syntax_error(G_PARSE, segment,
                       "unexpected opening bracket \"<\" encountered");
            }

            // ... otherwise add character to segment
            else
                segment.append(1, (char)c);
        }

    }

    // Process any pending segment
    if (segment.size() > 0) {
        if (in_markup)
            process_markup(&current, segment);
        else
            process_text(&current, segment);
    } // endif: there was a pending segment

    // Return
    return;
}


/***********************************************************************//**
 * @brief Process markup segment
 *
 * @param[in] current Handle to current node.
 * @param[in] segment Segment string.
 *
 * Process markup segment.
 ***************************************************************************/
void GXml::process_markup(GXmlNode** current, const std::string& segment)
{
    // Determine segment tag type
    MarkupType type = get_markuptype(segment);

    // Do tag specific processing
    switch (type) {

    // Handle element start tag
    case MT_ELEMENT_START:
        {
            // Create new element node, set it's parent, append it to the current
            // node and make it the current node
            GXmlElement* element = new GXmlElement(segment);
            element->parent(*current);
            (*current)->append(element);
            (*current) = element;
        }
        break;

    // Handle element end tag
    case MT_ELEMENT_END:
        {
            // Check if we expect an element end tag
            if ((*current)->type() != GXmlNode::NT_ELEMENT)
                throw GException::xml_syntax_error(G_PROCESS, segment,
                      "unexpected element end tag encountered");

            // Check if we have the correct end tag
            GXmlElement* element = (GXmlElement*)(*current);
            element->parse_stop(segment);

            // Set current node pointer back to parent of the current node
            (*current) = element->parent();
        }
        break;

    // Append empty-element tag
    case MT_ELEMENT_EMPTY:
        {
            GXmlElement* element = new GXmlElement(segment);
            element->parent(*current);
            (*current)->append(element);
        }
        break;

    // Append comment markup
    case MT_COMMENT:
        {
            GXmlComment* comment = new GXmlComment(segment);
            (*current)->append(comment);
        }
        break;

    // Declaration markup
    case MT_DECLARATION:
        {
            // Verify if declaration tag is allowed
            if (*current != &m_root)
                throw GException::xml_syntax_error(G_PROCESS, segment,
                      "unexpected declaration markup encountered");
            if (m_root.children() > 0)
                throw GException::xml_syntax_error(G_PROCESS, segment,
                      "declaration markup only allowed in first line");

            // Create temporary element to read in declaration attributes
            GXmlElement* element = new GXmlElement(segment);
            size_t       pos     = 5;
            while (pos != std::string::npos)
                element->parse_attribute(&pos, segment);

            // Set attribute values
            std::string version    = element->attribute("version");
            std::string encoding   = element->attribute("encoding");
            std::string standalone = element->attribute("standalone");
            if (version.length() > 0)
                m_root.version(version);
            if (encoding.length() > 0)
                m_root.encoding(encoding);
            if (standalone.length() > 0)
                m_root.standalone(standalone);

            // Delete temporary element
            delete element;
        }
        break;

    // Processing tag
    case MT_PROCESSING:
        {
            GXmlPI* pi = new GXmlPI(segment);
            (*current)->append(pi);
        }
        break;

    // Invalid tag, throw an error
    case MT_INVALID:
        throw GException::xml_syntax_error(G_PROCESS, segment, "invalid tag");
        break;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Process text segment
 *
 * @param[in] current Handle to current node.
 * @param[in] segment Segment string.
 *
 * Process text segment.
 ***************************************************************************/
void GXml::process_text(GXmlNode** current, const std::string& segment)
{
    // Continue only if text segment is not empty
    if (segment.size() > 0) {

        // Continue only if non whitespace characters are found
        size_t pos = segment.find_first_not_of("\x20\x09\x0d\x0a\x85");
        if (pos != std::string::npos) {

            // Allocate and append node
            GXmlText* node = new GXmlText(segment);
            (*current)->append(node);

        } // endif: there was not only whitespace

    } // endif: segment was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get Markup type of segment
 *
 * @param[in] segment Segment for which Markup Type should be determined.
 *
 * Returns Markup Type of segment.
 ***************************************************************************/
GXml::MarkupType GXml::get_markuptype(const std::string& segment) const
{
    // Initialise with invalid Markup Type
    MarkupType type = MT_INVALID;

    // Get length of segment
    int n = segment.length();

    // Check for comment
    if (n >= 7 && (segment.compare(0,4,"<!--") == 0) &&
             (segment.compare(n-3,3,"-->") == 0))
        type = MT_COMMENT;

    // Check for declaration
    else if (n >= 7 && (segment.compare(0,6,"<?xml ") == 0) &&
             (segment.compare(n-2,2,"?>") == 0))
        type = MT_DECLARATION;

    // Check for processing instruction
    else if (n >= 4 && (segment.compare(0,2,"<?") == 0) &&
             (segment.compare(n-2,2,"?>") == 0))
        type = MT_PROCESSING;

    // Check for empty element tag
    else if (n >= 3 && (segment.compare(0,1,"<") == 0) &&
             (segment.compare(n-2,2,"/>") == 0))
        type = MT_ELEMENT_EMPTY;

    // Check for element end tag
    else if (n >= 3 && (segment.compare(0,2,"</") == 0) &&
             (segment.compare(n-1,1,">") == 0))
        type = MT_ELEMENT_END;

    // Check for element start tag
    else if (n >= 2 && (segment.compare(0,1,"<") == 0) &&
             (segment.compare(n-1,1,">") == 0))
        type = MT_ELEMENT_START;

    // Return type
    return type;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] xml Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GXml& xml)
{
    // Put object in stream
    os << "=== GXml ===" << std::endl;
    xml.m_root.print(os, 0);

    // Return output stream
    return os;
}

