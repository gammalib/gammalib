/***************************************************************************
 *                     GXml.cpp - XML class implementation                 *
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
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                                     "GXml::load(std::string&)"
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
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] xml Object from which the instance should be built.
 ***************************************************************************/
GXml::GXml(const GXml& xml)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename XML file from which object should be constructed.
 ***************************************************************************/
GXml::GXml(const std::string& filename)
{
    // Initialise private members for clean destruction
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
 * @param[in] xml Object which should be assigned.
 ***************************************************************************/
GXml& GXml::operator= (const GXml& xml)
{
    // Execute only if object is not identical
    if (this != &xml) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
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
 * @brief Load XML file.
 *
 * @param[in] filename Name of file to be loaded.
 *
 * @exception GException::file_open_error
 *            Unable to open parameter file (read access requested).
 *
 * Loads XML file by reading all lines from the XML file.
 ***************************************************************************/
void GXml::load(const std::string& filename)
{
    // Clear object
    clear();

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
 * Save XML file.
 ***************************************************************************/
void GXml::save(const std::string& filename)
{
    
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
    // Free memory
    //if (m_par      != NULL) delete [] m_par;

    // Signal free pointers
    //m_par      = NULL;

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
    bool        in_tag  = false;
    std::string segment;
    GXmlNode*   current = &m_root;

    // Main parsing loop
    while ((c = fgetc(fptr)) != EOF) {

        // If we are not within a tag and if a tag is reached then add
        // the text segment to the nodes and switch to in_tag mode
        if (in_tag == false) {

            // Tag start reached?
            if (c == '<') {

                // Add text segment to nodes (ignores empty segments)
                process_text(&current, segment);

                // Prepare new segment and signal that we are within tag
                segment.clear();
                segment.append(1, (char)c);
                in_tag = true;
            }

            // Tag stop encountered?
            else if (c == '>') {
                 segment.append(1, (char)c);
                 throw GException::xml_syntax_error(G_PARSE, segment,
                                   "unexpected closing bracket '>' encountered");
            }

            // ... otherwise add character to segment
            else
                segment.append(1, (char)c);
        }

        // If we are within a tag and if a tag end is reached then process
        // the tag and switch to not in_tag mode
        else {

            // Tag stop reached?
            if (c == '>') {

                // Process tag
                segment.append(1, (char)c);
                process_tag(&current, segment);

                // Prepare new segment and signal that we are not within tag
                segment.clear();
                in_tag = false;
            }

            // Tag start encountered?
            else if (c == '<') {
                 segment.append(1, (char)c);
                 throw GException::xml_syntax_error(G_PARSE, segment,
                                   "unexpected opening bracket '<' encountered");
            }

            // ... otherwise add character to segment
            else
                segment.append(1, (char)c);
        }

    }

    // Process any pending segment
    if (segment.size() > 0) {
        if (in_tag)
            process_tag(&current, segment);
        else
            process_text(&current, segment);
    } // endif: there was a pending segment

    // Return
    return;
}


/***********************************************************************//**
 * @brief Process tag segment
 *
 * @param[in] current Handle to current node.
 * @param[in] segment Segment string.
 *
 * Process a tag segment.
 ***************************************************************************/
void GXml::process_tag(GXmlNode** current, const std::string& segment)
{
    // Initialise some variables
    GXmlElement* element;

    // Determine segment tag type
    TagType type = get_tagtype(segment);

    // Do tag specific processing
    switch (type) {

    // Element start tag
    case TT_ELEMENT_START:
        // Create new element node, set it's parent, append it to the current
        // node and make it the current node
        element = new GXmlElement(segment);
        element->parent(*current);
        (*current)->append(element);
        (*current) = element;
        std::cout << "START:" << segment << ":" << std::endl;
        break;

    // Element end tag
    case TT_ELEMENT_END:
        // Check if we expect an element end tag
        if ((*current)->type() != GXmlNode::NT_ELEMENT)
            throw GException::xml_syntax_error(G_PROCESS, segment,
                              "unexpected element end tag");

        // Check if we have the correct end tag
        element = (GXmlElement*)(*current);
        element->parse_stop(segment);
        
        // Set current node pointer back to parent of the current node
        (*current) = element->parent();
        std::cout << "STOP:" << segment << ":" << std::endl;
        break;

    // Empty-element tag
    case TT_ELEMENT_EMPTY:
        // Create new element node, set it's parent, and append it to the
        // current node
        element = new GXmlElement(segment);
        element->parent(*current);
        (*current)->append(element);
        std::cout << "EMPTY:" << segment << ":" << std::endl;
        break;

    // Comment tag
    case TT_COMMENT:
        std::cout << "COMMENT:" << segment << ":" << std::endl;
        // Just add comment
        break;

    // Declaration tag
    case TT_DECLARATION:
        std::cout << "DECL:" << segment << ":" << std::endl;
        // Only valid if we are in root node on first element (must be 1st line)
        break;

    // Processing tag
    case TT_PROCESSING:
        std::cout << "PI:" << segment << ":" << std::endl;
        // Just add unknown
        break;

    // Invalid tag, throw an error
    case TT_INVALID:
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
 * @brief Get Tag Type of segment
 *
 * @param[in] segement Segment for which Tag Type should be determined.
 *
 * Returns Tag Type of segment.
 ***************************************************************************/
GXml::TagType GXml::get_tagtype(const std::string& segment) const
{
    // Initialise with invalid Tag Type
    TagType type = TT_INVALID;
    
    // Get length of segment
    int n = segment.length();

    // Check for comment
    if (n >= 7 && (segment.compare(0,4,"<!--") == 0) &&
             (segment.compare(n-3,3,"-->") == 0))
        type = TT_COMMENT;

    // Check for declaration
    else if (n >= 7 && (segment.compare(0,5,"<?xml") == 0) &&
             (segment.compare(n-2,2,"?>") == 0))
        type = TT_DECLARATION;

    // Check for processing instruction
    else if (n >= 4 && (segment.compare(0,2,"<?") == 0) &&
             (segment.compare(n-2,2,"?>") == 0))
        type = TT_PROCESSING;

    // Check for empty element tag
    else if (n >= 3 && (segment.compare(0,1,"<") == 0) &&
             (segment.compare(n-2,2,"/>") == 0))
        type = TT_ELEMENT_EMPTY;

    // Check for element end tag
    else if (n >= 3 && (segment.compare(0,2,"</") == 0) &&
             (segment.compare(n-1,1,">") == 0))
        type = TT_ELEMENT_END;

    // Check for element start tag
    else if (n >= 2 && (segment.compare(0,1,"<") == 0) &&
             (segment.compare(n-1,1,">") == 0))
        type = TT_ELEMENT_START;

    // Return type
    return type;
}


/***********************************************************************//**
 * @brief Check if character is whitespace
 *
 * @param[in] c Character to check.
 *
 * Checks if character is a XML 1.1 whitespace character
 ***************************************************************************/
bool GXml::is_whitespace(const int& c) const
{
    // Set result
    bool result = (c == '\x20' ||    // Whitespace
                   c == '\x09' ||    // Tab
                   c == '\x0d' ||    // Carriage return
                   c == '\x0a' ||    // Line feed
                   c == '\x85' ||    // Next line
                   c == L'\x2028');  // Line separator

    // Return result
    return result;
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

