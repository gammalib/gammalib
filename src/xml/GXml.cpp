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
#include "GXml.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                                     "GXml::load(std::string&)"

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
    m_name.clear();

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
    m_name  = xml.m_name;

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
 ***************************************************************************/
void GXml::parse(FILE* fptr)
{
    // Define parser modes
    enum Mode { OPEN,          // Wait for next '<'
                CLOSE };       // Wait for next '>'

    // Initialise parser
    int  c;
    int  level = 0;
    Mode mode  = OPEN;

    // Main parsing loop
    while ((c = fgetc(fptr)) != EOF) {

        // Perform mode dependent action
        switch (mode) {
        case OPEN:
            if (c == '<') {
                mode = CLOSE;
                std::cout << "CLOSE" << level;
            }
            else if (!is_whitespace(c)) {
                std::cout << " ***ERROR*** ";
            }
            break;
        case CLOSE:
            if (c == '>') {
                mode = OPEN;
                std::cout << "OPEN" << level;
            }
            break;
        default:
            break;
        }

        // Get next character
        std::cout << (char)c;

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if character is whitespace
 *
 * @param[in] c Character to check.
 ***************************************************************************/
bool GXml::is_whitespace(const int& c)
{
    // Set result
    bool result = (c == ' ');

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

    // Return output stream
    return os;
}

