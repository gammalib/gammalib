/***************************************************************************
 *       GXmlAttribute.cpp - XML attribute node class implementation       *
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
 * @file GXmlAttribute.cpp
 * @brief XML attribute node class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GXmlAttribute.hpp"

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
GXmlAttribute::GXmlAttribute(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] attr Object from which the instance should be built.
 ***************************************************************************/
GXmlAttribute::GXmlAttribute(const GXmlAttribute& attr)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(attr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Name-Value pair constructor
 *
 * @param[in] name Name for instance building.
 * @param[in] value Value for instance building.
 ***************************************************************************/
GXmlAttribute::GXmlAttribute(const std::string& name, const std::string& value)
{
    // Initialise private members for clean destruction
    init_members();

    // Set attribute
    m_name  = name;
    m_value = value;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlAttribute::~GXmlAttribute(void)
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
 * @param[in] attr Object which should be assigned.
 ***************************************************************************/
GXmlAttribute& GXmlAttribute::operator= (const GXmlAttribute& attr)
{
    // Execute only if object is not identical
    if (this != &attr) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(attr);

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
void GXmlAttribute::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write attribute into file
 *
 * @param[in] fptr File pointer.
 ***************************************************************************/
void GXmlAttribute::write(FILE* fptr) const
{
    // Write attribute into file
    fprintf(fptr, " %s=%s", m_name.c_str(), m_value.c_str());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print attribute in output stream
 *
 * @param[in] os Output stream into which the node will be printed.
 ***************************************************************************/
void GXmlAttribute::print(std::ostream& os) const
{
    // Put attribute in output stream
    os << " " << m_name << "=" << m_value;

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
void GXmlAttribute::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_value.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node Object from which members which should be copied.
 ***************************************************************************/
void GXmlAttribute::copy_members(const GXmlAttribute& attr)
{
    // Copy attributes
    m_name  = attr.m_name;
    m_value = attr.m_value;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXmlAttribute::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GXmlAttribute* GXmlAttribute::clone(void) const
{
    return new GXmlAttribute(*this);
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

