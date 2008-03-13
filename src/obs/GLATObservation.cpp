/***************************************************************************
 *               GLATObservation.cpp  -  LAT Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GLATObservation.cpp
 * @brief GLATObservationclass implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GLATObservation.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                  GLATObservation constructors/destructors               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GLATObservation::GLATObservation() : GObservation()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] ft1name FT1 FITS filename.
 * @param[in] ft2name FT2 FITS filename.
 ***************************************************************************/
GLATObservation::GLATObservation(const std::string& ft1name, 
                                 const std::string& ft2name) : GObservation()
{
    // Initialise class members for clean destruction
    init_members();
    
    // Allocate FT1 & FT2 FITS files
    m_ft1 = new GFits;
    m_ft2 = new GFits;
    
    // Open FT1 & FT2 FITS files
    m_ft1->open(ft1name);
    m_ft2->open(ft2name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs Observation from which the instance should be built.
 ***************************************************************************/
GLATObservation::GLATObservation(const GLATObservation& obs) : GObservation(obs)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATObservation::~GLATObservation()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        GLATObservation operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] obs Observation which should be assigned.
 ***************************************************************************/
GLATObservation& GLATObservation::operator= (const GLATObservation& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Copy base class members
        this->GObservation::operator=(obs);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(obs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                     GLATObservation public methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return pointer to FT1 file.
 ***************************************************************************/
GFits* GLATObservation::ft1(void) const
{
    // Return FT1 pointer
    return m_ft1;
}


/***********************************************************************//**
 * @brief Return pointer to FT2 file.
 ***************************************************************************/
GFits* GLATObservation::ft2(void) const
{
    // Return FT2 pointer
    return m_ft2;
}


/*==========================================================================
 =                                                                         =
 =                     GLATObservation private methods                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATObservation::init_members(void)
{
    // Initialise members
    m_ft1    = NULL;
    m_ft2    = NULL;
    m_ltcube = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation to be copied
 ***************************************************************************/
void GLATObservation::copy_members(const GLATObservation& obs)
{
    // Copy members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATObservation::free_members(void)
{
    // Free members
    if (m_ft1    != NULL) delete m_ft1;
    if (m_ft2    != NULL) delete m_ft2;
    if (m_ltcube != NULL) delete m_ltcube;

    // Mark memory as free
    m_ft1    = NULL;
    m_ft2    = NULL;
    m_ltcube = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GLATObservation* GLATObservation::clone(void) const
{
    return new GLATObservation(*this);
}


/*==========================================================================
 =                                                                         =
 =                         GLATObservation friends                         =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                 Other functions used by GLATObservation                 =
 =                                                                         =
 ==========================================================================*/
