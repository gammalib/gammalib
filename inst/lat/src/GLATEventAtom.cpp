/***************************************************************************
 *               GLATEventAtom.cpp  -  LAT event atom class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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
 * @file GLATEventAtom.cpp
 * @brief GLATEventAtom class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GLATEventAtom.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COPY_MEMBERS    "GLATEventAtom::copy_members(const GLATEventAtom&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                   GLATEventAtom constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GLATEventAtom::GLATEventAtom() : GEventAtom()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] atom Event atom from which the instance should be built.
 ***************************************************************************/
GLATEventAtom::GLATEventAtom(const GLATEventAtom& atom) : GEventAtom(atom)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(atom);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATEventAtom::~GLATEventAtom()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GLATEventAtom operators                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] atom Event atom which should be assigned.
 ***************************************************************************/
GLATEventAtom& GLATEventAtom::operator= (const GLATEventAtom& atom)
{
    // Execute only if object is not identical
    if (this != &atom) {

        // Copy base class members
        this->GEventAtom::operator=(atom);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(atom);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                       GLATEventAtom public methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return model value
 *
 * @param[in] models Model descriptor.
 ***************************************************************************/
double GLATEventAtom::model(GModels& models)
{
    // DUMMY: Set model value
    double model = 1.0;
    
    // Return
    return model;
}


/***********************************************************************//**
 * @brief Return model value and gradient
 *
 * @param[in] models Model descriptor.
 * @param[out] gradient Pointer to gradient vector.
 ***************************************************************************/
double GLATEventAtom::model(GModels& models, GVector* gradient)
{
    // DUMMY: Set model value
    double model = 1.0;
    
    // Set default gradient
    for (int i = 0; i < gradient->size(); ++i)
        (*gradient)(i) = model;

    // Return
    return model;
}


/*==========================================================================
 =                                                                         =
 =                       GLATEventAtom private methods                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATEventAtom::init_members(void)
{
    // Initialise LAT data format attributes
    m_theta               = 0.0;
    m_phi                 = 0.0;
    m_zenith_angle        = 0.0;
    m_earth_azimuth_angle = 0.0;
    m_event_id            = 0;
    m_run_id              = 0;
    m_recon_version       = 0;
    m_calib_version[0]    = 0;
    m_calib_version[1]    = 0;
    m_calib_version[2]    = 0;
    m_event_class         = 0;
    m_conversion_type     = 0;
    m_livetime            = 0.0;
    m_difrsp              = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] atom GLATEventAtom members which should be copied.
 ***************************************************************************/
void GLATEventAtom::copy_members(const GLATEventAtom& atom)
{
    // Copy LAT data format attributes
    m_theta               = atom.m_theta;
    m_phi                 = atom.m_phi;
    m_zenith_angle        = atom.m_zenith_angle;
    m_earth_azimuth_angle = atom.m_earth_azimuth_angle;
    m_event_id            = atom.m_event_id;
    m_run_id              = atom.m_run_id;
    m_recon_version       = atom.m_recon_version;
    m_calib_version[0]    = atom.m_calib_version[0];
    m_calib_version[1]    = atom.m_calib_version[1];
    m_calib_version[2]    = atom.m_calib_version[2];
    m_event_class         = atom.m_event_class;
    m_conversion_type     = atom.m_conversion_type;
    m_livetime            = atom.m_livetime;

    // Copy other attributes
    m_num_difrsp = atom.m_num_difrsp;

    // If there are diffuse response components then copy them
    if (m_num_difrsp > 0 && atom.m_difrsp != NULL) {

        // Allocate memory for diffuse response components
        m_difrsp = new double[m_num_difrsp];
        if (m_difrsp == NULL)
            throw GException::mem_alloc(G_COPY_MEMBERS, m_num_difrsp);

        // Copy diffuse response components
        for (int i = 0; i < m_num_difrsp; ++i)
            m_difrsp[i] = atom.m_difrsp[i];

    } // endif: there were diffuse response components to copy

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEventAtom::free_members(void)
{
    // Free memory
    if (m_difrsp != NULL) delete [] m_difrsp;

    // Signal free pointers
    m_difrsp = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GLATEventAtom* GLATEventAtom::clone(void) const
{
    return new GLATEventAtom(*this);
}


/*==========================================================================
 =                                                                         =
 =                          GLATEventAtom friends                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                 Other functions used by GLATEventAtom                   =
 =                                                                         =
 ==========================================================================*/
