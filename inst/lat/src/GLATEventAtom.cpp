/***************************************************************************
 *               GLATEventAtom.cpp  -  LAT event atom class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATEventAtom.cpp
 * @brief GLATEventAtom class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GLATEventAtom.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MODEL                    "GLATEventAtom::model(GModels&, GVector*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GLATEventAtom::GLATEventAtom(void) : GEventAtom()
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
GLATEventAtom::~GLATEventAtom(void)
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
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return model value and gradient
 *
 * @param[in] models Model descriptor.
 * @param[out] gradient Pointer to gradient vector (NULL=not computed).
 *
 * @exception GException::gradient_par_mismatch
 *            Dimension of gradient vector mismatches number of parameters.
 *
 * Implements generic model and gradient evaluation for the LAT instrument.
 *
 * @todo Requires implementation of all model types (not only factorized
 *       point sources which are currently the only type that is
 *       supported)
 ***************************************************************************/
double GLATEventAtom::model(GModels& models, GVector* gradient) const
{
    // Verify that gradients vector has the same dimension than the
    // model has parameters
    #if defined(G_RANGE_CHECK)
    if (models.npars() != gradient->size())
        throw GException::gradient_par_mismatch(G_MODEL, gradient->size(), 
                                                models.npars());
    #endif

    // Integral over source direction, energy and time
    //for (srcTime = ...
    //for (srcEng = ...
    //for (srcDir = ...
    GTime   srcTime = *time();    // Assume no time dispersion
    GEnergy srcEng  = *energy();  // Assume no energy dispersion
    GSkyDir srcDir;               // Needs to be implemented

    // Get source term
    double source = models.eval_gradients(srcDir, srcEng, srcTime);
    
    // Get IRF
    double irf = rsp()->irf(*dir(), *energy(), *time(), 
                            srcDir, srcEng, srcTime, *pnt());
    
    // Evaluate model
    double model = source * irf;

    // Set gradient vector
    if (gradient != NULL) {
        for (int i = 0; i < gradient->size(); ++i)
            (*gradient)(i) = irf * models.par(i)->gradient();
    }

    //}
    //}
    //}

    // Return
    return model;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
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
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put atom into output stream
 *
 * @param[in] os Output stream into which the bin will be dumped
 * @param[in] bin Bin to be dumped
 *
 * @todo Needs to be implemented
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATEventAtom& atom)
{
    // Put atom in output stream
    os << "..." << " ";
        
    // Return output stream
    return os;
}
