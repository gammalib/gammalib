/***************************************************************************
 *               GCTAEventAtom.cpp  -  CTA event atom class                *
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
 * @file GCTAEventAtom.cpp
 * @brief GCTAEventAtom class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GCTAException.hpp"
#include "GCTAEventAtom.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MODEL                    "GCTAEventAtom::model(GModels&, GVector*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GCTAEventAtom::GCTAEventAtom() : GEventAtom()
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
GCTAEventAtom::GCTAEventAtom(const GCTAEventAtom& atom) : GEventAtom(atom)
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
GCTAEventAtom::~GCTAEventAtom()
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
GCTAEventAtom& GCTAEventAtom::operator= (const GCTAEventAtom& atom)
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
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return model value and gradient
 *
 * @param[in] models Model descriptor.
 * @param[out] gradient Pointer to gradient vector.
 *
 * @exception GCTAException::response_not_set
 *            Response function has not been set.
 * @exception GException::gradient_par_mismatch
 *            Gradient dimension mismatches number of parameters.
 *
 * Implements generic model and gradient evaluation for the CTA instrument.
 *
 * @todo Implement CTA specific model filtering.
 * @todo Implement correct gradient computation.
 ***************************************************************************/
double GCTAEventAtom::model(GModels& models, GVector* gradient) const
{
    // Make sure that response pointer exists
    if (rsp() == NULL)
        throw GCTAException::response_not_set(G_MODEL);

    // Verify that number of model parameter is identical to the dimension
    // of the gradient vector
    #if defined(G_RANGE_CHECK)
    if (models.npars() != gradient->size())
        throw GException::gradient_par_mismatch(G_MODEL, gradient->size(), 
                                                models.npars());
    #endif

    // Initialise model
    double model = 0.0;

    // Loop over models
    for (int i = 0; i < models.size(); ++i) {

        // Check if model is a CTA model and if it should be used for this
        // event
        // TO BE IMPLEMENTED

        // Add model
        model += models(i)->eval_gradients(*dir(), *energy(), *time(), *rsp(), *pnt());

    }

    // Optionally set gradient vector
    if (gradient != NULL) {
        for (int i = 0; i < gradient->size(); ++i) {
            double grad    = models.par(i)->gradient();
            (*gradient)(i) = (std::isinf(grad)) ? 0.0 : grad;
        }
    }

    // Return
    return model;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 *
 * @todo Need to implement GCTAInstDir:clear() and GCTAPointing:clear()
 ***************************************************************************/
void GCTAEventAtom::init_members(void)
{
    // Initialise CTA data format attributes
    //m_dir.clear();
    //m_pnt.clear();
    m_rsp         = NULL;
    m_event_id    = 0;
    m_flags       = 0;
    m_multip      = 0;
    m_telmask     = 0; 
    m_dir_err     = 0.0;
    m_detx        = 0.0;
    m_dety        = 0.0;
    m_alt_pnt     = 0.0;
    m_az_pnt      = 0.0;
    m_alt         = 0.0;
    m_az          = 0.0;
    m_corex       = 0.0;
    m_corey       = 0.0;
    m_core_err    = 0.0;
    m_xmax        = 0.0;
    m_xmax_err    = 0.0;
    m_energy_err  = 0.0;
    m_hil_msw     = 0.0;
    m_hil_msw_err = 0.0;
    m_hil_msl     = 0.0;
    m_hil_msl_err = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] atom GCTAEventAtom members which should be copied.
 ***************************************************************************/
void GCTAEventAtom::copy_members(const GCTAEventAtom& atom)
{
    // Copy CTA data format attributes
    m_dir         = atom.m_dir;
    m_pnt         = atom.m_pnt;
    m_rsp         = atom.m_rsp;
    m_event_id    = atom.m_event_id;
    m_flags       = atom.m_flags;
    m_multip      = atom.m_multip;
    m_telmask     = atom.m_telmask; 
    m_dir_err     = atom.m_dir_err;
    m_detx        = atom.m_detx;
    m_dety        = atom.m_dety;
    m_alt_pnt     = atom.m_alt_pnt;
    m_az_pnt      = atom.m_az_pnt;
    m_alt         = atom.m_alt;
    m_az          = atom.m_az;
    m_corex       = atom.m_corex;
    m_corey       = atom.m_corey;
    m_core_err    = atom.m_core_err;
    m_xmax        = atom.m_xmax;
    m_xmax_err    = atom.m_xmax_err;
    m_energy_err  = atom.m_energy_err;
    m_hil_msw     = atom.m_hil_msw;
    m_hil_msw_err = atom.m_hil_msw_err;
    m_hil_msl     = atom.m_hil_msl;
    m_hil_msl_err = atom.m_hil_msl_err;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventAtom::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GCTAEventAtom* GCTAEventAtom::clone(void) const
{
    return new GCTAEventAtom(*this);
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put atom into output stream
 *
 * @param[in] os Output stream into which the atom will be dumped
 * @param[in] atom Atom to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTAEventAtom& atom)
{
    // Put bin in output stream
    os.precision(3);
    os << std::fixed;
    os << "Time=" << atom.m_time;
    os << " Energy=" << atom.m_energy;
    os << atom.m_dir; 
        
    // Return output stream
    return os;
}
