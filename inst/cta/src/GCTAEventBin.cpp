/***************************************************************************
 *                GCTAEventBin.cpp  -  CTA event bin class                 *
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
 * @file GCTAEventBin.cpp
 * @brief GCTAEventBin class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GCTAEventBin.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MODEL                     "GCTAEventBin::model(GModels&, GVector*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAEventBin::GCTAEventBin(void) : GEventBin()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bin Event bin from which the instance should be built.
 ***************************************************************************/
GCTAEventBin::GCTAEventBin(const GCTAEventBin& bin) : GEventBin(bin)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(bin);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAEventBin::~GCTAEventBin(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] bin Event bin which should be assigned.
 ***************************************************************************/
GCTAEventBin& GCTAEventBin::operator= (const GCTAEventBin& bin)
{
    // Execute only if object is not identical
    if (this != &bin) {

        // Copy base class members
        this->GEventBin::operator=(bin);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(bin);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return model value and (optionally) gradient
 *
 * @param[in] models Model descriptor.
 * @param[out] gradient Pointer to gradient vector (NULL=not computed).
 *
 * @exception GException::gradient_par_mismatch
 *            Dimension of gradient vector mismatches number of parameters.
 *
 * Implements generic model and gradient evaluation for the CTA instrument.
 *
 * @todo Requires multiplication of model with solid angle, duration, and
 *       energy binsize. Current method not correct.
 * @todo Requires implementation of all model types (not only factorized
 *       point sources which are currently the only type that is
 *       supported)
 ***************************************************************************/
double GCTAEventBin::model(GModels& models, GVector* gradient) const
{
    // Verify that gradients vector has the same dimension than the
    // model has parameters
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

    // Set gradient vector
    if (gradient != NULL) {
        for (int i = 0; i < gradient->size(); ++i)
            (*gradient)(i) = models.par(i)->gradient();
    }

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
void GCTAEventBin::init_members(void)
{
    // Initialise CTA specific attributes
    m_dir = NULL;
    m_pnt = NULL;
    m_rsp = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bin GCTAEventBin members which should be copied.
 ***************************************************************************/
void GCTAEventBin::copy_members(const GCTAEventBin& bin)
{
    // Copy CTA specific attributes
    m_dir = bin.m_dir;
    m_pnt = bin.m_pnt;
    m_rsp = bin.m_rsp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventBin::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GCTAEventBin* GCTAEventBin::clone(void) const
{
    return new GCTAEventBin(*this);
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put bin into output stream
 *
 * @param[in] os Output stream into which the bin will be dumped
 * @param[in] bin Bin to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTAEventBin& bin)
{
    // Put bin in output stream
    os << bin.m_counts << " ";
        
    // Return output stream
    return os;
}
