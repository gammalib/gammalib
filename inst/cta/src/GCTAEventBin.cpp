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
#include <cmath>
#include "GException.hpp"
#include "GCTAException.hpp"
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
 * @brief Return size of event bin
 ***************************************************************************/
double GCTAEventBin::size(void) const
{
    // Compute bin size
    double size = *omega() * ewidth()->MeV() * *ontime();

    // Return bin size
    return size;
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
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAEventBin::init_members(void)
{
    // Initialise CTA specific attributes
    m_dir    = NULL;
    m_omega  = NULL;
    m_ewidth = NULL;
    m_ontime = NULL;

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
    m_dir    = bin.m_dir;
    m_omega  = bin.m_omega;
    m_ewidth = bin.m_ewidth;
    m_ontime = bin.m_ontime;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * This class does not allocate any memory but simply holds pointers. Hence
 * nothing has to be deallocated.
 ***************************************************************************/
void GCTAEventBin::free_members(void)
{
    // Return
    return;
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
