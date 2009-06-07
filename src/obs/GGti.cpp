/***************************************************************************
 *                 GGti.cpp  -  Good time interval class                   *
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
 * @file GGti.cpp
 * @brief Good time interval class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GGti.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COPY_MEMBERS      "GGti::copy_members(const GGti&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                        GGti constructors/destructors                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GGti::GGti()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] gti Good time interval from which the instance should be built.
 ***************************************************************************/
GGti::GGti(const GGti& gti)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(gti);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GGti::~GGti()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              GGti operators                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] gti Good time interval which should be assigned.
 ***************************************************************************/
GGti& GGti::operator= (const GGti& gti)
{
    // Execute only if object is not identical
    if (this != &gti) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(gti);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                           GGti public methods                           =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                           GGti private methods                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GGti::init_members(void)
{
    // Initialise members
    m_num    = 0;
	m_tstart = 0.0;
	m_tstop  = 0.0;
	m_ontime = 0.0;
	m_elapse = 0.0;
	m_start  = NULL;
	m_stop   = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] gti GGti members which should be copied.
 ***************************************************************************/
void GGti::copy_members(const GGti& gti)
{
    // Copy attributes
    m_num    = gti.m_num;
    m_tstart = gti.m_tstart;
    m_tstop  = gti.m_tstop;
    m_ontime = gti.m_ontime;
    m_elapse = gti.m_elapse;

    // Copy start/stop times
    if (m_num > 0) {
        m_start = new double[m_num];
        m_stop  = new double[m_num];
        if (m_start == NULL || m_stop == NULL)
            throw GException::mem_alloc(G_COPY_MEMBERS, m_num);
        for (int i = 0; i < m_num; ++i) {
            m_start[i] = gti.m_start[i];
            m_stop[i]  = gti.m_stop[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GGti::free_members(void)
{
    // Free memory
    if (m_start != NULL) delete [] m_start;
    if (m_stop  != NULL) delete [] m_stop;

    // Signal free pointers
    m_start = NULL;
    m_stop  = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               GGti friends                              =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                       Other functions used by GGti                      =
 =                                                                         =
 ==========================================================================*/
