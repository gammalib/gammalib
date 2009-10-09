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
#include "GFits.hpp"

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


/***********************************************************************//**
 * @brief Clone GTI
 *
 * Cloning provides a copy of the GTIs. Cloning is used to allocate
 * derived classes into a base class pointer.
 ***************************************************************************/
GGti* GGti::clone(void) const
{
    // Clone this image
    return new GGti(*this);
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

/***********************************************************************//**
 * @brief Load GTI intervals from file.
 ***************************************************************************/
void GGti::load(const std::string& filename)
{
	// Free members
	free_members();

	// Initialise attributes
	init_members();
	
	// Allocate FITS file
	GFits file;
	
	// Open FITS file
	file.open(filename);
	
	// Get GTI HDU
	GFitsHDU* hdu = file.hdu("GTI");
	
	// Extract GTI information from FITS file
	m_num = hdu->card("NAXIS2")->integer();
	if (m_num > 0) {
	
		// Get GTI intervals
		m_start = new double[m_num];
		m_stop  = new double[m_num];
		for (int i = 0; i < m_num; ++i) {
			m_start[i] = hdu->column("START")->real(i);
			m_stop[i]  = hdu->column("STOP")->real(i);
			m_ontime  += (m_stop[i] - m_start[i]);
		}
		
		// Set attributes
		m_tstart = m_start[0];
		m_tstop  = m_stop[m_num-1];
		m_elapse = m_tstop - m_tstart;

	}
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return start time
 ***************************************************************************/
double GGti::tstart(void)
{
    // Return
    return m_tstart;
}


/***********************************************************************//**
 * @brief Return stop time
 ***************************************************************************/
double GGti::tstop(void)
{
    // Return
    return m_tstop;
}


/***********************************************************************//**
 * @brief Return ontime
 ***************************************************************************/
double GGti::ontime(void)
{
    // Return
    return m_ontime;
}


/***********************************************************************//**
 * @brief Return elapsed time
 ***************************************************************************/
double GGti::elapse(void)
{
    // Return
    return m_elapse;
}


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

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream into which the GTI will be dumped
 * @param[in] gti GTI to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GGti& gti)
{
    // Put header in stream
    os << "=== GGti ===" << std::endl;
    os << " Number of intervals .......: " << gti.m_num << std::endl;
    os << " Start time ................: " << gti.m_tstart << std::endl;
    os << " Stop time .................: " << gti.m_tstop << std::endl;
    os << " Ontime ....................: " << gti.m_ontime << " sec" << std::endl;
    os << " Elapsed time ..............: " << gti.m_elapse << " sec" << std::endl;
	
    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                       Other functions used by GGti                      =
 =                                                                         =
 ==========================================================================*/
