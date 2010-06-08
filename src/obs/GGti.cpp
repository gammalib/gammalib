/***************************************************************************
 *                 GGti.cpp  -  Good time interval class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GGti.cpp
 * @brief Good time interval class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <iomanip.h>
#include "GException.hpp"
#include "GGti.hpp"
#include "GFits.hpp"

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
GGti::GGti(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] gti Object from which the instance should be built.
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
GGti::~GGti(void)
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear region of interest
 ***************************************************************************/
void GGti::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Good Time Interval
 *
 * @param[in] tstart Start time of interval to be appended.
 * @param[in] tstop Stop time of interval to be appended.
 *
 * This method appends a new GTI to the object. If the specified time
 * interval overlaps with an existing time interval, the existing time
 * interval will be extended. If overlap with more than a single interval
 * exists, the intervals will be merged. Otherwise a new interval will
 * be inserted, respecting the time ordering of the intervals. If the
 * time interval is not valid (tstop <= tstart), nothing is done.
 *
 * @todo Throw error if specified time interval is not valid.
 * @todo Method not yet implemented.
 ***************************************************************************/
void GGti::append(const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {
    
        // If GTI is empty then append interval ...
        if (m_num < 1)
            insert(0, tstart, tstop);

        // ... otherwise check for overlaps and perform proper action
        else {
        }
    
    } // endif: Time interval was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load GTI intervals from file.
 *
 * @param[in] filename Name of file from which GTIs are to be loaded.
 *
 * @todo Method assumes that times are in MET. Probably this method should
 * not exist but be instrument specific.
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
	
		// Set GTIs
		m_start = new GTime[m_num];
		m_stop  = new GTime[m_num];
		for (int i = 0; i < m_num; ++i) {
			m_start[i].met(hdu->column("START")->real(i));
			m_stop[i].met(hdu->column("STOP")->real(i));
		}
		
		// Set attributes
        set_attributes();

	}
	
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
void GGti::init_members(void)
{
    // Initialise members
    m_num     = 0;
	m_tstart.clear();
	m_tstop.clear();
	m_ontime  = 0.0;
	m_telapse = 0.0;
	m_start   = NULL;
	m_stop    = NULL;

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
    m_num     = gti.m_num;
    m_tstart  = gti.m_tstart;
    m_tstop   = gti.m_tstop;
    m_ontime  = gti.m_ontime;
    m_telapse = gti.m_telapse;

    // Copy start/stop times
    if (m_num > 0) {
        m_start = new GTime[m_num];
        m_stop  = new GTime[m_num];
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


/***********************************************************************//**
 * @brief Set class attributes
 ***************************************************************************/
void GGti::set_attributes(void)
{
    // Continue only if there are GTIs
    if (m_num > 0) {
    
        // Set attributes
        m_tstart  = m_start[0];
        m_tstop   = m_stop[m_num-1];
        m_telapse = m_tstop.met() - m_tstart.met();
        m_ontime  = 0.0;
		for (int i = 0; i < m_num; ++i)
			m_ontime += (m_stop[i].met() - m_start[i].met());

    } // endif: there were GTIs
    
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


/***********************************************************************//**
 * @brief Insert Good Time Interval
 *
 * @param[in] inx Index at which GTI is to be inserted.
 * @param[in] tstart Start time of interval to be appended.
 * @param[in] tstop Stop time of interval to be appended.
 *
 * @todo Check than inx is valid.
 ***************************************************************************/
void GGti::insert(int inx, const GTime& tstart, const GTime& tstop)
{
    // Continue only if time interval is valid
    if (tstop > tstart) {
    
        // Allocate new intervals
        int    num   = m_num+1;
        GTime* start = new GTime[m_num];
        GTime* stop  = new GTime[m_num];
        
        // Copy intervals before GTI to be inserted
        for (int i = 0; i < inx; ++i) {
            start[i] = m_start[i];
            stop[i]  = m_stop[i];
        }
        
        // Insert GTI
        start[inx] = tstart;
        stop[inx]  = tstop;

        // Copy intervals after GTI to be inserted
        for (int i = inx+1; i < num; ++i) {
            start[i] = m_start[i-1];
            stop[i]  = m_stop[i-1];
        }
        
        // Free memory
        if (m_start != NULL) delete [] m_start;
        if (m_stop  != NULL) delete [] m_stop;
        
        // Set new memory
        m_start = start;
        m_stop  = stop;
        
        // Set attributes
        m_num = num;
        set_attributes();
    
    } // endif: Time interval was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
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
    // Put GTIs in stream
    os << "=== GGti ===" << std::endl;
    os << " Number of intervals .......: " << gti.m_num << std::endl;
    os << " Ontime ....................: " << gti.ontime() << " sec" << std::endl;
    os << " Elapsed time ..............: " << gti.telapse() << " sec" << std::endl;
    os << " Time range ................: " << gti.tstart() << " - " << gti.tstop();
	
    // Return output stream
    return os;
}
