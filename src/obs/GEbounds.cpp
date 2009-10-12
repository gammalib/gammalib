/***************************************************************************
 *                 GEbounds.cpp  -  Energy boundary class                  *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
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
 * @file GEbounds.cpp
 * @brief Energy boundary class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GEbounds.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COPY_MEMBERS              "GEbounds::copy_members(const GEbounds&)"
#define G_LOAD       "GEbounds::load(const std::string&, const std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                      GEbounds constructors/destructors                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GEbounds::GEbounds()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ebds Energy boundaries from which the instance should be built.
 ***************************************************************************/
GEbounds::GEbounds(const GEbounds& ebds)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(ebds);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEbounds::~GEbounds()
{
    // Free members
    free_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone GEbounds
 *
 * Cloning provides a copy of the GEbounds. Cloning is used to allocate
 * derived classes into a base class pointer.
 ***************************************************************************/
GEbounds* GEbounds::clone(void) const
{
    // Clone this image
    return new GEbounds(*this);
}


/*==========================================================================
 =                                                                         =
 =                            GEbounds operators                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] ebds Energy boundaries to be assigned.
 ***************************************************************************/
GEbounds& GEbounds::operator= (const GEbounds& ebds)
{
    // Execute only if object is not identical
    if (this != &ebds) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(ebds);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GEbounds public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Load energy boundaries from file.
 *
 * @param[in] filename FITS filename from which GEbounds is to be loaded.
 * @param[in] extname FITS extension name of the energy boundaries.
 ***************************************************************************/
void GEbounds::load(const std::string& filename, const std::string& extname)
{
	// Allocate FITS file
	GFits file;
	
	// Open FITS file
	file.open(filename);
	
	// Get energy boundary HDU
	GFitsHDU* hdu = file.hdu(extname);

    // Load energy boundaries
    load(hdu);

	// Close FITS file
	file.close();
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy boundaries from file.
 *
 * @param[in] hdu Pointer to FITS HDU from which GEbounds are loaded.
 ***************************************************************************/
void GEbounds::load(GFitsHDU* hdu)
{
	// Free members
	free_members();

	// Initialise attributes
	init_members();
	
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Extract energy boundary information from FITS file
        m_num = hdu->card("NAXIS2")->integer();
        if (m_num > 0) {
	
            // Allocate memory
            m_channel = new int[m_num];
            m_emin    = new double[m_num];
            m_emax    = new double[m_num];
            if (m_channel == NULL || m_emin == NULL || m_emax == NULL)
                throw GException::mem_alloc(G_LOAD, m_num);

            // Copy information
            for (int i = 0; i < m_num; ++i) {
                m_channel[i] = hdu->column("CHANNEL")->integer(i);
                m_emin[i]    = hdu->column("E_MIN")->real(i);
                m_emax[i]    = hdu->column("E_MAX")->real(i);
            }
        
        } // endif: there were channels to read
    
        // Get OGIP keywords
        m_telescope   = hdu->card("TELESCOP")->string();
        m_instrument  = hdu->card("INSTRUME")->string();
        m_filter      = hdu->card("FILTER")->string();
        m_chantype    = hdu->card("CHANTYPE")->string();
        m_detchannels = hdu->card("DETCHANS")->string();
    
        // Fix energy scale
        m_escale = 1.0;
    
    }
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get minimum energy.
 ***************************************************************************/
double GEbounds::emin(void) const
{
    // Determine minimum energy
    double emin = (m_emin != NULL) ? m_emin[0] : 0.0;
    
    // Return
    return emin;
}


/***********************************************************************//**
 * @brief Get maximum energy.
 ***************************************************************************/
double GEbounds::emax(void) const
{
    // Determine minimum energy
    double emax = (m_emax != NULL) ? m_emax[m_num-1] : 0.0;
    
    // Return
    return emax;
}


/*==========================================================================
 =                                                                         =
 =                         GEbounds private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEbounds::init_members(void)
{
    // Initialise members
    m_num     = 0;
	m_channel = NULL;
	m_emin    = NULL;
	m_emax    = NULL;
    m_escale  = 1.0;
    m_telescope.clear();
    m_instrument.clear();
    m_filter.clear();
    m_chantype.clear();
    m_detchannels.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ebds GEbounds members which should be copied.
 ***************************************************************************/
void GEbounds::copy_members(const GEbounds& ebds)
{
    // Copy attributes
    m_num         = ebds.m_num;
    m_escale      = ebds.m_escale;
    m_telescope   = ebds.m_telescope;
    m_instrument  = ebds.m_instrument;
    m_filter      = ebds.m_filter;
    m_chantype    = ebds.m_chantype;
    m_detchannels = ebds.m_detchannels;

    // Copy arrays
    if (m_num > 0) {
        m_channel = new int[m_num];
        m_emin    = new double[m_num];
        m_emax    = new double[m_num];
        if (m_channel == NULL || m_emin == NULL || m_emax == NULL)
            throw GException::mem_alloc(G_COPY_MEMBERS, m_num);
        for (int i = 0; i < m_num; ++i) {
            m_channel[i] = ebds.m_channel[i];
            m_emin[i]    = ebds.m_emin[i];
            m_emax[i]    = ebds.m_emax[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEbounds::free_members(void)
{
    // Free memory
    if (m_channel != NULL) delete [] m_channel;
    if (m_emin    != NULL) delete [] m_emin;
    if (m_emax    != NULL) delete [] m_emax;

    // Signal free pointers
	m_channel = NULL;
	m_emin    = NULL;
	m_emax    = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GEbounds friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put energy boundaries in output stream
 *
 * @param[in] os Output stream into which the GEbounds will be dumped
 * @param[in] ebds GEbounds to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GEbounds& ebds)
{
    // Put GEbounds in stream
    os.precision(3);
    os << "=== GEbounds ===" << std::endl;
    os << " Number of boundaries ......: " << ebds.m_num << std::endl;
    os << " Energy range ..............: " << std::fixed << 
            ebds.emin() << " - " << ebds.emax() << " keV" << std::endl;
    os << " Channel type ..............: " << ebds.m_chantype << std::endl;
    os << " Telescope .................: " << ebds.m_telescope << std::endl;
    os << " Instrument ................: " << ebds.m_instrument;
	
    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GEbounds                    =
 =                                                                         =
 ==========================================================================*/
