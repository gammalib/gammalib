/***************************************************************************
 *                GCTAEventCube.cpp  -  CTA event cube class               *
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
 * @file GCTAEventCube.cpp
 * @brief GCTAEventCube class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GCTAEventCube.hpp"
#include "GFitsImageDbl.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_POINTER                               "GCTAEventCube::pointer(int)"

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
GCTAEventCube::GCTAEventCube(void) : GEventCube()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Event cube from which the instance should be built.
 ***************************************************************************/
GCTAEventCube::GCTAEventCube(const GCTAEventCube& cube) : GEventCube(cube)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAEventCube::~GCTAEventCube(void)
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
 * @param[in] cube Event cube to be assigned.
 ***************************************************************************/
GCTAEventCube& GCTAEventCube::operator= (const GCTAEventCube& cube)
{
    // Execute only if object is not identical
    if (this != &cube) {

        // Copy base class members
        this->GEventCube::operator=(cube);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(cube);

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
 * @brief Load CTA counts map from FITS file.
 *
 * @param[in] filename Counts map FITS filename to be loaded.
 *
 * It is assumed that the counts map resides in the primary extension of the
 * FITS file, the energy boundaries reside in the EBOUNDS extension and the
 * Good Time Intervals reside in the GTI extension.
 ***************************************************************************/
void GCTAEventCube::load(const std::string& filename)
{
    // Free and initialise base class members
    this->GEvents::free_members();
    this->GEvents::init_members();

    // Free and initialise base class members
    this->GEventCube::free_members();
    this->GEventCube::init_members();

    // Free and initialise class members
    free_members();
    init_members();

    // Allocate FITS file
    GFits file;

    // Open counts map FITS file
    file.open(filename);

    // Get HDUs
    GFitsHDU* hdu_cntmap  = file.hdu("Primary");
    GFitsHDU* hdu_ebounds = file.hdu("EBOUNDS");
    GFitsHDU* hdu_gti     = file.hdu("GTI");

    // Load counts map
    load_cntmap(hdu_cntmap);

    // Load energy boundaries
    load_ebds(hdu_ebounds);

    // Load energy boundaries
    load_gti(hdu_gti);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get pointer to element
 *
 * @param[in] index Event index for which pointer will be returned.
 *
 * This method provides the event attributes to the event bin. It 
 *
 * @todo Implement conversion routine from event cube index to direction.
 ***************************************************************************/
GCTAEventBin* GCTAEventCube::pointer(int index)
{
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_elements)
        throw GException::out_of_range(G_POINTER, index, 0, m_elements-1);
    #endif

    // Get pixel and energy bin
    int ipix = index / ebins();
    int ieng = index % pixels();

    // Set GEventBin pointers
    m_bin.m_counts = &(m_counts[index]);
    m_bin.m_time   = &m_time;
    m_bin.m_energy = &(m_energies[ieng]);

    // Set GCTAEventBin pointers
    m_bin.m_dir    = &(m_dirs[ipix]);
    m_bin.m_pnt    = &m_pnt;
    m_bin.m_rsp    = &m_rsp;
    
    // Return pointer
    return (GCTAEventBin*)&(m_bin);
}


/***********************************************************************//**
 * @brief Return number of events in cube
 ***************************************************************************/
int GCTAEventCube::number(void) const
{
    // Initialise result
    double number = 0.0;

    // Sum event cube
    if (m_elements > 0 && m_counts != NULL) {
        for (int i = 0; i < m_elements; ++i)
            number += m_counts[i];
    }

    // Return
    return int(number+0.5);
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 *
 * @todo Implement GSkymap.clear(), GCTAEventBin.clear(), and 
 *       GCTAResponse.clear() methods
 ***************************************************************************/
void GCTAEventCube::init_members(void)
{
    // Initialise members
    //m_bin.clear();
    //m_map.clear();
    m_bin      = GCTAEventBin();
    m_map      = GSkymap();
    m_counts   = NULL;
    m_dirs     = NULL;
    m_energies = NULL;
    m_time.clear();
    m_pnt.clear();
    //m_rsp.clear();
    m_rsp      = GCTAResponse();
    m_ebds.clear();
    m_gti.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube GCTAEventCube members which should be copied.
 ***************************************************************************/
void GCTAEventCube::copy_members(const GCTAEventCube& cube)
{
    // Copy CTA specific attributes
    m_bin  = cube.m_bin;
    m_map  = cube.m_map;
    m_time = cube.m_time;
    m_pnt  = cube.m_pnt;
    m_rsp  = cube.m_rsp;
    m_ebds = cube.m_ebds;
    m_gti  = cube.m_gti;

    // Set counter to copied skymap pixels
    m_counts = m_map.pixels();

    // Copy bin energies
    if (cube.ebins() > 0 && cube.m_energies != NULL) {
        m_energies = new GEnergy[cube.ebins()];
        for (int i = 0; i < cube.ebins(); ++i)
            m_energies[i] = cube.m_energies[i];
    }

    // Copy sky directions
    if (cube.pixels() > 0 && cube.m_dirs != NULL) {
        m_dirs = new GCTAInstDir[cube.pixels()];
        for (int i = 0; i < cube.pixels(); ++i)
            m_dirs[i] = cube.m_dirs[i];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventCube::free_members(void)
{
    // Free memory
    if (m_energies != NULL) delete [] m_energies;
    if (m_dirs     != NULL) delete [] m_dirs;

    // Signal free pointers
    m_energies = NULL;
    m_dirs     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GCTAEventCube* GCTAEventCube::clone(void) const
{
    return new GCTAEventCube(*this);
}


/***********************************************************************//**
 * @brief Load CTA counts map from HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which events are loaded.
 *
 * This method load a CTA counts map from a FITS HDU. The counts map is
 * stored in a GSkymap object, and a pointer is set up to access the pixels
 * individually. Recall that skymap pixels are stored in the order
 * (ebin,ix,iy), i.e. the energy axis is the most rapidely varying axis,
 * while the counts map is stored in the order (ix,iy,ebin), i.e. the x
 * axis is the most rapidely varying axis.
 ***************************************************************************/
void GCTAEventCube::load_cntmap(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Load counts map as sky map
        m_map.read(hdu);
        
        // Set GEventCube attributes
        m_dim      = (m_map.nmaps() > 1) ? 3 : 2;
        m_naxis    = new int[m_dim];
        m_naxis[0] = m_map.nx();
        m_naxis[1] = m_map.ny();
        m_elements = m_map.npix();
        if (m_dim == 3) {
            m_naxis[2]  = m_map.nmaps();
            m_elements *= m_naxis[2];
        }
        
        // Set GCTAEventCube attributes
        m_counts = m_map.pixels();
        
        // Set pixel directions
        if (m_dirs != NULL) delete [] m_dirs;
        m_dirs = new GCTAInstDir[m_map.npix()];
        for (int i = 0, iy = 0; iy < m_naxis[1]; ++iy) {
            for (int ix = 0; ix < m_naxis[0]; ++ix, ++i)
                m_dirs[i] = GCTAInstDir(m_map.xy2dir(GSkyPixel(double(ix), double(iy))));
        }

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy boundaries from HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which energy boundaries are loaded.
 ***************************************************************************/
void GCTAEventCube::load_ebds(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Load energy boundaries
        m_ebds.load(hdu);

        // Set log mean energy of all bins
        if (m_ebds.size() > 0) {

            // Setup bin energies
            if (m_energies != NULL) delete [] m_energies;
            m_energies = new GEnergy[m_ebds.size()];

            // Set log mean energies
            for (int i = 0; i < m_ebds.size(); ++i)
                m_energies[i] = m_ebds.elogmean(i);
        
        } // endif: there were energy bins

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load GTIs from HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which GTIs are loaded.
 *
 * @todo Implement method.
 ***************************************************************************/
void GCTAEventCube::load_gti(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

    } // endif: HDU was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GCTAEventCube friends                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put CTA event cube in output stream
 *
 * @param[in] os Output stream into which the event cube will be dumped
 * @param[in] cube Event cube to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTAEventCube& cube)
{
    // Put CTA event list in output stream
    os << "=== GCTAEventCube ===" << std::endl;
    os << " Number of elements ........: " << cube.size() << std::endl;
    os << " Number of pixels ..........: " << cube.pixels() << std::endl;
    os << " Number of energy bins .....: " << cube.ebins() << std::endl;
    os << " Number of events ..........: " << cube.number() << std::endl;
    os << cube.m_ebds;

    // Return output stream
    return os;
}
