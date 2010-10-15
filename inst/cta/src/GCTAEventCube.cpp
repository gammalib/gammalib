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
#include "GCTAException.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponse.hpp"
#include "GFitsImageDbl.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_POINTER                               "GCTAEventCube::pointer(int)"
#define G_SET_DIRECTIONS                    "GCTAEventCube::set_directions()"
#define G_SET_ENERGIES                        "GCTAEventCube::set_energies()"
#define G_SET_TIME                                "GCTAEventCube::set_time()"

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
 * @brief Clear object.
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GCTAEventCube::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GEventCube::free_members();
    this->GEvents::free_members();

    // Initialise members
    this->GEvents::init_members();
    this->GEventCube::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load CTA counts map from FITS file.
 *
 * @param[in] filename Counts map FITS filename to be loaded.
 *
 * It is assumed that the counts map resides in the primary extension of the
 * FITS file, the energy boundaries reside in the EBOUNDS extension and the
 * Good Time Intervals reside in the GTI extension.  The method clears the
 * object before loading, thus any events residing in the object before
 * loading will be lost.
 ***************************************************************************/
void GCTAEventCube::load(const std::string& filename)
{
    // Clear object
    clear();

    // Allocate FITS file
    GFits file;

    // Open counts map FITS file
    file.open(filename);

    // Get HDUs
    GFitsHDU* hdu_cntmap  = file.hdu("Primary");
    GFitsHDU* hdu_ebounds = file.hdu("EBOUNDS");
    GFitsHDU* hdu_gti     = file.hdu("GTI");

    // Load counts map
    read_cntmap(hdu_cntmap);

    // Load energy boundaries
    read_ebds(hdu_ebounds);

    // Load energy boundaries
    read_gti(hdu_gti);

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
 * This method provides the event attributes to the event bin. The event bin
 * is in fact physically stored in the event cube, and only a single event
 * bin is indeed allocated. This method sets up the pointers in the event
 * bin so that a client can easily access the information of individual bins
 * as if they were stored in an array.
 * The method returns a NULL pointer if the index is out of the valid range.
 *
 * @todo Static pointers could be set by init_members().
 * @todo Should we really return a NULL pointer in case that the index
 *       is not valid? Should we not better throw an exception? 
 ***************************************************************************/
GCTAEventBin* GCTAEventCube::pointer(int index)
{
//    #if defined(G_RANGE_CHECK)
//    if (index < 0 || index >= m_elements)
//        throw GException::out_of_range(G_POINTER, index, 0, m_elements-1);
//    #endif

    // Preset pointer with NULL
    GCTAEventBin* ptr = NULL;

    // Set pointer if index is in range
    if (index >=0 && index < m_elements) {

        // Set pointer to static element
        ptr = (GCTAEventBin*)&m_bin;

        // Get pixel and energy bin indices. Note that in GSkymap that holds
        // the counts, the energy axis is the most rapidely varying axis.
        int ipix = index / ebins();
        int ieng = index % ebins();

        // Set GEventBin pointers
        ptr->m_counts = &(m_counts[index]);
        ptr->m_time   = &m_time;
        ptr->m_energy = &(m_energies[ieng]);

        // Set GCTAEventBin pointers
        ptr->m_dir    = &(m_dirs[ipix]);
        ptr->m_pnt    = &m_pnt;
        ptr->m_rsp    = (m_obs != NULL) ?
                        (GCTAResponse*)((GObservation*)m_obs)->response() : NULL;
        ptr->m_omega  = &(m_omega[ipix]);
        ptr->m_ewidth = &(m_ewidth[ieng]);
        ptr->m_ontime = &m_ontime;

    } // endif: valid index

    // Return pointer
    return ptr;
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
 * @todo Implement GSkymap.clear(), GCTAEventBin.clear() methods
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
    m_omega    = NULL;
    m_energies = NULL;
    m_ewidth   = NULL;
    m_obs      = NULL;
    m_time.clear();
    m_ontime   = 0.0;
    m_pnt.clear();
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
    m_bin    = cube.m_bin;
    m_map    = cube.m_map;
    m_time   = cube.m_time;
    m_ontime = cube.m_ontime;
    m_pnt    = cube.m_pnt;
    m_ebds   = cube.m_ebds;
    m_gti    = cube.m_gti;
    m_obs    = cube.m_obs;

    // Set counter to copied skymap pixels
    m_counts = m_map.pixels();

    // Copy sky directions and solid angles
    if (cube.npix() > 0) {
        if (cube.m_dirs != NULL) {
            m_dirs = new GCTAInstDir[cube.npix()];
            for (int i = 0; i < cube.npix(); ++i)
                m_dirs[i] = cube.m_dirs[i];
        }
        if (cube.m_omega != NULL) {
            m_omega = new double[cube.npix()];
            for (int i = 0; i < cube.npix(); ++i)
                m_omega[i] = cube.m_omega[i];
        }
    }

    // Copy bin energies and widths
    if (cube.ebins() > 0) {
        if (cube.m_energies != NULL) {
            m_energies = new GEnergy[cube.ebins()];
            for (int i = 0; i < cube.ebins(); ++i)
                m_energies[i] = cube.m_energies[i];
        }
        if (cube.m_ewidth != NULL) {
            m_ewidth = new GEnergy[cube.ebins()];
            for (int i = 0; i < cube.ebins(); ++i)
                m_ewidth[i] = cube.m_ewidth[i];
        }
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
    if (m_ewidth   != NULL) delete [] m_ewidth;
    if (m_energies != NULL) delete [] m_energies;
    if (m_omega    != NULL) delete [] m_omega;
    if (m_dirs     != NULL) delete [] m_dirs;

    // Signal free pointers
    m_dirs     = NULL;
    m_omega    = NULL;
    m_energies = NULL;
    m_ewidth   = NULL;

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
 * @brief Read CTA counts map from HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which events are loaded.
 *
 * This method reads a CTA counts map from a FITS HDU. The counts map is
 * stored in a GSkymap object, and a pointer is set up to access the pixels
 * individually. Recall that skymap pixels are stored in the order
 * (ebin,ix,iy), i.e. the energy axis is the most rapidely varying axis,
 * while the counts map is stored in the order (ix,iy,ebin), i.e. the x
 * axis is the most rapidely varying axis.
 ***************************************************************************/
void GCTAEventCube::read_cntmap(GFitsHDU* hdu)
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

        // Set sky directions
        set_directions();

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy boundaries from HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which energy boundaries are loaded.
 *
 * Read the energy boundaries from the HDU.
 ***************************************************************************/
void GCTAEventCube::read_ebds(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Read energy boundaries
        m_ebds.read(hdu);

        // Set log mean energies and energy widths
        set_energies();

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read GTIs from HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which GTIs are loaded.
 *
 * Reads the Good Time Intervals from the GTI extension.
 ***************************************************************************/
void GCTAEventCube::read_gti(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Read Good Time Intervals
        m_gti.read(hdu);

        // Set time
        set_time();

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sky directions and solid angles of events cube.
 *
 * @exception GCTAException::no_sky
 *            No sky pixels found in event cube.
 *
 * This method computes the sky directions and solid angles for all event
 * cube pixels. Sky directions are stored in an array of GCTAInstDir objects
 * while solid angles are stored in units of sr in an array of double
 * precision variables. In addition, this method also sets the event cube
 * pixel pointer.
 ***************************************************************************/
void GCTAEventCube::set_directions(void)
{
    // Throw an error if we have no sky pixels
    if (npix() < 1)
        throw GCTAException::no_sky(G_SET_DIRECTIONS, "Every CTA event cube"
                                   " needs a definiton of the sky pixels.");

    // Set pixel pointer
    m_counts = m_map.pixels();

    // Delete old pixel directions and solid angles
    if (m_omega != NULL) delete [] m_omega;
    if (m_dirs  != NULL) delete [] m_dirs;

    // Set pixel directions and solid angles
    m_dirs  = new GCTAInstDir[npix()];
    m_omega = new double[npix()];
    for (int i = 0, iy = 0; iy < ny(); ++iy) {
        for (int ix = 0; ix < nx(); ++ix, ++i) {
            GSkyPixel pixel = GSkyPixel(double(ix), double(iy));
            m_dirs[i]       = GCTAInstDir(m_map.xy2dir(pixel));
            m_omega[i]      = m_map.omega(pixel);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log mean energies and energy widths of event cube.
 *
 * @exception GCTAException::no_ebds
 *            No energy boundaries found in event cube.
 *
 * This method computes the log mean energies and the energy widths of the
 * event cube. The log mean energies and energy widths are stored unit
 * independent in arrays of GEnergy objects.
 ***************************************************************************/
void GCTAEventCube::set_energies(void)
{
    // Throw an error if we have no energy bins
    if (ebins() < 1)
        throw GCTAException::no_ebds(G_SET_ENERGIES, "Every CTA event cube"
                             " needs a definiton of the energy boundaries.");

    // Delete old bin energies and energy widths
    if (m_ewidth   != NULL) delete [] m_ewidth;
    if (m_energies != NULL) delete [] m_energies;

    // Setup bin energies and energy widths
    m_energies = new GEnergy[ebins()];
    m_ewidth   = new GEnergy[ebins()];
    for (int i = 0; i < ebins(); ++i) {
        m_energies[i] = m_ebds.elogmean(i);
        m_ewidth[i]   = m_ebds.emax(i) - m_ebds.emin(i);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set mean event time and ontime of event cube.
 *
 * @exception GCTAException::no_gti
 *            No Good Time Intervals found in event cube.
 *
 * This method computes the mean event time and the ontime of the event
 * cube. The mean event time is the average between the start and the stop
 * time. The ontime is the sum of all Good Time Intervals.
 *
 * @todo Could add a more sophisticated mean event time computation that
 *       weights by the length of the GTIs, yet so far we do not really use
 *       the mean event time, hence there is no rush to implement this.
 ***************************************************************************/
void GCTAEventCube::set_time(void)
{
    // Throw an error if GTI is empty
    if (m_gti.size() < 1)
        throw GCTAException::no_gti(G_SET_TIME, "Every CTA event cube needs"
                  " associated GTIs to allow the computation of the ontime.");

    // Compute mean time
    m_time = 0.5 * (m_gti.tstart() + m_gti.tstop());

    // Set ontime
    m_ontime = m_gti.ontime();

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
    os << " Number of pixels ..........: " << cube.npix() << std::endl;
    os << " Number of energy bins .....: " << cube.ebins() << std::endl;
    os << " Number of events ..........: " << cube.number() << std::endl;
    os << cube.m_ebds;

    // Return output stream
    return os;
}
