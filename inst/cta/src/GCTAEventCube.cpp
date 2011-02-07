/***************************************************************************
 *                GCTAEventCube.cpp  -  CTA event cube class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
#include "GTools.hpp"
#include "GCTAException.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAEventBin.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponse.hpp"
#include "GSkymap.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsImage.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAXIS                                   "GCTAEventCube::naxis(int)"
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
 * @brief Void constructor
 ***************************************************************************/
GCTAEventCube::GCTAEventCube(void) : GEventCube()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky map constructor
 *
 * @param[in] map Sky map.
 *
 * Construct instance of events cube from sky map.
 ***************************************************************************/
GCTAEventCube::GCTAEventCube(const GSkymap& map) : GEventCube()
{
    // Initialise members
    init_members();

    // Set sky map
    m_map = map;

    // Set sky directions
    set_directions();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Event cube.
 ***************************************************************************/
GCTAEventCube::GCTAEventCube(const GCTAEventCube& cube) : GEventCube(cube)
{
    // Initialise members
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
 * @param[in] cube Event cube.
 ***************************************************************************/
GCTAEventCube& GCTAEventCube::operator= (const GCTAEventCube& cube)
{
    // Execute only if object is not identical
    if (this != &cube) {

        // Copy base class members
        this->GEventCube::operator=(cube);

        // Free members
        free_members();

        // Initialise members
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
 * @brief Clear instance
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
 * @brief Clone instance
 ***************************************************************************/
GCTAEventCube* GCTAEventCube::clone(void) const
{
    return new GCTAEventCube(*this);
}


/***********************************************************************//**
 * @brief Return number of bins in event cube
 ***************************************************************************/
int GCTAEventCube::size(void) const
{
    // Compute number of bins
    int nbins = m_map.npix() * m_map.nmaps(); 

    // Return number of bins
    return nbins;
}


/***********************************************************************//**
 * @brief Return dimension of event cube
 ***************************************************************************/
int GCTAEventCube::dim(void) const
{
    // Compute dimension from sky map
    int dim = (m_map.nmaps() > 1) ? 3 : 2;

    // Return dimension
    return dim;
}


/***********************************************************************//**
 * @brief Return number of bins in axis
 *
 * @param[in] axis Axis.
 *
 * @exception GException::out_of_range
 *            Axis is out of range.
 *
 * Returns the number of bins along a given event cube axis.
 ***************************************************************************/
int GCTAEventCube::naxis(int axis) const
{
    // Optionally check if the axis is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= dim())
        throw GException::out_of_range(G_NAXIS, axis, 0, dim()-1);
    #endif

    // Set result
    int naxis = 0;
    switch (axis) {
    case 0:
        naxis = m_map.nx();
        break;
    case 1:
        naxis = m_map.ny();
        break;
    case 2:
        naxis = m_map.npix();
        break;
    }

    // Return result
    return naxis;
}


/***********************************************************************//**
 * @brief Load CTA counts map from FITS file
 *
 * @param[in] filename FITS file name of counts map.
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

    // Open counts map FITS file
    GFits file(filename);

    // Get HDUs
    GFitsImage* hdu_cntmap  = file.image("Primary");
    GFitsTable* hdu_ebounds = file.table("EBOUNDS");
    GFitsTable* hdu_gti     = file.table("GTI");

    // Load counts map
    read_cntmap(hdu_cntmap);

    // Load energy boundaries
    read_ebds(hdu_ebounds);

    // Load GTIs
    read_gti(hdu_gti);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA event cube into FITS file.
 *
 * @param[in] file FITS file.
 ***************************************************************************/
void GCTAEventCube::write(GFits* file) const
{
    // Write cube
    m_map.write(file);

    // Write energy boundaries
    m_ebds.write(file);

    // Write Good Time intervals
    m_gti.write(file);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to an event
 *
 * @param[in] index Event index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Event index not in valid range.
 *
 * This method provides the event attributes to the event bin. The event bin
 * is in fact physically stored in the event cube, and only a single event
 * bin is indeed allocated. This method sets up the pointers in the event
 * bin so that a client can easily access the information of individual bins
 * as if they were stored in an array.
 ***************************************************************************/
GCTAEventBin* GCTAEventCube::pointer(int index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_POINTER, index, 0, size()-1);
    #endif

    // Get pixel and energy bin indices.
    int ipix = index % npix();
    int ieng = index / npix();

    // Set pointers
    m_bin.m_counts = &(m_map.pixels()[index]);
    m_bin.m_energy = &(m_energies[ieng]);
    m_bin.m_time   = &m_time;
    m_bin.m_dir    = &(m_dirs[ipix]);
    m_bin.m_omega  = &(m_omega[ipix]);
    m_bin.m_ewidth = &(m_ewidth[ieng]);
    m_bin.m_ontime = &m_ontime;

    // Return pointer
    return &m_bin;
}


/***********************************************************************//**
 * @brief Return number of events in cube
 ***************************************************************************/
int GCTAEventCube::number(void) const
{
    // Initialise result
    double number = 0.0;

    // Get pointer on skymap pixels
    double* pixels = m_map.pixels();

    // Sum event cube
    if (size() > 0 && pixels != NULL) {
        for (int i = 0; i < size(); ++i)
            number += pixels[i];
    }

    // Return
    return int(number+0.5);
}


/***********************************************************************//**
 * @brief Print event cube information
 ***************************************************************************/
std::string GCTAEventCube::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAEventCube ===");
    result.append("\n"+parformat("Number of elements")+str(size()));
    result.append("\n"+parformat("Number of pixels")+str(npix()));
    result.append("\n"+parformat("Number of energy bins")+str(ebins()));
    result.append("\n"+parformat("Number of events")+str(number()));

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAEventCube::init_members(void)
{
    // Initialise members
    m_bin.clear();
    m_map.clear();
    m_time.clear();
    m_ebds.clear();
    m_gti.clear();
    m_dirs     = NULL;
    m_omega    = NULL;
    m_energies = NULL;
    m_ewidth   = NULL;
    m_ontime   = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Event cube.
 ***************************************************************************/
void GCTAEventCube::copy_members(const GCTAEventCube& cube)
{
    // Copy CTA specific attributes
    m_bin    = cube.m_bin;
    m_map    = cube.m_map;
    m_time   = cube.m_time;
    m_ontime = cube.m_ontime;
    m_ebds   = cube.m_ebds;
    m_gti    = cube.m_gti;

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
 * @brief Read CTA counts map from HDU.
 *
 * @param[in] hdu Pointer to image HDU.
 *
 * This method reads a CTA counts map from a FITS HDU. The counts map is
 * stored in a GSkymap object.
 ***************************************************************************/
void GCTAEventCube::read_cntmap(GFitsImage* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Load counts map as sky map
        m_map.read(hdu);

        // Set sky directions
        set_directions();

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy boundaries from HDU.
 *
 * @param[in] hdu Pointer to energy boundaries table.
 *
 * Read the energy boundaries from the HDU.
 ***************************************************************************/
void GCTAEventCube::read_ebds(GFitsTable* hdu)
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
 * @param[in] hdu Pointer to GTI table.
 *
 * Reads the Good Time Intervals from the GTI extension.
 ***************************************************************************/
void GCTAEventCube::read_gti(GFitsTable* hdu)
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
 * precision variables.
 ***************************************************************************/
void GCTAEventCube::set_directions(void)
{
    // Throw an error if we have no sky pixels
    if (npix() < 1)
        throw GCTAException::no_sky(G_SET_DIRECTIONS, "Every CTA event cube"
                                   " needs a definiton of the sky pixels.");

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
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
