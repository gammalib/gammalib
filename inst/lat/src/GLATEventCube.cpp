/***************************************************************************
 *                GLATEventCube.cpp  -  LAT event cube class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATEventCube.cpp
 * @brief GLATEventCube class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATEventCube.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_POINTER                               "GLATEventCube::pointer(int)"
#define G_DIFFNAME                            "GLATEventCube::diffname(int&)"
#define G_DIFFRSP                              "GLATEventCube::diffrsp(int&)"
#define G_READ_SRCMAP               "GLATEventCube::read_srcmap(GFitsImage*)"
#define G_SET_DIRECTIONS                    "GLATEventCube::set_directions()"
#define G_SET_ENERGIES                        "GLATEventCube::set_energies()"
#define G_SET_TIME                                "GLATEventCube::set_time()"

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
GLATEventCube::GLATEventCube(void) : GEventCube()
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
GLATEventCube::GLATEventCube(const GLATEventCube& cube) : GEventCube(cube)
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
GLATEventCube::~GLATEventCube(void)
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
GLATEventCube& GLATEventCube::operator= (const GLATEventCube& cube)
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
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATEventCube::clear(void)
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
GLATEventCube* GLATEventCube::clone(void) const
{
    return new GLATEventCube(*this);
}


/***********************************************************************//**
 * @brief Load LAT counts map from FITS file.
 *
 * @param[in] filename Counts map FITS filename to be loaded.
 *
 * It is assumed that the counts map resides in the primary extension of the
 * FITS file, the energy boundaries reside in the EBOUNDS extension and the
 * Good Time Intervals reside in the GTI extension.  The method clears the
 * object before loading, thus any events residing in the object before
 * loading will be lost.
 ***************************************************************************/
void GLATEventCube::load(const std::string& filename)
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

    // Load additional source maps
    for (int i = 1; i < file.size(); ++i) {
        if (file.hdu(i)->exttype() == GFitsHDU::HT_IMAGE) {
            GFitsImage* hdu_srcmap = file.image(i);
            read_srcmap(hdu_srcmap);
        }
    }

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get pointer to element
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
GLATEventBin* GLATEventCube::pointer(int index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_elements)
        throw GException::out_of_range(G_POINTER, index, 0, m_elements-1);
    #endif

    // Get pixel and energy bin indices. Note that in GSkymap that holds
    // the counts, the energy axis is the most rapidely varying axis.
    m_bin.m_index = index;
    m_bin.m_ipix  = index / ebins();
    m_bin.m_ieng  = index % ebins();

    // Set pointers
    m_bin.m_cube   = this;
    m_bin.m_counts = &(m_counts[index]);
    m_bin.m_energy = &(m_energies[m_bin.m_ieng]);
    m_bin.m_time   = &m_time;
    m_bin.m_dir    = &(m_dirs[m_bin.m_ipix]);
    m_bin.m_omega  = &(m_omega[m_bin.m_ipix]);
    m_bin.m_ewidth = &(m_ewidth[m_bin.m_ieng]);
    m_bin.m_ontime = &m_ontime;

    // Return pointer
    return &m_bin;
}


/***********************************************************************//**
 * @brief Return number of events in cube
 ***************************************************************************/
int GLATEventCube::number(void) const
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


/***********************************************************************//**
 * @brief Print event cube information
 ***************************************************************************/
std::string GLATEventCube::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATEventCube ===\n");
    result.append(parformat("Number of elements")+str(size())+"\n");
    result.append(parformat("Number of pixels"));
    result.append(str(m_map.nx())+" x "+str(m_map.ny())+"\n");
    result.append(parformat("Number of energy bins")+str(ebins())+"\n");
    result.append(parformat("Number of events")+str(number())+"\n");

    // Append source maps
    result.append(parformat("Number of source maps")+str(m_srcmap.size()));
    for (int i = 0; i < m_srcmap.size(); ++i) {
        result.append("\n"+parformat(" "+m_srcmap_names[i]));
        result.append(str(m_srcmap[i]->nx()));
        result.append(" x ");
        result.append(str(m_srcmap[i]->ny()));
        result.append(" x ");
        result.append(str(m_srcmap[i]->nmaps()));
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return name of diffuse model
 *
 * @param[in] index Diffuse model index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Model index out of valid range.
 *
 * Returns name of diffuse model.
 ***************************************************************************/
std::string GLATEventCube::diffname(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= ndiffrsp())
        throw GException::out_of_range(G_DIFFNAME, index, 0, ndiffrsp()-1);
    #endif

    // Return
    return m_srcmap_names[index];
}


/***********************************************************************//**
 * @brief Return diffuse response map
 *
 * @param[in] index Diffuse model index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Model index out of valid range.
 *
 * Returns pointer to diffuse model sky map.
 ***************************************************************************/
GSkymap* GLATEventCube::diffrsp(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= ndiffrsp())
        throw GException::out_of_range(G_DIFFRSP, index, 0, ndiffrsp()-1);
    #endif

    // Return
    return m_srcmap[index];
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 *
 * @todo Implement GSkymap.clear() method
 ***************************************************************************/
void GLATEventCube::init_members(void)
{
    // Initialise members
    m_bin.clear();
    //m_map.clear();
    m_map      = GSkymap();
    m_counts   = NULL;
    m_dirs     = NULL;
    m_omega    = NULL;
    m_energies = NULL;
    m_ewidth   = NULL;
    m_time.clear();
    m_ontime   = 0.0;
    m_ebds.clear();
    m_gti.clear();
    m_srcmap.clear();
    m_srcmap_names.clear();
    m_enodes.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube GLATEventCube members which should be copied.
 ***************************************************************************/
void GLATEventCube::copy_members(const GLATEventCube& cube)
{
    // Copy LAT specific attributes
    m_bin          = cube.m_bin;
    m_map          = cube.m_map;
    m_time         = cube.m_time;
    m_ontime       = cube.m_ontime;
    m_ebds         = cube.m_ebds;
    m_gti          = cube.m_gti;
    m_srcmap       = cube.m_srcmap;
    m_srcmap_names = cube.m_srcmap_names;
    m_enodes       = cube.m_enodes;

    // Set counter to copied skymap pixels
    m_counts = m_map.pixels();

    // Copy sky directions and solid angles
    if (cube.npix() > 0) {
        if (cube.m_dirs != NULL) {
            m_dirs = new GLATInstDir[cube.npix()];
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
void GLATEventCube::free_members(void)
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
 * @brief Read LAT counts map from HDU.
 *
 * @param[in] hdu Pointer to image HDU.
 *
 * This method reads a LAT counts map from a FITS image. The counts map is
 * stored in a GSkymap object, and a pointer is set up to access the pixels
 * individually. Recall that skymap pixels are stored in the order
 * (ebin,ix,iy), i.e. the energy axis is the most rapidely varying axis,
 * while the counts map is stored in the order (ix,iy,ebin), i.e. the x
 * axis is the most rapidely varying axis.
 ***************************************************************************/
void GLATEventCube::read_cntmap(GFitsImage* hdu)
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

        // Set pixel pointer
        m_counts = m_map.pixels();

        // Set sky directions
        set_directions();

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read LAT source map from HDU.
 *
 * @param[in] hdu Pointer to image HDU.
 *
 * @exception GLATException::wcs_incompatible
 *            Source map not compatible with sky map
 *
 * This method reads a LAT source map from a FITS image. The source map is
 * stored in a GSkymap object. Recall that skymap pixels are stored in the
 * order (ebin,ix,iy), i.e. the energy axis is the most rapidely varying
 * axis.
 *
 * @todo Add WCS comparison method to check whether source map is consistent
 *       with counts map.
 ***************************************************************************/
void GLATEventCube::read_srcmap(GFitsImage* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Allocate skymap
        GSkymap* map = new GSkymap;

        // Read skymap
        map->read(hdu);

        // Check that source map is consistent with counts map (to be replaced
        // later by WCS comparison method)
        if (m_map.nx() != map->nx() ||
            m_map.ny() != map->ny())
            throw GLATException::wcs_incompatible(G_READ_SRCMAP, hdu->extname());

        // Check that source map has required number of energy bins
        if (m_map.nmaps()+1 != map->nmaps())
            throw GLATException::wcs_incompatible(G_READ_SRCMAP, hdu->extname());

        // Append source map to list of maps
        m_srcmap.push_back(map);
        m_srcmap_names.push_back(hdu->extname());

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
void GLATEventCube::read_ebds(GFitsTable* hdu)
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
void GLATEventCube::read_gti(GFitsTable* hdu)
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
 * @exception GLATException::no_sky
 *            No sky pixels found in event cube.
 *
 * This method computes the sky directions and solid angles for all event
 * cube pixels. Sky directions are stored in an array of GLATInstDir objects
 * while solid angles are stored in units of sr in a double precision array.
 ***************************************************************************/
void GLATEventCube::set_directions(void)
{
    // Throw an error if we have no sky pixels
    if (npix() < 1)
        throw GLATException::no_sky(G_SET_DIRECTIONS, "Every LAT event cube"
                                   " needs a definiton of the sky pixels.");

    // Delete old pixel directions and solid angles
    if (m_omega != NULL) delete [] m_omega;
    if (m_dirs  != NULL) delete [] m_dirs;

    // Set pixel directions and solid angles
    m_dirs  = new GLATInstDir[npix()];
    m_omega = new double[npix()];
    for (int i = 0, iy = 0; iy < ny(); ++iy) {
        for (int ix = 0; ix < nx(); ++ix, ++i) {
            GSkyPixel pixel = GSkyPixel(double(ix), double(iy));
            m_dirs[i]       = GLATInstDir(m_map.xy2dir(pixel));
            m_omega[i]      = m_map.omega(pixel);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log mean energies and energy widths of event cube.
 *
 * @exception GLATException::no_ebds
 *            No energy boundaries found in event cube.
 *
 * This method computes the log mean energies and the energy widths of the
 * event cube. The log mean energies and energy widths are stored unit
 * independent in arrays of GEnergy objects.
 ***************************************************************************/
void GLATEventCube::set_energies(void)
{
    // Throw an error if we have no energy bins
    if (ebins() < 1)
        throw GLATException::no_ebds(G_SET_ENERGIES, "Every LAT event cube"
                             " needs a definiton of the energy boundaries.");

    // Delete old bin energies and energy widths
    if (m_ewidth   != NULL) delete [] m_ewidth;
    if (m_energies != NULL) delete [] m_energies;

    // Clear energy nodes
    m_enodes.clear();

    // Setup bin energies, energy widths and energy nodes
    m_energies = new GEnergy[ebins()];
    m_ewidth   = new GEnergy[ebins()];
    for (int i = 0; i < ebins(); ++i) {
        m_energies[i] = m_ebds.elogmean(i);
        m_ewidth[i]   = m_ebds.emax(i) - m_ebds.emin(i);
        m_enodes.append(log10(m_ebds.emin(i).MeV()));
    }
    m_enodes.append(log10(m_ebds.emax(ebins()-1).MeV()));
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set mean event time and ontime of event cube.
 *
 * @exception GLATException::no_gti
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
void GLATEventCube::set_time(void)
{
    // Throw an error if GTI is empty
    if (m_gti.size() < 1)
        throw GLATException::no_gti(G_SET_TIME, "Every LAT event cube needs"
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
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
