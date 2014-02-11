/***************************************************************************
 *             GLATEventCube.cpp - Fermi/LAT event cube class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATEventCube.cpp
 * @brief Fermi/LAT event cube class definition
 * @author Juergen Knoedlseder
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
#define G_NAXIS                                   "GLATEventCube::naxis(int)"
#define G_DIFFNAME                            "GLATEventCube::diffname(int&)"
#define G_DIFFRSP                              "GLATEventCube::diffrsp(int&)"
#define G_READ_SRCMAP               "GLATEventCube::read_srcmap(GFitsImage&)"
#define G_SET_DIRECTIONS                    "GLATEventCube::set_directions()"
#define G_SET_ENERGIES                        "GLATEventCube::set_energies()"
#define G_SET_TIMES                              "GLATEventCube::set_times()"
#define G_SET_BIN                              "GLATEventCube::set_bin(int&)"

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
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube LAT event cube.
 ***************************************************************************/
GLATEventCube::GLATEventCube(const GLATEventCube& cube) : GEventCube(cube)
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
 * @param[in] cube LAT event cube.
 * @return LAT event cube.
 ***************************************************************************/
GLATEventCube& GLATEventCube::operator= (const GLATEventCube& cube)
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


/***********************************************************************//**
 * @brief Event bin access operator
 *
 * @param[in] index Event index [0,...,size()-1].
 * @return Pointer to event bin.
 *
 * Returns pointer to an event bin.
 ***************************************************************************/
GLATEventBin* GLATEventCube::operator[](const int& index)
{
    // Set event bin
    set_bin(index);

    // Return pointer
    return (&m_bin);
}


/***********************************************************************//**
 * @brief Event bin access operator (const version)
 *
 * @param[in] index Event index [0,...,size()-1].
 * @return Pointer to event bin.
 *
 * Returns pointer to an event bin.
 ***************************************************************************/
const GLATEventBin* GLATEventCube::operator[](const int& index) const
{
    // Set event bin (circumvent const correctness)
    const_cast<GLATEventCube*>(this)->set_bin(index);

    // Return pointer
    return (&m_bin);
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
 *
 * @return Deep copy of LAT event cube.
 ***************************************************************************/
GLATEventCube* GLATEventCube::clone(void) const
{
    return new GLATEventCube(*this);
}


/***********************************************************************//**
 * @brief Return number of bins in event cube
 *
 * @return Number of bins in event cube.
 ***************************************************************************/
int GLATEventCube::size(void) const
{
    // Compute number of bins
    int nbins = m_map.npix() * m_map.nmaps(); 

    // Return number of bins
    return nbins;
}


/***********************************************************************//**
 * @brief Return dimension of event cube
 *
 * @return Number of dimensions in event cube (either 2 or 3).
 ***************************************************************************/
int GLATEventCube::dim(void) const
{
    // Compute dimension from sky map
    int dim = (m_map.nmaps() > 1) ? 3 : 2;

    // Return dimension
    return dim;
}


/***********************************************************************//**
 * @brief Return number of bins in axis
 *
 * @param[in] axis Axis [0,...,dim()-1]
 * @return Number of bins in specified axis
 *
 * @exception GException::out_of_range
 *            Axis is out of range.
 *
 * Returns the number of bins along a given event cube @p axis.
 ***************************************************************************/
int GLATEventCube::naxis(const int& axis) const
{
    // Optionally check if the axis is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= dim()) {
        throw GException::out_of_range(G_NAXIS, "LAT event cube axis", axis, dim());
    }
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
 * @brief Load LAT event cube from FITS file
 *
 * @param[in] filename FITS file name.
 ***************************************************************************/
void GLATEventCube::load(const std::string& filename)
{
    // Clear object
    clear();

    // Open FITS file
    GFits file(filename);

    // Read counts map
    read(file);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save LAT event cube into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing FITS file? (default=false)
 *
 * Save the LAT event cube into FITS file.
 ***************************************************************************/
void GLATEventCube::save(const std::string& filename,
                         const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write event cube
    write(fits);
    
    // Save FITS file
    fits.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read LAT event cube from FITS file.
 *
 * @param[in] fits FITS file.
 *
 * It is assumed that the counts map resides in the primary extension of the
 * FITS file, the energy boundaries reside in the EBOUNDS extension and the
 * Good Time Intervals reside in the GTI extension.  The method clears the
 * object before loading, thus any events residing in the object before
 * loading will be lost.
 ***************************************************************************/
void GLATEventCube::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get HDUs
    const GFitsImage& hdu_cntmap  = *fits.image("Primary");
    const GFitsTable& hdu_ebounds = *fits.table("EBOUNDS");
    const GFitsTable& hdu_gti     = *fits.table("GTI");

    // Load counts map
    read_cntmap(hdu_cntmap);

    // Load energy boundaries
    read_ebds(hdu_ebounds);

    // Load GTIs
    read_gti(hdu_gti);

    // Load additional source maps
    for (int i = 1; i < fits.size(); ++i) {
        if (fits.at(i)->exttype() == GFitsHDU::HT_IMAGE) {
            const GFitsImage& hdu_srcmap = *fits.image(i);
            read_srcmap(hdu_srcmap);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write LAT event cube into FITS file
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GLATEventCube::write(GFits& fits) const
{
    // Write cube
    m_map.write(fits);

    // Write energy boundaries
    ebounds().write(fits);

    // Write Good Time intervals
    gti().write(fits);

    // Write additional source maps
    for (int i = 0; i < m_srcmap.size(); ++i) {
        m_srcmap[i]->write(fits);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of events in cube
 *
 * @return Total number of events in event cube.
 ***************************************************************************/
int GLATEventCube::number(void) const
{
    // Initialise result
    double number = 0.0;

    // Get pointer on skymap pixels
    const double* pixels = m_map.pixels();

    // Sum event cube
    if (size() > 0 && pixels != NULL) {
        for (int i = 0; i < size(); ++i) {
            number += pixels[i];
        }
    }

    // Return
    return int(number+0.5);
}


/***********************************************************************//**
 * @brief Set event cube from sky map
 ***************************************************************************/
void GLATEventCube::map(const GSkymap& map)
{
    // Store sky map
    m_map = map;

    // Compute sky directions
    set_directions();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print event cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing event cube information.
 ***************************************************************************/
std::string GLATEventCube::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATEventCube ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of elements") +
                      gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of pixels"));
        result.append(gammalib::str(m_map.nx()) +
                      " x " +
                      gammalib::str(m_map.ny()));
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(ebins()));
        result.append("\n"+gammalib::parformat("Number of events") +
                      gammalib::str(number()));

        // Append time interval
        result.append("\n"+gammalib::parformat("Time interval"));
        if (gti().size() > 0) {
            result.append(gammalib::str(tstart().secs()) +
                          " - " +
                          gammalib::str(tstop().secs())+" sec");
        }
        else {
            result.append("not defined");
        }

        // Append energy range
        result.append("\n"+gammalib::parformat("Energy range"));
        if (ebounds().size() > 0) {
            result.append(emin().print()+" - "+emax().print(chatter));
        }
        else {
            result.append("not defined");
        }

        // Append detailed information
        GChatter reduced_chatter = gammalib::reduce(chatter);
        if (reduced_chatter > SILENT) {

            // Append sky projection
            if (m_map.projection() != NULL) {
                result.append("\n"+m_map.projection()->print(reduced_chatter));
            }

            // Append source maps
            result.append("\n"+gammalib::parformat("Number of source maps"));
            result.append(gammalib::str(m_srcmap.size()));
            for (int i = 0; i < m_srcmap.size(); ++i) {
                result.append("\n"+gammalib::parformat(" "+m_srcmap_names[i]));
                result.append(gammalib::str(m_srcmap[i]->nx()));
                result.append(" x ");
                result.append(gammalib::str(m_srcmap[i]->ny()));
                result.append(" x ");
                result.append(gammalib::str(m_srcmap[i]->nmaps()));
            }

        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return name of diffuse model
 *
 * @param[in] index Diffuse model index [0,...,ndiffrsp()-1].
 * @return Name of diffuse model.
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
    if (index < 0 || index >= ndiffrsp()) {
        throw GException::out_of_range(G_DIFFNAME, index, 0, ndiffrsp()-1);
    }
    #endif

    // Return
    return m_srcmap_names[index];
}


/***********************************************************************//**
 * @brief Return diffuse response map
 *
 * @param[in] index Diffuse model index [0,...,ndiffrsp()-1].
 * @return Pointer to diffuse response map.
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
    if (index < 0 || index >= ndiffrsp()) {
        throw GException::out_of_range(G_DIFFRSP, index, 0, ndiffrsp()-1);
    }
    #endif

    // Return
    return m_srcmap[index];
}


/***********************************************************************//**
 * @brief Computes the maximum radius (in degrees) around a given source
 *        direction that fits spatially into the event cube
 *
 * @param[in] srcDir Source direction.
 * @return Maximum radius in degrees that fully fits into event cube.
 *
 * By computing the sky directions of the event cube boundaries, the maximum
 * radius is computed that fits fully within the event cube. This method is
 * used for PSF normalization.
 ***************************************************************************/
double GLATEventCube::maxrad(const GSkyDir& srcDir) const
{
    // Initialise radius
    double radius = 180.0;

    // Move along upper edge in longitude
    int iy = 0;
    for (int ix = 0; ix < nx(); ++ix) {
        GSkyPixel pixel    = GSkyPixel(double(ix), double(iy));
        double    distance = m_map.pix2dir(pixel).dist_deg(srcDir);
        if (distance < radius) {
            radius = distance;
        }
    }

    // Move along lower edge in longitude
    iy = ny()-1;
    for (int ix = 0; ix < nx(); ++ix) {
        GSkyPixel pixel    = GSkyPixel(double(ix), double(iy));
        double    distance = m_map.pix2dir(pixel).dist_deg(srcDir);
        if (distance < radius) {
            radius = distance;
        }
    }

    // Move along left edge in latitude
    int ix = 0;
    for (int iy = 0; iy < ny(); ++iy) {
        GSkyPixel pixel    = GSkyPixel(double(ix), double(iy));
        double    distance = m_map.pix2dir(pixel).dist_deg(srcDir);
        if (distance < radius) {
            radius = distance;
        }
    }

    // Move along right edge in latitude
    ix = nx()-1;
    for (int iy = 0; iy < ny(); ++iy) {
        GSkyPixel pixel    = GSkyPixel(double(ix), double(iy));
        double    distance = m_map.pix2dir(pixel).dist_deg(srcDir);
        if (distance < radius) {
            radius = distance;
        }
    }

    // Return radius
    return radius;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATEventCube::init_members(void)
{
    // Initialise members
    m_bin.clear();
    m_map.clear();
    m_time.clear();
    m_srcmap.clear();
    m_srcmap_names.clear();
    m_enodes.clear();
    m_dirs.clear();
    m_solidangle.clear();
    m_energies.clear(); 
    m_ewidth.clear(); 
    m_ontime = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Fermi/LAT event cube.
 ***************************************************************************/
void GLATEventCube::copy_members(const GLATEventCube& cube)
{
    // Copy Fermi/LAT specific attributes
    m_bin          = cube.m_bin;
    m_map          = cube.m_map;
    m_time         = cube.m_time;
    m_ontime       = cube.m_ontime;
    m_enodes       = cube.m_enodes;
    m_dirs         = cube.m_dirs;
    m_solidangle   = cube.m_solidangle;
    m_energies     = cube.m_energies;
    m_ewidth       = cube.m_ewidth;

    // Clone source maps and copy source map names
    m_srcmap.clear();
    for (int i = 0; i < cube.m_srcmap.size(); ++i) {
        m_srcmap.push_back(cube.m_srcmap[i]->clone());
    }
    m_srcmap_names = cube.m_srcmap_names;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEventCube::free_members(void)
{
    // Release source maps
    for (int i = 0; i < m_srcmap.size(); ++i) {
        if (m_srcmap[i] != NULL) delete m_srcmap[i];
        m_srcmap[i] = NULL;
    }
    m_srcmap.clear();
    m_srcmap_names.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Fermi/LAT counts map from HDU.
 *
 * @param[in] hdu Image HDU.
 *
 * This method reads a Fermi/LAT counts map from a FITS image. The counts map
 * is stored in a GSkymap object, and a pointer is set up to access the
 * pixels individually. Recall that skymap pixels are stored in the order
 * (ix,iy,ebin).
 ***************************************************************************/
void GLATEventCube::read_cntmap(const GFitsImage& hdu)
{
    // Load counts map as sky map
    m_map.read(hdu);

    // Set sky directions
    set_directions();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read LAT source map from HDU.
 *
 * @param[in] hdu Image HDU.
 *
 * @exception GLATException::wcs_incompatible
 *            Source map not compatible with sky map
 *
 * This method reads a LAT source map from a FITS image. The source map is
 * stored in a GSkymap object and is given in units of counts/pixel/MeV.
 ***************************************************************************/
void GLATEventCube::read_srcmap(const GFitsImage& hdu)
{
    // Allocate skymap
    GSkymap* map = new GSkymap;

    // Read skymap
    map->read(hdu);

    // Check that source map sky projection is consistent with counts
    // map sky projection
    if (*(m_map.projection()) != *(map->projection())) {
        throw GLATException::wcs_incompatible(G_READ_SRCMAP, hdu.extname());
    }

    // Check that source map dimension is consistent with counts map
    // dimension
    if (m_map.nx() != map->nx() ||
        m_map.ny() != map->ny()) {
        throw GLATException::wcs_incompatible(G_READ_SRCMAP, hdu.extname());
    }

    // Check that source map has required number of energy bins
    if (m_map.nmaps()+1 != map->nmaps()) {
        throw GLATException::wcs_incompatible(G_READ_SRCMAP, hdu.extname());
    }

    // Append source map to list of maps
    m_srcmap.push_back(map);
    m_srcmap_names.push_back(hdu.extname());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy boundaries from HDU.
 *
 * @param[in] hdu Energy boundaries table.
 *
 * Read the energy boundaries from the HDU.
 *
 * @todo Energy bounds read method should take const GFitsTable* as argument
 ***************************************************************************/
void GLATEventCube::read_ebds(const GFitsTable& hdu)
{
    // Read energy boundaries
    m_ebounds.read(hdu);

    // Set log mean energies and energy widths
    set_energies();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read GTIs from HDU.
 *
 * @param[in] hdu GTI table.
 *
 * Reads the Good Time Intervals from the GTI extension.
 *
 * @todo GTI read method should take const GFitsTable* as argument
 ***************************************************************************/
void GLATEventCube::read_gti(const GFitsTable& hdu)
{
    // Read Good Time Intervals
    m_gti.read(hdu);

    // Set time
    set_times();

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
    if (npix() < 1) {
        throw GLATException::no_sky(G_SET_DIRECTIONS, "Every LAT event cube"
                                   " needs a definition of the sky pixels.");
    }

    // Clear old pixel directions and solid angle
    m_dirs.clear();
    m_solidangle.clear();

    // Reserve space for pixel directions and solid angles
    m_dirs.reserve(npix());
    m_solidangle.reserve(npix());

    // Set pixel directions and solid angles
    for (int iy = 0; iy < ny(); ++iy) {
        for (int ix = 0; ix < nx(); ++ix) {
            GSkyPixel pixel = GSkyPixel(double(ix), double(iy));
            m_dirs.push_back(GLATInstDir(m_map.pix2dir(pixel)));
            m_solidangle.push_back(m_map.solidangle(pixel));
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
    if (ebins() < 1) {
        throw GLATException::no_ebds(G_SET_ENERGIES, "Every LAT event cube"
                             " needs a definition of the energy boundaries.");
    }

    // Clear old bin energies and energy widths
    m_energies.clear();
    m_ewidth.clear();
    m_enodes.clear();

    // Reserve space for bin energies and energy widths
    m_energies.reserve(ebins());
    m_ewidth.reserve(ebins());

    // Setup bin energies, energy widths and energy nodes
    for (int i = 0; i < ebins(); ++i) {
        m_energies.push_back(ebounds().elogmean(i));
        m_ewidth.push_back(ebounds().emax(i) -  ebounds().emin(i));
        m_enodes.append(log10(ebounds().emin(i).MeV()));
    }
    m_enodes.append(log10(ebounds().emax(ebins()-1).MeV()));
    
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
void GLATEventCube::set_times(void)
{
    // Throw an error if GTI is empty
    if (m_gti.size() < 1) {
        throw GLATException::no_gti(G_SET_TIMES, "Every LAT event cube needs"
                  " associated GTIs to allow the computation of the ontime.");
    }

    // Compute mean time
    m_time = 0.5 * (gti().tstart() + gti().tstop());

    // Set ontime
    m_ontime = m_gti.ontime();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set event bin
 *
 * @param[in] index Event index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Event index is outside valid range.
 * @exception GLATException::no_energies
 *            Energy vectors have not been set up.
 * @exception GLATException::no_dirs
 *            Sky directions and solid angles vectors have not been set up.
 *
 * This method provides the event attributes to the event bin. The event bin
 * is in fact physically stored in the event cube, and only a single event
 * bin is indeed allocated. This method sets up the pointers in the event
 * bin so that a client can easily access the information of individual bins
 * as if they were stored in an array.
 ***************************************************************************/
void GLATEventCube::set_bin(const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET_BIN, index, 0, size()-1);
    }
    #endif

    // Check for the existence of energies and energy widths
    if (m_energies.size() != ebins() || m_ewidth.size() != ebins()) {
        throw GLATException::no_energies(G_SET_BIN);
    }

    // Check for the existence of sky directions and solid angles
    if (m_dirs.size() != npix() || m_solidangle.size() != npix()) {
        throw GLATException::no_dirs(G_SET_BIN);
    }

    // Get pixel and energy bin indices.
    m_bin.m_index = index;
    m_bin.m_ipix  = index % npix();
    m_bin.m_ieng  = index / npix();

    // Set pointers
    m_bin.m_cube       = this;
    m_bin.m_counts     = const_cast<double*>(&(m_map.pixels()[index]));
    m_bin.m_energy     = &(m_energies[m_bin.m_ieng]);
    m_bin.m_time       = &m_time;
    m_bin.m_dir        = &(m_dirs[m_bin.m_ipix]);
    m_bin.m_solidangle = &(m_solidangle[m_bin.m_ipix]);
    m_bin.m_ewidth     = &(m_ewidth[m_bin.m_ieng]);
    m_bin.m_ontime     = &m_ontime;

    // Return
    return;
}
