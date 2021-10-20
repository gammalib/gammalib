/***************************************************************************
 *             GLATEventCube.cpp - Fermi/LAT event cube class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2021 by Juergen Knoedlseder                         *
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
 * @brief Fermi/LAT event cube class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GFilename.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"
#include "GLATEventCube.hpp"

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
#define G_MISSING_WCS_KLUDGE         //!< Handle source maps with missing WCS

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
 * @brief File name constructor
 *
 * @param[in] filename Counts cube filename.
 *
 * Construct event cube object by loading the events from a FITS file.
 ***************************************************************************/
GLATEventCube::GLATEventCube(const GFilename& filename) : GEventCube()
{
    // Initialise members
    init_members();

    // Load counts cube
    load(filename);

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
GLATEventCube& GLATEventCube::operator=(const GLATEventCube& cube)
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
        throw GException::out_of_range(G_NAXIS, "LAT event cube axis",
                                       axis, dim());
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
void GLATEventCube::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read event cube from FITS file
    read(fits);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save LAT event cube into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing FITS file? (default: false)
 *
 * Save the LAT event cube into FITS file.
 ***************************************************************************/
void GLATEventCube::save(const GFilename& filename,
                         const bool&      clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write event cube into FITS file
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
 * FITS file, the energy boundaries reside in the `EBOUNDS` extension and the
 * Good Time Intervals reside in the `GTI` extension.  The method clears the
 * object before loading, thus any events residing in the object before
 * loading will be lost.
 ***************************************************************************/
void GLATEventCube::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get HDUs
    const GFitsImage& hdu_cntmap  = *fits.image("Primary");
    const GFitsTable& hdu_ebounds = *fits.table(gammalib::extname_ebounds);
    const GFitsTable& hdu_gti     = *fits.table(gammalib::extname_gti);

    // Load counts map
    read_cntmap(hdu_cntmap);

    // Load energy boundaries
    read_ebds(hdu_ebounds);

    // Load GTIs
    read_gti(hdu_gti);

    // Load additional source maps
    for (int i = 1; i < fits.size(); ++i) {

        // Only consider image extensions
        if (fits.at(i)->exttype() == GFitsHDU::HT_IMAGE) {

            // Get FITS image
            const GFitsImage& hdu_srcmap = *fits.image(i);

            // Handle files that do not have WCS in their additional
            // source map extensions
            #if defined(G_MISSING_WCS_KLUDGE)
            std::string keys[] = {"CTYPE1", "CTYPE2", "CTYPE3",
                                  "CRPIX1", "CRPIX2", "CRPIX3",
                                  "CRVAL1", "CRVAL2", "CRVAL3",
                                  "CDELT1", "CDELT2", "CDELT3",
                                  "CUNIT1", "CUNIT2", "CUNIT3",
                                  "CROTA2"};
            for (int i = 0; i < 16; ++i) {
                if (!hdu_srcmap.has_card(keys[i])) {
                    const_cast<GFitsImage&>(hdu_srcmap).card(hdu_cntmap.card(keys[i]));
                }
            }
            #endif

            // Read source map
            read_srcmap(hdu_srcmap);

        } // endif: extension is an image

    } // endfor: looped over additional source maps

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
void GLATEventCube::map(const GSkyMap& map)
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
 * @param[in] chatter Chattiness.
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
        throw GException::out_of_range(G_DIFFNAME, "Diffuse model index",
                                       index, ndiffrsp());
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
GSkyMap* GLATEventCube::diffrsp(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= ndiffrsp()) {
        throw GException::out_of_range(G_DIFFRSP, "Diffuse model index",
                                       index, ndiffrsp());
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
    double radius = 0.0;

    // Continue only if sky direction is within sky map
    if (m_map.contains(srcDir)) {

        // Set to largest possible radius
        radius = 180.0;

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
    
    } // endif: sky direction within sky map

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
    m_polarization.clear();
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
    m_polarization = cube.m_polarization;
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
 * is stored in a GSkyMap object, and a pointer is set up to access the
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
 * @exception GException::invalid_argument
 *            Source map in @p hdu is not compatible with event cube.
 *
 * This method reads a LAT source map from a FITS image. The source map is
 * stored in a GSkyMap object and is given in units of counts/pixel/MeV.
 ***************************************************************************/
void GLATEventCube::read_srcmap(const GFitsImage& hdu)
{
    // Allocate skymap
    GSkyMap* map = new GSkyMap;

    // Read skymap
    map->read(hdu);

    // Check that source map sky projection is consistent with counts
    // map sky projection
    if (*(m_map.projection()) != *(map->projection())) {
        std::string msg = "Specified source map projection in HDU \""+
                          hdu.extname()+"\" is incompatible with map projection "
                          "of the event cube. Please specify a compatible "
                          "map projection.";
        throw GException::invalid_argument(G_READ_SRCMAP, msg);
    }

    // Check that source map dimension is consistent with counts map
    // dimension
    if (m_map.nx() != map->nx() ||
        m_map.ny() != map->ny()) {
        std::string msg = "Number of pixels ("+gammalib::str(map->nx())+","+
                          gammalib::str(map->ny())+") in source map in HDU \""+
                          hdu.extname()+"\" is not compatible with number of "
                          "pixels ("+gammalib::str(m_map.nx())+","+
                          gammalib::str(m_map.ny())+") in the event cube. "
                          "Please specify a source map of compatible size.";
        throw GException::invalid_argument(G_READ_SRCMAP, msg);
    }

    // Check that source map has required number of energy bins
    if (m_map.nmaps()+1 != map->nmaps()) {
        std::string msg = "Number of maps ("+gammalib::str(map->nmaps())+")"+
                          " in source map in HDU \""+hdu.extname()+"\" is not "
                          "compatible with number of energy bin boundaries ("+
                          gammalib::str(m_map.nmaps()+1)+") of the event cube. "
                          "Please specify a source map of compatible size.";
        throw GException::invalid_argument(G_READ_SRCMAP, msg);
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
 * Reads the Good Time Intervals from the GTI extension. Since the Fermi
 * LAT Science Tools do not set corrently the time reference for source
 * maps, the method automatically adds this missing information so that
 * the time reference is set correctly. The time reference that is assumed
 * for Fermi LAT is
 *
 *      MJDREFI 51910
 *      MJDREFF 0.00074287037037037
 * 
 ***************************************************************************/
void GLATEventCube::read_gti(const GFitsTable& hdu)
{
    // Work on a local copy of the HDU to make the kluge work
    GFitsTable* hdu_local = hdu.clone();

    // Kluge: modify HDU table in case that the MJDREF header keyword is
    // blank. This happens for Fermi LAT source maps since the Science
    // Tools do not properly write out the time reference. We hard-code
    // here the Fermi LAT time reference to circumvent the problem.
    if (hdu.has_card("MJDREF")) {
        if (gammalib::strip_whitespace(hdu.string("MJDREF")).empty()) {
            const_cast<GFitsHeader&>(hdu_local->header()).remove("MJDREF");
            hdu_local->card("MJDREFI", 51910,
                            "Integer part of MJD reference");
            hdu_local->card("MJDREFF", 0.00074287037037037,
                            "Fractional part of MJD reference");
        }
    }

    // Read Good Time Intervals
    m_gti.read(*hdu_local);

    // Set time
    set_times();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sky directions and solid angles of events cube.
 *
 * @exception GException::invalid_value
 *            No sky pixels found in event cube.
 *
 * This method computes the sky directions and solid angles for all event
 * cube pixels. Sky directions are stored in an array of GLATInstDir objects
 * while solid angles are stored in units of sr in a double precision array.
 ***************************************************************************/
void GLATEventCube::set_directions(void)
{
    // Throw an exception if we have no sky pixels
    if (npix() < 1) {
        std::string msg = "LAT event cube contains no sky pixels. Please "
                          "provide a valid LAT event cube.";
        throw GException::invalid_value(G_SET_DIRECTIONS, msg);
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
 * @exception GException::invalid_value
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
        std::string msg = "LAT event cube contains no energy bins. Please "
                          "provide a valid LAT event cube.";
        throw GException::invalid_value(G_SET_ENERGIES, msg);
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
 * @brief Set mean event time and ontime of event cube
 *
 * @exception GException::invalid_value
 *            No Good Time Intervals found.
 *
 * Computes the mean time of the event cube by taking the mean between start
 * and stop time. Computes also the ontime by summing up of all good time
 * intervals.
 *
 * @todo Could add a more sophisticated mean event time computation that
 *       weights by the length of the GTIs, yet so far we do not really use
 *       the mean event time, hence there is no rush to implement this.
 ***************************************************************************/
void GLATEventCube::set_times(void)
{
    // Throw an error if GTI is empty
    if (m_gti.size() < 1) {
        std::string msg = "No Good Time Intervals have been found in event "
                          "cube. Every LAT event cube needs a definition "
                          "of the Good Time Intervals.";
        throw GException::invalid_value(G_SET_TIMES, msg);
    }

    // Compute mean time
    m_time = m_gti.tstart() + 0.5 * (m_gti.tstop() - m_gti.tstart());

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
 * @exception GException::invalid_value
 *            Energy vectors have not been set up.
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
        throw GException::out_of_range(G_SET_BIN, "Event index", index, size());
    }
    #endif

    // Check for the existence of energies and energy widths
    if (m_energies.size() != ebins() || m_ewidth.size() != ebins()) {
        std::string msg = "Number of energy bins ("+gammalib::str(ebins())+
                          ") is incompatible with size of energies ("+
                          gammalib::str(m_energies.size())+") an/or size of "
                          "energy widths ("+gammalib::str(m_ewidth.size())+
                          "). Please provide an correctly initialised event "
                          "cube.";
        throw GException::invalid_value(G_SET_BIN, msg);
    }

    // Check for the existence of sky directions and solid angles
    if (m_dirs.size() != npix() || m_solidangle.size() != npix()) {
        std::string msg = "Number of sky pixels ("+gammalib::str(npix())+
                          ") is incompatible with size of sky directions ("+
                          gammalib::str(m_dirs.size())+") and/or size of "
                          "solid angles ("+gammalib::str(m_solidangle.size())+
                          "). Please provide an correctly initialised event "
                          "cube.";
        throw GException::invalid_value(G_SET_BIN, msg);
    }

    // Get pixel and energy bin indices.
    m_bin.m_index = index;
    m_bin.m_ipix  = index % npix();
    m_bin.m_ieng  = index / npix();

    // Set pointers
    m_bin.m_cube         = this;
    m_bin.m_counts       = const_cast<double*>(&(m_map.pixels()[index]));
    m_bin.m_energy       = &(m_energies[m_bin.m_ieng]);
    m_bin.m_time         = &m_time;
    m_bin.m_polarization = &m_polarization;
    m_bin.m_dir          = &(m_dirs[m_bin.m_ipix]);
    m_bin.m_solidangle   = &(m_solidangle[m_bin.m_ipix]);
    m_bin.m_ewidth       = &(m_ewidth[m_bin.m_ieng]);
    m_bin.m_ontime       = &m_ontime;

    // Return
    return;
}
