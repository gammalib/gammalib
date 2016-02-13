/***************************************************************************
 *          GCOMEventCube.cpp - COMPTEL event bin container class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file GCOMEventCube.cpp
 * @brief COMPTEL event bin container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFilename.hpp"
#include "GCOMSupport.hpp"
#include "GCOMEventCube.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAXIS                                   "GCOMEventCube::naxis(int)"
#define G_SET_SCATTER_DIRECTIONS    "GCOMEventCube::set_scatter_directions()"
#define G_SET_ENERGIES                        "GCOMEventCube::set_energies()"
#define G_SET_TIMES                              "GCOMEventCube::set_times()"
#define G_SET_BIN                              "GCOMEventCube::set_bin(int&)"

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
 *
 * Constructs an empty event cube.
 ***************************************************************************/
GCOMEventCube::GCOMEventCube(void) : GEventCube()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load constructor
 *
 * @param[in] filename Event cube FITS filename.
 *
 * Construct an event cube by loading information from a FITS file.
 ***************************************************************************/
GCOMEventCube::GCOMEventCube(const GFilename& filename) : GEventCube()
{
    // Initialise members
    init_members();

    // Load event cube
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] map Sky map.
 * @param[in] ebds Energy boundaries.
 * @param[in] gti Good Time intervals.
 * @param[in] phimin Minimum scatter angle.
 * @param[in] dphi Scatter angle step size.
 *
 * Construct an event cube from sky map, energy boundaries, Good Time
 * Intervals and the scatter angle minimum value and step size.
 ***************************************************************************/
GCOMEventCube::GCOMEventCube(const GSkyMap&  map,
                             const GEbounds& ebds,
                             const GGti&     gti,
                             const double&   phimin,
                             const double&   dphi) : GEventCube()
{
    // Initialise members
    init_members();

    // Set sky map, energy boundaries and GTI
    m_map = map;
    this->ebounds(ebds);
    this->gti(gti);

    // Set scatter directions
    set_scatter_directions();

    // Set scatter angles
    set_scatter_angles(phimin, dphi);

    // Set energies
    set_energies();

    // Set times
    set_times();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Event cube.
 ***************************************************************************/
GCOMEventCube::GCOMEventCube(const GCOMEventCube& cube) : GEventCube(cube)
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
GCOMEventCube::~GCOMEventCube(void)
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
 * @return Event cube.
 ***************************************************************************/
GCOMEventCube& GCOMEventCube::operator=(const GCOMEventCube& cube)
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
 * Returns pointer to an event bin. Note that the returned pointer is in
 * fact always the same, but the method sets the pointers within the
 * event bin so that they point to the appropriate information.
 ***************************************************************************/
GCOMEventBin* GCOMEventCube::operator[](const int& index)
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
 * @return Const pointer to event bin.
 *
 * Returns pointer to an event bin. Note that the returned pointer is in
 * fact always the same, but the method sets the pointers within the
 * event bin so that they point to the appropriate information.
 ***************************************************************************/
const GCOMEventBin* GCOMEventCube::operator[](const int& index) const
{
    // Set event bin (circumvent const correctness)
    ((GCOMEventCube*)this)->set_bin(index);

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
void GCOMEventCube::clear(void)
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
 * @return Pointer to deep copy of event cube.
 ***************************************************************************/
GCOMEventCube* GCOMEventCube::clone(void) const
{
    return new GCOMEventCube(*this);
}


/***********************************************************************//**
 * @brief Return number of bins in event cube
 *
 * @return Number of bins in event cube.
 *
 * The number of bins in the event cube is the product of the number of bins
 * along the three cube axes.
 ***************************************************************************/
int GCOMEventCube::size(void) const
{
    // Compute number of bins
    int nbins = m_map.npix() * m_map.nmaps(); 

    // Return number of bins
    return nbins;
}


/***********************************************************************//**
 * @brief Return dimension of event cube
 *
 * @return Number of dimensions in event cube.
 *
 * The dimension of the cube is either 2 or 3, depending on whether several
 * scatter angle layers exist or not.
 ***************************************************************************/
int GCOMEventCube::dim(void) const
{
    // Compute dimension from sky map
    int dim = (m_map.nmaps() > 1) ? 3 : 2;

    // Return dimension
    return dim;
}


/***********************************************************************//**
 * @brief Return number of bins in axis
 *
 * @param[in] axis Axis [0,...,dim()-1].
 * @return Number of bins in axis.
 *
 * @exception GException::out_of_range
 *            Axis is out of range.
 *
 * Returns the number of bins along a given event cube @p axis.
 ***************************************************************************/
int GCOMEventCube::naxis(const int& axis) const
{
    // Optionally check if the axis is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= dim()) {
        throw GException::out_of_range(G_NAXIS, "COMPTEL event cube axis", axis, dim());
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
        naxis = m_map.nmaps();
        break;
    }

    // Return result
    return naxis;
}


/***********************************************************************//**
 * @brief Load COMPTEL event cube from FITS file
 *
 * @param[in] filename FITS filename.
 *
 * The method clears the object before loading, thus any events residing in
 * the object before loading will be lost.
 ***************************************************************************/
void GCOMEventCube::load(const GFilename& filename)
{
    // Clear object
    clear();

    // Open DRE FITS file
    GFits file(filename.filename());

    // Load DRE cube
    read(file);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save COMPTEL event cube into FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite existing FITS file (default=false).
 *
 * Save the COMPTEL event cube into FITS file.
 ***************************************************************************/
void GCOMEventCube::save(const GFilename& filename, const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write event cube
    write(fits);
    
    // Save FITS file
    fits.saveto(filename.filename(), clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL event cube from FITS file
 *
 * @param[in] fits FITS file.
 *
 * Reads an COMPTEL event cube from a DRE FITS file. The event cube is loaded
 * into a skymap and the header keywords are analyzed to determine the
 * scatter angle axis, the energy range and the time range.
 ***************************************************************************/
void GCOMEventCube::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get image
    const GFitsImage& image = *fits.image("Primary");

    // Load counts map as sky map
    m_map.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    com_wcs_mer2car(m_map);

    // Set scatter directions
    set_scatter_directions();

    // Extract Phibar minimum and step size
    double crval3 = image.real("CRVAL3");
    double crpix3 = image.real("CRPIX3");
    double dphi   = image.real("CDELT3");
    double phimin = crval3 + (crpix3-1.0) * dphi;

    // Set scatter angles
    set_scatter_angles(phimin, dphi);

    // Extract energy range
    GEnergy emin;
    GEnergy emax;
    emin.MeV(image.real("E_MIN"));
    emax.MeV(image.real("E_MAX"));

    // Set energy boundaries
    m_ebounds.clear();
    m_ebounds.append(emin, emax);

    // Set energies
    set_energies();

    // Extract time range
    GTime  tstart;
    GTime  tstop;
    tstart.jd(image.real("TSTART"));
    tstop.jd(image.real("TSTOP"));

    // Set GTIs
    m_gti.clear();
    m_gti.append(tstart, tstop);

    // Set times
    set_times();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write COMPTEL event cube into FITS file.
 *
 * @param[in] file FITS file.
 *
 * @todo Write energy boundaries and GTIs into FITS header
 ***************************************************************************/
void GCOMEventCube::write(GFits& file) const
{
    // Write cube
    m_map.write(file);

    // Write energy boundaries

    // Write Good Time intervals

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of events in cube
 *
 * @return Number of events in event cube.
 *
 * This method returns the number of events in the event cube rounded to the
 * nearest integer.
 ***************************************************************************/
int GCOMEventCube::number(void) const
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
 *
 * @param[in] map Sky map.
 * @param[in] phimin Minimum scatter angle.
 * @param[in] dphi Scatter angle step size.
 *
 * Sets an event cube using a skymap and information about the scatter
 * angle minimum and step size.
 *
 * @todo Check if we really need this method.
 ***************************************************************************/
void GCOMEventCube::map(const GSkyMap& map, const double& phimin,
                        const double& dphi)
{
    // Store sky map
    m_map = map;

    // Compute sky directions
    set_scatter_directions();

    // Set scatter angles
    set_scatter_angles(phimin, dphi);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print event cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing event cube information.
 ***************************************************************************/
std::string GCOMEventCube::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMEventCube ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of events"));
        result.append(gammalib::str(number()));
        result.append("\n"+gammalib::parformat("Number of elements"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Size (Chi x Psi x Phi)"));
        result.append(gammalib::str(nchi())+" x ");
        result.append(gammalib::str(npsi())+" x ");
        result.append(gammalib::str(nphi()));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(gammalib::str(emin().MeV())+" - ");
        result.append(gammalib::str(emax().MeV())+" MeV");
        result.append("\n"+gammalib::parformat("Mean energy"));
        result.append(m_energy.print(chatter));
        result.append("\n"+gammalib::parformat("Energy bin width"));
        result.append(m_ewidth.print(chatter));
        result.append("\n"+gammalib::parformat("Time interval"));
        result.append(gammalib::str(tstart().jd())+" - ");
        result.append(gammalib::str(tstop().jd())+" Julian days");
        result.append("\n"+gammalib::parformat("Mean time"));
        result.append(m_time.print(chatter));
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(m_ontime)+" s");

        // Append skymap definition
        result.append("\n"+m_map.print(chatter));

    } // endif: chatter was not silent

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
 *
 * The method initialises the class members to a well defined state. It also
 * prepares the event bin member that will be returned in case of an operator
 * access to the class. As the event bin member manipulates pointers, it will
 * be possible to directly write to the memory that has been associated with
 * a bin.
 ***************************************************************************/
void GCOMEventCube::init_members(void)
{
    // Initialise members
    m_bin.clear();
    m_dir.clear();
    m_map.clear();
    m_time.clear();
    m_ontime = 0.0;
    m_energy.clear();
    m_ewidth.clear();
    m_dirs.clear();
    m_solidangle.clear();
    m_phi.clear();
    m_dphi.clear();

    // Prepare event bin
    init_bin();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Event cube.
 *
 * This method copies the class members from another event cube in the actual
 * object. It also prepares the event bin member that will be returned in
 * case of an operator access to the class.
 ***************************************************************************/
void GCOMEventCube::copy_members(const GCOMEventCube& cube)
{
    // Copy members
    m_dir        = cube.m_dir;
    m_map        = cube.m_map;
    m_time       = cube.m_time;
    m_ontime     = cube.m_ontime;
    m_energy     = cube.m_energy;
    m_ewidth     = cube.m_ewidth;
    m_dirs       = cube.m_dirs;
    m_solidangle = cube.m_solidangle;
    m_phi        = cube.m_phi;
    m_dphi       = cube.m_dphi;

    // Prepare event bin
    init_bin();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMEventCube::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sky directions and solid angles of events cube
 *
 * @exception GException::invalid_value
 *            No sky pixels have been defined.
 *
 * This method computes the sky directions and solid angles for all (Chi,Psi)
 * values of the event cube. Sky directions are stored in an array of GSkyDir
 * objects while solid angles are stored in units of sr in an array of double
 * precision variables.
 ***************************************************************************/
void GCOMEventCube::set_scatter_directions(void)
{
    // Throw an error if we have no sky pixels
    if (npix() < 1) {
        std::string msg = "No sky pixels have been found in event cube. "
                          "Every COMPTEL event cube needs a definition of "
                          "sky pixels.";
        throw GException::invalid_value(G_SET_SCATTER_DIRECTIONS, msg);
    }

    // Clear vectors
    m_dirs.clear();
    m_solidangle.clear();

    // Reserve space for pixel directions and solid angles
    m_dirs.reserve(npix());
    m_solidangle.reserve(npix());

    // Set pixel directions and solid angles
    for (int iy = 0; iy < npsi(); ++iy) {
        for (int ix = 0; ix < nchi(); ++ix) {
            GSkyPixel pixel = GSkyPixel(double(ix), double(iy));
            m_dirs.push_back(m_map.pix2dir(pixel));
            m_solidangle.push_back(m_map.solidangle(pixel));
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log mean energy and energy width of event cube
 *
 * @param[in] phimin Scatter angle of first layer.
 * @param[in] dphi Scatter angle difference between layers.
 *
 * This method sets a vector of scatter directions occuring in the event
 * cube. It assumes that scatter directions are linearly space.
 ***************************************************************************/
void GCOMEventCube::set_scatter_angles(const double& phimin, const double& dphi)
{
    // Clear vectors
    m_phi.clear();
    m_dphi.clear();

    // Reserve space for pixel directions and solid angles
    m_phi.reserve(nphi());
    m_dphi.reserve(nphi());

    // Set scatter angles
    for (int iz = 0; iz < nphi(); ++iz) {
        double phi = phimin + iz*dphi;
        m_phi.push_back(phi);
        m_dphi.push_back(dphi);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log mean energy and energy width of event cube
 *
 * @exception GException::invalid_value
 *            No energy boundaries found.
 ***************************************************************************/
void GCOMEventCube::set_energies(void)
{
    // Throw an error if GTI is empty
    if (m_ebounds.size() < 1) {
        std::string msg = "No energy boundaries have been found in event "
                          "cube. Every COMPTEL event cube needs a definition "
                          "of the energy boundaries.";
        throw GException::invalid_value(G_SET_ENERGIES, msg);
    }

    // Compute the logarithmic mean energy
    m_energy = m_ebounds.elogmean(0);

    // Compute the energy bin size
    m_ewidth = m_ebounds.emax(0) - m_ebounds.emin(0);

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
 ***************************************************************************/
void GCOMEventCube::set_times(void)
{
    // Throw an error if GTI is empty
    if (m_gti.size() < 1) {
        std::string msg = "No Good Time Intervals have been found in event "
                          "cube. Every COMPTEL event cube needs a definition "
                          "of the Good Time Intervals.";
        throw GException::invalid_value(G_SET_TIMES, msg);
    }
    
    // Compute mean time
    m_time = 0.5 * (m_gti.tstart() + m_gti.tstop());

    // Set ontime
    m_ontime = m_gti.ontime();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise event bin
 *
 * This method initialises the event bin. The event bin is cleared and all
 * fixed pointers are set. Only the m_counts and the m_solidangle member of the
 * event bin will be set to NULL, but these will be set by the set_bin method
 * which is called before any event bin access.
 ***************************************************************************/
void GCOMEventCube::init_bin(void)
{
    // Prepare event bin
    m_bin.free_members();
    m_bin.m_counts     = NULL;      //!< Will be set by set_bin method
    m_bin.m_dir        = &m_dir;    //!< Content will be set by set_bin method
    m_bin.m_solidangle = NULL;      //!< Will be set by set_bin method
    m_bin.m_time       = &m_time;   //!< Fixed content
    m_bin.m_ontime     = &m_ontime; //!< Fixed content
    m_bin.m_energy     = &m_energy; //!< Fixed content
    m_bin.m_ewidth     = &m_ewidth; //!< Fixed content

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
 *            Sky directions and solid angles vectors have not been set up.
 *
 * This method provides the event attributes to the event bin. The event bin
 * is in fact physically stored in the event cube, and only a single event
 * bin is indeed allocated. This method sets up the pointers in the event
 * bin so that a client can easily access the information of individual bins
 * as if they were stored in an array.
 ***************************************************************************/
void GCOMEventCube::set_bin(const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET_BIN, index, 0, size()-1);
    }
    #endif

    // Check for the existence of sky directions and solid angles
    if (m_dirs.size() != npix() || m_solidangle.size() != npix()) {
        std::string msg = "Sky direction vector has not been setup for "
                          "event cube. Every COMPTEL event cube needs an "
                          "internal vector of sky directions.";
        throw GException::invalid_value(G_SET_BIN, msg);
    }

    // Get pixel and energy bin indices.
    int ipix = index % npix();
    int iphi = index / npix();

    // Set indices
    m_bin.m_index = index;

    // Set instrument direction
    m_dir.dir(m_dirs[ipix]);
    m_dir.phibar(m_phi[iphi]);
    
    // Set pointers
    m_bin.m_counts      = const_cast<double*>(&(m_map.pixels()[index]));
    m_bin.m_solidangle  = &(m_solidangle[ipix]);

    // Return
    return;
}
