/***************************************************************************
 *           GCTAEventCube.cpp  -  CTA event bin container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAEventCube.cpp
 * @brief CTA event bin container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GFits.hpp"
#include "GCTAEventCube.hpp"
#include "GCTATypemaps.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAXIS                                   "GCTAEventCube::naxis(int)"
#define G_ENERGY                                "GCTAEventCube::energy(int&)"
#define G_SET_DIRECTIONS                    "GCTAEventCube::set_directions()"
#define G_SET_ENERGIES                        "GCTAEventCube::set_energies()"
#define G_SET_TIMES                              "GCTAEventCube::set_times()"
#define G_SET_BIN                              "GCTAEventCube::set_bin(int&)"

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
 * @brief File name constructor
 *
 * @param[in] filename Counts cube filename.
 *
 * Constructs instance of events cube from a counts cube file.
 ***************************************************************************/
GCTAEventCube::GCTAEventCube(const GFilename& filename) : GEventCube()
{
    // Initialise members
    init_members();

    // Load counts cube
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
 *
 * Constructs instance of events cube from a sky map, energy boundaries and
 * Good Time Intervals. All event cube weights are set to unity.
 ***************************************************************************/
GCTAEventCube::GCTAEventCube(const GSkyMap&  map,
                             const GEbounds& ebds,
                             const GGti&     gti) : GEventCube()
{
    // Initialise members
    init_members();

    // Set sky map, energy boundaries and GTI
    m_map = map;
    this->ebounds(ebds);
    this->gti(gti);

    // Set sky directions
    set_directions();

    // Set energies
    set_energies();

    // Set times
    set_times();

    // Set all weights to unity
    m_weights = m_map;
    m_weights = 1.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] map Sky map.
 * @param[in] weights Event cube weights.
 * @param[in] ebds Energy boundaries.
 * @param[in] gti Good Time intervals.
 *
 * Constructs instance of events cube from a sky map, a map of weights,
 * energy boundaries and Good Time Intervals.
 ***************************************************************************/
GCTAEventCube::GCTAEventCube(const GSkyMap&  map,
                             const GSkyMap&  weights,
                             const GEbounds& ebds,
                             const GGti&     gti) : GEventCube()
{
    // Initialise members
    init_members();

    // Set sky map, energy boundaries and GTI
    m_map = map;
    this->ebounds(ebds);
    this->gti(gti);

    // Set weight map
    m_weights = weights;

    // Set sky directions
    set_directions();

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
 * @return Event cube
 ***************************************************************************/
GCTAEventCube& GCTAEventCube::operator=(const GCTAEventCube& cube)
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
 *
 * Returns pointer to an event bin.
 ***************************************************************************/
GCTAEventBin* GCTAEventCube::operator[](const int& index)
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
 *
 * Returns pointer to an event bin.
 ***************************************************************************/
const GCTAEventBin* GCTAEventCube::operator[](const int& index) const
{
    // Set event bin (circumvent const correctness)
    (const_cast<GCTAEventCube*>(this))->set_bin(index);

    // Return pointer
    return (&m_bin);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear CTA event cube
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
 * @brief Clone CTA event cube
 *
 * @return Pointer to deep copy of CTA event cube.
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
int GCTAEventCube::naxis(const int& axis) const
{
    // Optionally check if the axis is valid
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= dim()) {
        throw GException::out_of_range(G_NAXIS, "CTA event cube axis", axis, dim());
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
 * @brief Load CTA event cube from FITS file
 *
 * @param[in] filename FITS filename.
 *
 * Loads the event cube from a FITS file. See the read() method for more
 * information about the structure of the FITS file.
 *
 * The method clears the object before loading, thus any events residing in
 * the object before loading will be lost.
 ***************************************************************************/
void GCTAEventCube::load(const GFilename& filename)
{
    // Clear object
    clear();

    // Open counts map FITS file
    GFits fits(filename);

    // Load counts map
    read(fits);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save CTA event cube into FITS file
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite existing FITS file.
 *
 * Save the CTA event cube into FITS file. See the write() method for more
 * information about the structure of the FITS file.
 ***************************************************************************/
void GCTAEventCube::save(const GFilename& filename,
                         const bool&      clobber) const
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
 * @brief Read CTA event cube from FITS file
 *
 * @param[in] fits FITS file.
 *
 * Read an event cube from a FITS file. The following HDUs will be read
 *
 *      COUNTS - Counts cube (or primary extension if COUNTS does not exist)
 *      WEIGHTS - Weights for each counts cube bin (optional)
 *      EBOUNDS - Energy boundaries
 *      GTI - Good Time Intervals
 *
 * The method clears the event cube before reading, thus any events residing
 * in the event cube will be lost.
 ***************************************************************************/
void GCTAEventCube::read(const GFits& fits)
{
    // Clear object
    clear();

    // Set counts cube HDU
    std::string counts_hdu("Primary");
    if (fits.contains(gammalib::extname_cta_counts)) {
        counts_hdu = gammalib::extname_cta_counts;
    }
    
    // Get HDUs
    const GFitsImage& hdu_cntmap  = *fits.image(counts_hdu);
    const GFitsTable& hdu_ebounds = *fits.table(gammalib::extname_ebounds);
    const GFitsTable& hdu_gti     = *fits.table(gammalib::extname_gti);

    // Load counts map
    read_cntmap(hdu_cntmap);

    // Load energy boundaries
    read_ebds(hdu_ebounds);

    // Load GTIs
    read_gti(hdu_gti);

    // If a WEIGHTS HDU exist then load it from the FITS file ...
    if (fits.contains(gammalib::extname_cta_weights)) {

        // Get WEIGHTS HDU
        const GFitsImage& hdu_weights = *fits.image(gammalib::extname_cta_weights);

        // Read HDU
        m_weights.read(hdu_weights);

    }

    // ... otherwise set the weight map to unity
    else {
        m_weights = m_map;
        m_weights = 1.0;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA event cube into FITS file.
 *
 * @param[in] fits FITS file.
 *
 * Writes CTA event cube into a FITS file. The following HDUs will be written
 *
 *      COUNTS - Counts cube
 *      WEIGHTS - Weights for each counts cube bin
 *      EBOUNDS - Energy boundaries
 *      GTI - Good Time Intervals
 *
 * The counts cube will be written as a double precision sky map. The
 * weighing cube contains the weight for each counts cube, computed as the
 * fraction of the event bin that has been selected in filling the counts
 * cube. This allows to account for proper stacking of observations with
 * different energy thresholds, or different regions of interest. The energy
 * boundaries for all counts cube layers are also written, as well as the
 * Good Time Intervals that have been used in generating the counts cube.
 ***************************************************************************/
void GCTAEventCube::write(GFits& fits) const
{
    // Remove HDUs if they exist already
    if (fits.contains(gammalib::extname_gti)) {
        fits.remove(gammalib::extname_gti);
    }
    if (fits.contains(gammalib::extname_ebounds)) {
        fits.remove(gammalib::extname_ebounds);
    }
    if (fits.contains(gammalib::extname_cta_weights)) {
        fits.remove(gammalib::extname_cta_weights);
    }
    if (fits.contains(gammalib::extname_cta_counts)) {
        fits.remove(gammalib::extname_cta_counts);
    }
    if (fits.contains("Primary")) {
        fits.remove("Primary");
    }

    // Write counts cube
    m_map.write(fits, gammalib::extname_cta_counts);

    // Write cube weighting
    m_weights.write(fits, gammalib::extname_cta_weights);

    // Write energy boundaries
    ebounds().write(fits);

    // Write Good Time intervals
    gti().write(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return number of events in cube
 *
 * @return Number of events in cube, rounded to nearest integer value.
 *
 * Returns the total number of events in the cube, rounded to the nearest
 * integer value. All cube bins with a negative content will be excluded
 * from the sum.
 ***************************************************************************/
int GCTAEventCube::number(void) const
{
    // Initialise result
    double number = 0.0;

    // Get pointer on skymap pixels
    const double* pixels = m_map.pixels();

    // Sum event cube
    if (size() > 0 && pixels != NULL) {
        for (int i = 0; i < size(); ++i) {
            if (pixels[i] > 0.0) {
                number += pixels[i];
            }
        }
    }

    // Return
    return int(number+0.5);
}


/***********************************************************************//**
 * @brief Set event cube counts from sky map
 *
 * @param[in] counts Event cube counts sky map.
 *
 * Sets event cube counts from sky map. The methods also sets all weights to
 * unity.
 ***************************************************************************/
void GCTAEventCube::counts(const GSkyMap& counts)
{
    // Store sky map
    m_map = counts;

    // Compute sky directions
    set_directions();

    // If event cube has pointing then set DETX and DETY coordinates of
    // instrument direction
    if (m_has_pnt) {
        set_detxy(m_pnt);
    }

    // Set all weights to unity
    m_weights = m_map;
    m_weights = 1.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return energy of cube layer
 *
 * @param[in] index Event cube layer index [0,...,ebins()-1]
 * @return Energy of event cube layer
 *
 * @exception GException::out_of_range
 *            Layer index is out of range.
 *
 * Returns the energy of the event cube layer specified by @p index.
 ***************************************************************************/
const GEnergy& GCTAEventCube::energy(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= ebins()) {
        throw GException::out_of_range(G_ENERGY,
              "CTA event cube energy axis", index, ebins());
    }
    #endif

    // Return energy
    return (m_energies[index]);
}


/***********************************************************************//**
 * @brief Print event cube information
 *
 * @param[in] chatter Chattiness.
 * @return String containing event cube information.
 ***************************************************************************/
std::string GCTAEventCube::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAEventCube ===");
        result.append("\n"+gammalib::parformat("Number of events") +
                      gammalib::str(number()));
        result.append("\n"+gammalib::parformat("Number of elements") +
                      gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Number of pixels") +
                      gammalib::str(npix()));
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(ebins()));

        // Append GTI interval
        result.append("\n"+gammalib::parformat("Time interval"));
        if (gti().size() > 0) {
            result.append(gammalib::str(tstart().mjd()));
            result.append(" - ");
            result.append(gammalib::str(tstop().mjd())+" days");
        }
        else {
            result.append("not defined");
        }

        // Append pointing
        result.append("\n"+gammalib::parformat("Pointing"));
        if (m_has_pnt) {
            result.append(m_pnt.dir().print());
        }
        else {
            result.append("not defined");
        }
    
        // Append energy intervals
        if (gammalib::reduce(chatter) > SILENT) {
            if (ebounds().size() > 0) {
                result.append("\n"+ebounds().print(gammalib::reduce(chatter)));
            }
            else {
                result.append("\n"+gammalib::parformat("Energy intervals") +
                              "not defined");
            }
        }
        else {
            result.append("\n"+gammalib::parformat("Energy interval"));
            if (ebounds().size() > 0) {
                result.append(gammalib::str(emin().TeV()));
                result.append(" - ");
                result.append(gammalib::str(emax().TeV())+" TeV");
            }
            else {
                result.append("not defined");
            }
        }

        // Append skymap definition
        if (gammalib::reduce(chatter) > SILENT) {
            result.append("\n"+m_map.print(gammalib::reduce(chatter)));    
        }

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
 ***************************************************************************/
void GCTAEventCube::init_members(void)
{
    // Initialise members
    m_map.clear();
    m_weights.clear();
    m_bin.clear();
    m_time.clear();
    m_polarization.clear();
    m_pnt.clear();
    m_has_pnt = false;
    m_dirs.clear();
    m_solidangle.clear();
    m_energies.clear();
    m_ewidth.clear();
    m_ontime = 0.0;

    // Prepare event bin
    init_bin();

    // Set CTA time reference for GTIs
    m_gti.reference(GTimeReference(G_CTA_MJDREF, "s", "TT", "LOCAL"));

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
void GCTAEventCube::copy_members(const GCTAEventCube& cube)
{
    // Copy members. Note that the event bin is not copied as it will
    // be initialised later. The event bin serves just as a container of
    // pointers, hence we do not want to copy over the pointers from the
    // original class.
    m_map          = cube.m_map;
    m_weights      = cube.m_weights;
    m_time         = cube.m_time;
    m_polarization = cube.m_polarization;
    m_pnt          = cube.m_pnt;
    m_has_pnt      = cube.m_has_pnt;
    m_dirs         = cube.m_dirs;
    m_solidangle   = cube.m_solidangle;
    m_energies     = cube.m_energies;
    m_ewidth       = cube.m_ewidth;
    m_ontime       = cube.m_ontime;

    // Copy GTIs
    m_gti = cube.m_gti;

    // Prepare event bin
    init_bin();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventCube::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA counts map from HDU.
 *
 * @param[in] hdu Image HDU.
 *
 * This method reads a CTA counts map from a FITS HDU. The counts map is
 * stored in a GSkyMap object.
 ***************************************************************************/
void GCTAEventCube::read_cntmap(const GFitsImage& hdu)
{
    // Load counts map as sky map
    m_map.read(hdu);

    // Set sky directions
    set_directions();

    // If header contains pointing direction then also set DETX and DETY
    // coordinates
    if (hdu.has_card("RA_PNT") && hdu.has_card("DEC_PNT")) {

        // Read pointing direction
        double ra_pnt  = hdu.real("RA_PNT");
        double dec_pnt = hdu.real("DEC_PNT");
        GSkyDir pnt;
        pnt.radec_deg(ra_pnt, dec_pnt);
        m_pnt.dir(pnt);
        m_has_pnt = true;

        // Set DETX and DETY coordinates of instrument direction
        set_detxy(m_pnt);

    } // endif: header contained pointing direction

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy boundaries from HDU.
 *
 * @param[in] hdu Energy boundaries table.
 *
 * Read the energy boundaries from the HDU.
 ***************************************************************************/
void GCTAEventCube::read_ebds(const GFitsTable& hdu)
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
 * Reads the Good Time Intervals from the HDU.
 ***************************************************************************/
void GCTAEventCube::read_gti(const GFitsTable& hdu)
{
    // Read Good Time Intervals
    m_gti.read(hdu);

    // Set times
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
 * cube pixels. Sky directions are stored in an array of GCTAInstDir objects
 * while solid angles are stored in units of sr in an array of double
 * precision variables.
 *
 * A kluge has been introduced that handles invalid pixels. Invalid pixels
 * may occur if a Hammer-Aitoff projection is used. In this case, pixels may
 * lie outside the valid sky region. As invalid pixels lead to exceptions
 * in the WCS classes, we simply need to catch the exceptions here. Invalid
 * pixels are signaled by setting the solid angle of the pixel to 0.
 ***************************************************************************/
void GCTAEventCube::set_directions(void)
{
    // Throw an error if we have no sky pixels
    if (npix() < 1) {
        std::string msg = "CTA event cube contains no sky pixels. Please "
                          "provide a valid CTA event cube.";
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
            try {
                GSkyPixel pixel = GSkyPixel(double(ix), double(iy));
                m_dirs.push_back(GCTAInstDir(m_map.pix2dir(pixel)));
                m_solidangle.push_back(m_map.solidangle(pixel));
            }
            catch (GException::invalid_argument& e) {
                m_dirs.push_back(GCTAInstDir());
                m_solidangle.push_back(0.0);
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set DETX and DETY coordinates
 *
 * @param[in] pnt CTA pointing.
 *
 * Computes DETX and DETY coordinates for a given pointing direction.
 ***************************************************************************/
void GCTAEventCube::set_detxy(const GCTAPointing& pnt)
{
    // Get number of instrument directions
    int size = m_dirs.size();

    // Loop over all intstrument directions
    for (int i = 0; i < size; ++i) {

        // Get sky direction
        GSkyDir skydir = m_dirs[i].dir();

        // Set instrument direction from sky direction
        m_dirs[i] = pnt.instdir(skydir);

    } // endfor: looped over all instrument directions

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
void GCTAEventCube::set_energies(void)
{
    // Determine number of energy bins
    int ebins = ebounds().size();

    // Throw an error if we have no energy bins
    if (ebins < 1) {
        std::string msg = "CTA event cube contains no energy boundaries. "
                          "Please provide a valid CTA event cube.";
        throw GException::invalid_value(G_SET_ENERGIES, msg);
    }

    // Clear old bin energies and energy widths
    m_energies.clear();
    m_ewidth.clear();

    // Reserve space for bin energies and energy widths
    m_energies.reserve(ebins);
    m_ewidth.reserve(ebins);

    // Setup bin energies and energy widths
    for (int i = 0; i < ebins; ++i) {
        m_energies.push_back(ebounds().elogmean(i));
        m_ewidth.push_back(ebounds().emax(i) -  ebounds().emin(i));
    }

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
void GCTAEventCube::set_times(void)
{
    // Throw an error if GTI is empty
    if (m_gti.size() < 1) {
        std::string msg = "No Good Time Intervals have been found in event "
                          "cube. Every CTA event cube needs a definition "
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
 * @brief Initialise event bin
 *
 * This method initialises the event bin. The event bin is cleared and all
 * fixed pointers are set.
 ***************************************************************************/
void GCTAEventCube::init_bin(void)
{
    // Free any existing memory
    m_bin.free_members();

    // Set fixed pointers (those will not be set in set_bin)
    m_bin.m_time         = &m_time;
    m_bin.m_ontime       = &m_ontime;
    m_bin.m_polarization = &m_polarization;

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
void GCTAEventCube::set_bin(const int& index)
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

    // Set pixel and energy bin indices.
    m_bin.m_ipix = index % npix();
    m_bin.m_ieng = index / npix();

    // Set pointers
    m_bin.m_counts     = const_cast<double*>(&(m_map.pixels()[index]));
    m_bin.m_energy     = &(m_energies[m_bin.m_ieng]);
    m_bin.m_dir        = &(m_dirs[m_bin.m_ipix]);
    m_bin.m_solidangle = &(m_solidangle[m_bin.m_ipix]);
    m_bin.m_ewidth     = &(m_ewidth[m_bin.m_ieng]);
    m_bin.m_weight     = const_cast<double*>(&(m_weights.pixels()[index]));

    // Return
    return;
}
