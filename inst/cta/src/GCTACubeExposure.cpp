/***************************************************************************
 *          GCTACubeExposure.cpp - CTA cube analysis exposure class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2018 by Chia-Chun Lu                                *
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
 * @file GCTACubeExposure.cpp
 * @brief CTA cube analysis exposure class implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GLog.hpp"
#include "GObservations.hpp"
#include "GSkyRegionCircle.hpp"
#include "GCTALib.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_FILL_CUBE    "GCTACubeExposure::fill_cube(GCTAObservation&, GLog*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_CLIP_ETRUES         //!< Clip true energies to reconstructed ones

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTACubeExposure::GCTACubeExposure(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Exposure cube.
 ***************************************************************************/
GCTACubeExposure::GCTACubeExposure(const GCTACubeExposure& cube)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Exposure cube filename.
 *
 * Construct exposure cube by loading the information from an exposure cube
 * file.
 ***************************************************************************/
GCTACubeExposure::GCTACubeExposure(const GFilename& filename)
{
    // Initialise class members
    init_members();

    // Load exposure cube from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Event cube constructor
 *
 * @param[in] cube Event cube.
 *
 * Construct exposure cube using the same binning and sky projection that is
 * used for the event cube.
 ***************************************************************************/
GCTACubeExposure::GCTACubeExposure(const GCTAEventCube& cube)
{
    // Initialise class members
    init_members();

    // Store energy boundaries
    m_energies.set(cube.ebounds());

    // Set GNodeArray used for interpolation
    set_eng_axis();

    // Set exposure cube to event cube
    m_cube = cube.counts();
    m_cube.nmaps(m_energies.size());

    // Set all exposure cube pixels to zero as we want to have a clean map
    // upon construction
    m_cube = 0.0;

    // Return
    return;

}


/***********************************************************************//**
 * @brief Exposure cube constructor
 *
 * @param[in] wcs      World Coordinate System.
 * @param[in] coords   Coordinate System (CEL or GAL).
 * @param[in] x        X coordinate of sky map centre (deg).
 * @param[in] y        Y coordinate of sky map centre (deg).
 * @param[in] dx       Pixel size in x direction at centre (deg/pixel).
 * @param[in] dy       Pixel size in y direction at centre (deg/pixel).
 * @param[in] nx       Number of pixels in x direction.
 * @param[in] ny       Number of pixels in y direction.
 * @param[in] energies Energies.
 *
 * Constructs an exposure cube by specifying the sky map grid and the
 * energies.
 ***************************************************************************/
GCTACubeExposure::GCTACubeExposure(const std::string&   wcs,
                                   const std::string&   coords,
                                   const double&        x,
                                   const double&        y,
                                   const double&        dx,
                                   const double&        dy,
                                   const int&           nx,
                                   const int&           ny,
                                   const GEnergies&     energies)
{
    // Initialise class members
    init_members();

    // Store energies
    m_energies = energies;

    // Set GNodeArray used for interpolation
    set_eng_axis();

    // Create sky map
    m_cube = GSkyMap(wcs, coords, x, y, dx, dy, nx, ny, m_energies.size());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTACubeExposure::~GCTACubeExposure(void)
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
 * @param[in] cube Exposure cube.
 * @return Exposure cube.
 ***************************************************************************/
GCTACubeExposure& GCTACubeExposure::operator=(const GCTACubeExposure& cube)
{
    // Execute only if object is not identical
    if (this != &cube) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(cube);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return exposure (in units of cm2 s)
 *
 * @param[in] dir Coordinate of the true photon position.
 * @param[in] energy Energy of the true photon.
 * @return Exposure (in units of cm2 s)
 ***************************************************************************/
double GCTACubeExposure::operator()(const GSkyDir& dir, const GEnergy& energy) const
{
    // Set indices and weighting factors for interpolation
    update(energy.log10TeV());

    // Perform interpolation
    double exposure = m_wgt_left  * m_cube(dir, m_inx_left) +
                      m_wgt_right * m_cube(dir, m_inx_right);

    // Make sure that exposure does not become negative
    if (exposure < 0.0) {
        exposure = 0.0;
    }

    // Return exposure
    return exposure;
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
void GCTACubeExposure::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone exposure cube
 *
 * @return Deep copy of exposure cube instance.
 ***************************************************************************/
GCTACubeExposure* GCTACubeExposure::clone(void) const
{
    return new GCTACubeExposure(*this);
}


/***********************************************************************//**
 * @brief Set exposure cube for one CTA observation
 *
 * @param[in] obs CTA observation.
 *
 * Set the exposure cube for one CTA observation. The cube pixel values are
 * computed as product of the effective area and the livetime.
 ***************************************************************************/
void GCTACubeExposure::set(const GCTAObservation& obs)
{
    // Clear GTIs, reset livetime and exposure cube pixels
    m_gti.clear();
    m_livetime = 0.0;
    m_cube     = 0.0;

    // Fill exposure cube
    fill_cube(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill exposure cube from observation container
 *
 * @param[in] obs Observation container.
 * @param[in] log Pointer to logger.
 *
 * Set the exposure cube by summing the exposure for all CTA observations in
 * an observation container. The cube pixel values are computed as the sum
 * over the products of the effective area and the livetime.
 ***************************************************************************/
void GCTACubeExposure::fill(const GObservations& obs, GLog* log)
{
    // Clear GTIs, reset livetime and exposure cube pixels
    m_gti.clear();
    m_livetime = 0.0;
    m_cube     = 0.0;

    // Loop over all observations in container
    for (int i = 0; i < obs.size(); ++i) {

        // Get observation and continue only if it is a CTA observation
        const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>
                                     (obs[i]);

        // Skip observation if it's not CTA
        if (cta == NULL) {
            if (log != NULL) {
                *log << "Skipping ";
                *log << obs[i]->instrument();
                *log << " observation ";
                *log << "\"" << obs[i]->name() << "\"";
                *log << " (id=" << obs[i]->id() << ")" << std::endl;
            }
            continue;
        }

        // Skip observation if we have a binned observation
        if (cta->eventtype() == "CountsCube") {
            if (log != NULL) {
                *log << "Skipping binned ";
                *log << cta->instrument();
                *log << " observation ";
                *log << "\"" << cta->name() << "\"";
                *log << " (id=" << cta->id() << ")" << std::endl;
            }
            continue;
        }

        // Fill exposure cube
        fill_cube(*cta, log);

    } // endfor: looped over observations

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read exposure cube from FITS object
 *
 * @param[in] fits FITS object.
 *
 * Read the exposure cube from a FITS object.
 ***************************************************************************/
void GCTACubeExposure::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get HDUs
    const GFitsImage& hdu_expcube  = *fits.image("Primary");
    const GFitsTable& hdu_energies = *fits.table(gammalib::extname_energies);
    const GFitsTable& hdu_gti      = *fits.table(gammalib::extname_gti);

    // Read cube
    m_cube.read(hdu_expcube);

    // Read cube attributes
    read_attributes(hdu_expcube);

    // Read energies
    m_energies.read(hdu_energies);

    // Read GTIs
    m_gti.read(hdu_gti);

    // Set energy node array
    set_eng_axis();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA exposure cube into FITS file.
 *
 * @param[in] fits FITS file.
 *
 * Writes the exposure cube image, the energies and the Good Time Intervals
 * into the FITS file.
 ***************************************************************************/
void GCTACubeExposure::write(GFits& fits) const
{
    // Write cube
    m_cube.write(fits);

    // Get last HDU and write attributes
    if (fits.size() > 0) {
        GFitsHDU& hdu = *fits[fits.size()-1];
        write_attributes(hdu);
    }

    // Write energies
    m_energies.write(fits);

    // Write GTIs
    m_gti.write(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load exposure cube from FITS file
 *
 * @param[in] filename Performance table file name.
 *
 * Loads the exposure cube from a FITS file into the object.
 ***************************************************************************/
void GCTACubeExposure::load(const GFilename& filename)
{
    // Put into OpenMP criticial zone
    #pragma omp critical(GCTACubeExposure_load)
    {

    // Open FITS file
    GFits fits(filename);

    // Read PSF cube
    read(fits);

    // Close FITS file
    fits.close();

    } // end of OpenMP critical zone

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save exposure cube into FITS file
 *
 * @param[in] filename Exposure cube FITS file name.
 * @param[in] clobber Overwrite existing file?
 *
 * Save the exposure cube into a FITS file.
 ***************************************************************************/
void GCTACubeExposure::save(const GFilename& filename, const bool& clobber) const
{
    // Put into OpenMP criticial zone
    #pragma omp critical(GCTACubeExposure_save)
    {

    // Create FITS file
    GFits fits;

    // Write exposure cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Close Edisp file
    fits.close();

    } // end of OpenMP critical zone

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print exposure cube information
 *
 * @param[in] chatter Chattiness.
 * @return String containing exposure cube information.
 *
 * @todo Add content
 ***************************************************************************/
std::string GCTACubeExposure::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTACubeExposure ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Livetime"));
        result.append(gammalib::str(m_livetime)+" sec");

        // Append energies
        if (m_energies.size() > 0) {
            result.append("\n"+m_energies.print(chatter));
        }
        else {
            result.append("\n"+gammalib::parformat("Energies") +
                          "not defined");
        }

        // Append GTIs
        if (m_gti.size() > 0) {
            result.append("\n"+m_gti.print(chatter));
        }
        else {
            result.append("\n"+gammalib::parformat("Good Time Intervals") +
                          "not defined");
        }

        // Append skymap definition
        result.append("\n"+m_cube.print(chatter));

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTACubeExposure::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_cube.clear();
    m_energies.clear();
    m_elogmeans.clear();
    m_gti.clear();
    m_livetime = 0.0;

    // Initialise cache
    m_inx_left  = 0;
    m_inx_right = 0;
    m_wgt_left  = 0.0;
    m_wgt_right = 0.0;

    // Set CTA time reference for GTIs
    m_gti.reference(GTimeReference(G_CTA_MJDREF, "s", "TT", "LOCAL"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Exposure cube
 ***************************************************************************/
void GCTACubeExposure::copy_members(const GCTACubeExposure& cube)
{
    // Copy members
    m_filename  = cube.m_filename;
    m_cube      = cube.m_cube;
    m_energies  = cube.m_energies;
    m_elogmeans = cube.m_elogmeans;
    m_gti       = cube.m_gti;
    m_livetime  = cube.m_livetime;

    // Copy cache
    m_inx_left  = cube.m_inx_left;
    m_inx_right = cube.m_inx_right;
    m_wgt_left  = cube.m_wgt_left;
    m_wgt_right = cube.m_wgt_right;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTACubeExposure::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill exposure cube for one observation
 *
 * @param[in] obs Observation.
 * @param[in] log Pointer to logger.
 *
 * @exception GException::invalid_value
 *            No RoI or response found in CTA observation.
 *
 * Fill the exposure cube from one CTA observations. The exposure cube pixel
 * values are computed as the product of the effective area and the livetime.
 ***************************************************************************/
void GCTACubeExposure::fill_cube(const GCTAObservation& obs, GLog* log)
{
    // Set energy margin
    static const GEnergy margin(1.0, "MeV");

    // Only continue if we have an event list
    if (obs.eventtype() == "EventList") {

        // Extract pointing direction, energy boundaries and ROI from
        // observation
        GSkyDir  pnt         = obs.pointing().dir();
        GEbounds obs_ebounds = obs.ebounds();
        GCTARoi  roi         = obs.roi();

        // Check for RoI sanity
        if (!roi.is_valid()) {
            std::string msg = "No RoI information found in "+obs.instrument()+
                              " observation \""+obs.name()+"\". Run ctselect "
                              "to specify an RoI for this observation.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Convert RoI into a circular region for overlap checking
        GSkyRegionCircle roi_reg(roi.centre().dir(), roi.radius());

        // Extract response from observation
        const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>
                                     (obs.response());
        if (rsp == NULL) {
            std::string msg = "No valid instrument response function found in "+
                              obs.instrument()+" observation \""+obs.name()+
                              "\". Please specify the instrument response "
                              "function for this observation.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Skip observation if livetime is zero
        if (obs.livetime() == 0.0) {
            if (log != NULL) {
                *log << "Skipping unbinned ";
                *log << obs.instrument();
                *log << " observation ";
                *log << "\"" << obs.name() << "\"";
                *log << " (id=" << obs.id() << ") due to zero livetime";
                *log << std::endl;
            }
        }

        // Skip observation if observation is outside the bounds of the cube
        else if (!m_cube.overlaps(roi_reg)) {
            if (log != NULL) {
                *log << "Skipping unbinned ";
                *log << obs.instrument();
                *log << " observation ";
                *log << "\"" << obs.name() << "\"";
                *log << " (id=" << obs.id() << ") since it does not overlap ";
                *log << "with the exposure cube.";
                *log << std::endl;
            }
        }

        // ... otherwise continue
        else {

            // Announce observation usage
            if (log != NULL) {
                *log << "Including ";
                *log << obs.instrument();
                *log << " observation \"" << obs.name();
                *log << "\" (id=" << obs.id() << ")";
                *log << " in exposure cube computation." << std::endl;
            }

            // Loop over all pixels in sky map
            for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

                // Get pixel sky direction
                GSkyDir dir = m_cube.inx2dir(pixel);

                // Skip pixel if it is outside the RoI
                if (roi.centre().dir().dist_deg(dir) > roi.radius()) {
                    continue;
                }

                // Compute theta angle with respect to pointing direction in
                // radians
                double theta = pnt.dist(dir);

                // Loop over all exposure cube energies
                for (int iebin = 0; iebin < m_energies.size(); ++iebin){

                    // Skip exposure cube energy if the energy is outside the
                    // observation energy. The exposure cube energies are true
                    // energies while the observation energy boundaries are
                    // reconstructed energies, hence this is only an
                    // approximation, but probably the only we can really do.
                    // We allow here for a small margin in case of rounding
                    // errors in the energy boundaries.
                    #if defined(G_CLIP_ETRUES)
                    if (!(obs_ebounds.contains(m_energies[iebin])        ||
                          obs_ebounds.contains(m_energies[iebin]-margin) ||
                          obs_ebounds.contains(m_energies[iebin]+margin))) {
                        continue;
                    }
                    #endif

                    // Get logE/TeV
                    double logE = m_energies[iebin].log10TeV();

                    // Add to exposure cube (effective area * livetime)
                    m_cube(pixel, iebin) += rsp->aeff(theta, 0.0, 0.0, 0.0, logE) *
                                            obs.livetime();

                } // endfor: looped over energy bins

            } // endfor: looped over all pixels

            // If GTI is empty then set its time reference from the
            // observation. From then on we keep that time reference
            if (m_gti.is_empty()) {
                m_gti.reference(obs.gti().reference());
            }

            // Append GTIs and increment livetime
            m_gti.extend(obs.gti());
            m_livetime += obs.livetime();

        } // endelse: livetime was not zero

    } // endif: observation contained an event list

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update 1D cache
 *
 * @param[in] logE Log10 energy in TeV.
 *
 * Updates the 1D interpolation cache. The interpolation cache is composed
 * of two indices and weights that define 2 data values of the 2D skymap
 * that are used for linear interpolation.
 *
 * @todo Write down formula
 ***************************************************************************/
void GCTACubeExposure::update(const double& logE) const
{
    // Set value for node array
    m_elogmeans.set_value(logE);

    // Set indices and weighting factors for interpolation
    m_inx_left  = m_elogmeans.inx_left();
    m_inx_right = m_elogmeans.inx_right();
    m_wgt_left  = m_elogmeans.wgt_left();
    m_wgt_right = m_elogmeans.wgt_right();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for a logarithmic (base 10) energy axis
 *
 * Set axis nodes so that each node is the logarithm of the energy values.
 ***************************************************************************/
void GCTACubeExposure::set_eng_axis(void)
{
    // Get number of bins
    int bins = m_energies.size();

    // Clear node array
    m_elogmeans.clear();

    // Compute nodes
    for (int i = 0; i < bins; ++i) {

        // Get logE/TeV
        m_elogmeans.append(m_energies[i].log10TeV());

    }  // endfor: looped over energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read exposure attributes
 *
 * @param[in] hdu FITS HDU.
 *
 * Reads CTA exposure attributes from the HDU.
 ***************************************************************************/
void GCTACubeExposure::read_attributes(const GFitsHDU& hdu)
{
    // Read mandatory attributes
    m_livetime = (hdu.has_card("LIVETIME")) ? hdu.real("LIVETIME") : 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write attributes to exposure extension
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GCTACubeExposure::write_attributes(GFitsHDU& hdu) const
{
    // Compute some attributes
    double      tstart   = m_gti.tstart().convert(m_gti.reference());
    double      tstop    = m_gti.tstop().convert(m_gti.reference());
    double      telapse  = m_gti.telapse();
    double      ontime   = m_gti.ontime();
    double      deadc    = (ontime > 0.0 && m_livetime > 0.0) ? 
                           m_livetime / ontime : 1.0;
    std::string utc_obs  = m_gti.tstart().utc();
    std::string utc_end  = m_gti.tstop().utc();
    std::string date_obs = utc_obs.substr(0, 10);
    std::string time_obs = utc_obs.substr(11, 8);
    std::string date_end = utc_end.substr(0, 10);
    std::string time_end = utc_end.substr(11, 8);

    // Set observation information
    hdu.card("CREATOR",  "GammaLib", "Program which created the file");
    hdu.card("TELESCOP", "unknown",  "Telescope");
    hdu.card("OBS_ID",   "unknown",  "Observation identifier");
    hdu.card("DATE-OBS", date_obs,   "Observation start date");
    hdu.card("TIME-OBS", time_obs,   "Observation start time");
    hdu.card("DATE-END", date_end,   "Observation end date");
    hdu.card("TIME-END", time_end,   "Observation end time");

    // Set observation time information
    hdu.card("TSTART",   tstart, "[s] Mission time of start of observation");
    hdu.card("TSTOP",    tstop, "[s] Mission time of end of observation");
    m_gti.reference().write(hdu);
    hdu.card("TELAPSE",  telapse, "[s] Mission elapsed time");
    hdu.card("ONTIME",   ontime, "[s] Total good time including deadtime");
    hdu.card("LIVETIME", m_livetime, "[s] Total livetime");
    hdu.card("DEADC",    deadc, "Deadtime correction factor");
    hdu.card("TIMEDEL",  1.0, "Time resolution");

    // Return
    return;
}
