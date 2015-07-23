/***************************************************************************
 *          GCTACubeExposure.cpp - CTA cube analysis exposure class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Chia-Chun Lu                                *
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
#include "GCTACubeExposure.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAEventList.hpp"
#include "GMath.hpp"
#include "GTools.hpp"
#include "GLog.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SET                       "GCTACubeExposure::set(GCTAObservation&)"
#define G_FILL                "GCTACubeExposure::fill(GObservations&, GLog*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

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
GCTACubeExposure::GCTACubeExposure(const std::string& filename)
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
    m_ebounds = cube.ebounds();

    // Set GNodeArray used for interpolation
    set_eng_axis();

    // Set exposure cube to event cube
    m_cube = cube.map();

    // Set all exposure cube pixels to zero as we want to have a clean map
    // upon construction
    m_cube = 0.0;

    // Return
    return;

}


/***********************************************************************//**
 * @brief Exposure cube constructor
 *
 * @param[in] wcs     World Coordinate System.
 * @param[in] coords  Coordinate System (CEL or GAL).
 * @param[in] x       X coordinate of sky map centre (deg).
 * @param[in] y       Y coordinate of sky map centre (deg).
 * @param[in] dx      Pixel size in x direction at centre (deg/pixel).
 * @param[in] dy      Pixel size in y direction at centre (deg/pixel).
 * @param[in] nx      Number of pixels in x direction.
 * @param[in] ny      Number of pixels in y direction.
 * @param[in] ebounds Energy boundaries.
 *
 * Constructs an exposure cube by specifying the sky map grid and the energy
 * boundaries.
 ***************************************************************************/
GCTACubeExposure::GCTACubeExposure(const std::string&   wcs,
                                   const std::string&   coords,
                                   const double&        x,
                                   const double&        y,
                                   const double&        dx,
                                   const double&        dy,
                                   const int&           nx,
                                   const int&           ny,
                                   const GEbounds&      ebounds)
{
    // Initialise class members
    init_members();

    // Store energy boundaries
    m_ebounds = ebounds;

    // Set GNodeArray used for interpolation
    set_eng_axis();

    // Create sky map
    m_cube = GSkyMap(wcs, coords, x, y, dx, dy, nx, ny, m_ebounds.size());

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
 * @brief Set exposure cube from one CTA observation
 *
 * @param[in] obs CTA observation.
 *
 * Set the exposure cube from a single CTA observations. The cube pixel
 * values are computed as product of the effective area and the livetime.
 *
 * @todo: Throw an exception if response is not valid
 ***************************************************************************/
void GCTACubeExposure::set(const GCTAObservation& obs)
{

    // Only continue if we have an unbinned observation
    if (obs.eventtype() == "EventList") {

        // Clear GTIs, reset livetime and exposure cube pixels
        m_gti.clear();
        m_livetime = 0.0;
        m_cube     = 0.0;

        // Extract region of interest from CTA observation
        GCTARoi roi = obs.roi();

        // Check for RoI sanity
        if (!roi.is_valid()) {
            std::string msg = "No RoI information found in input observation "
                              "\""+obs.name()+"\". Run ctselect to specify "
                              "an RoI for this observation";
            throw GException::invalid_value(G_SET, msg);
        }

        // Get references on CTA response and pointing direction
        const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(obs.response());
        const GSkyDir&         pnt = obs.pointing().dir();

        // Continue only if response is valid
        if (rsp != NULL) {

            // Loop over all pixels in sky map
            for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

                // Get pixel sky direction
                GSkyDir dir = m_cube.inx2dir(pixel);

                // Continue only if pixel is within RoI
                if (roi.centre().dir().dist_deg(dir) <= roi.radius()) {

                    // Compute theta angle with respect to pointing direction
                    // in radians
                    double theta = pnt.dist(dir);

                    // Loop over all exposure cube energy bins
                    for (int iebin = 0; iebin < m_ebounds.size(); ++iebin){

                        // Get logE/TeV
                        double logE = m_ebounds.elogmean(iebin).log10TeV();

                        // Set exposure cube (effective area * lifetime)
                        m_cube(pixel, iebin) = rsp->aeff(theta, 0.0, 0.0, 0.0, logE) *
                                               obs.livetime();

                    } // endfor: looped over energy bins

                } // endif: pixel was within RoI
    
            } // endfor: looped over all pixels
    
            // Append GTIs and increment livetime
            m_gti.extend(obs.gti());
            m_livetime += obs.livetime();

        } // endif: response was valid

    } // endif: observation was unbinned

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill exposure cube from observation container
 *
 * @param[in] obs Observation container.
 * @param[in] log Pointer to logger (optional).
 *
 * @exception GException::invalid_value
 *            No event list found in CTA observations.
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
        const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(obs[i]);

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

        // Extract region of interest from CTA observation
        GCTARoi roi = cta->roi();

        // Extract energy boundaries from CTA observation
        GEbounds obs_ebounds = cta->ebounds();

        // Check for RoI sanity
        if (!roi.is_valid()) {
            std::string msg = "No RoI information found in input observation "
                              "\""+cta->name()+"\". Run ctselect to specify "
                              "an RoI for this observation";
            throw GException::invalid_value(G_FILL, msg);
        }

        // Get references on CTA response and pointing direction
        const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>(cta->response());
        const GSkyDir&         pnt = cta->pointing().dir();

        // Skip observation if we don't have an unbinned observation
        if (rsp == NULL) {
            if (log != NULL) {
                *log << "WARNING: ";
                *log << cta->instrument();
                *log << " observation \"" << cta->name();
                *log << "\" (id=" << cta->id() << ")";
                *log << " contains no IRF response.";
                *log << " Skipping this observation." << std::endl;
            }
            continue;
        }

        // Announce observation usage
        if (log != NULL) {
            *log << "Including ";
            *log << cta->instrument();
            *log << " observation \"" << cta->name();
            *log << "\" (id=" << cta->id() << ")";
            *log << " in exposure cube computation." << std::endl;
        }
            
        // Loop over all pixels in sky map
        for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

            // Get pixel sky direction
            GSkyDir dir = m_cube.inx2dir(pixel);
                    
            // Continue only if pixel is within RoI
            if (roi.centre().dir().dist_deg(dir) <= roi.radius()) {

                // Compute theta angle with respect to pointing
                // direction in radians
                double theta = pnt.dist(dir);

                // Loop over all exposure cube energy bins
                for (int iebin = 0; iebin < m_ebounds.size(); ++iebin){

                    // Only add exposure to pixel if energy bin is inside 
                    // observation energy boundaries
                    if (obs_ebounds.contains(m_ebounds.emin(iebin),
                                             m_ebounds.emax(iebin))) {

                        // Get logE/TeV
                        double logE = m_ebounds.elogmean(iebin).log10TeV();

                        // Add to exposure cube (effective area * livetime)
                        m_cube(pixel, iebin) += rsp->aeff(theta, 0.0, 0.0, 0.0, logE) *
                                                cta->livetime();

                    } // endif: energy bin was inside observation energy boundaries

                } // endfor: looped over energy bins

            } // endif: pixel within RoI

        } // endfor: looped over all pixels

        // Append GTIs and increment livetime
        m_gti.extend(cta->gti());
        m_livetime += cta->livetime();
            
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
    const GFitsImage& hdu_expcube = *fits.image("Primary");
    const GFitsTable& hdu_ebounds = *fits.table("EBOUNDS");
    const GFitsTable& hdu_gti     = *fits.table("GTI");

    // Read cube
    m_cube.read(hdu_expcube);

    // Read cube attributes
    read_attributes(hdu_expcube);

    // Read energy boundaries
    m_ebounds.read(hdu_ebounds);

    // Read GTIs
    m_gti.read(hdu_gti);

    // Set energy node array
    set_eng_axis();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA exposure cube into FITS object.
 *
 * @param[in] fits FITS file.
 *
 * Writes the exposure cube image, the energy boundaries and the Good Time
 * Intervals into the FITS object.
 ***************************************************************************/
void GCTACubeExposure::write(GFits& fits) const
{
    // Write cube
    m_cube.write(fits);

    // Get last HDU and write attributes
    GFitsHDU& hdu = *fits[fits.size()-1];
    write_attributes(hdu);

    // Write energy boundaries
    m_ebounds.write(fits);

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
void GCTACubeExposure::load(const std::string& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read PSF cube
    read(fits);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save exposure cube into FITS file
 *
 * @param[in] filename Exposure cube FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * Save the exposure cube into a FITS file.
 *
 * @todo Implement method
 ***************************************************************************/
void GCTACubeExposure::save(const std::string& filename, const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write exposure cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print exposure cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
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

        // Append energy intervals
        if (m_ebounds.size() > 0) {
            result.append("\n"+m_ebounds.print(chatter));
        }
        else {
            result.append("\n"+gammalib::parformat("Energy intervals") +
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
    m_ebounds.clear();
    m_elogmeans.clear();
    m_gti.clear();
    m_livetime = 0.0;

    // Initialise cache
    m_inx_left  = 0;
    m_inx_right = 0;
    m_wgt_left  = 0.0;
    m_wgt_right = 0.0;
   
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
    m_ebounds   = cube.m_ebounds;
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
 *
 * Set axis nodes so that each node is the logarithmic mean of the lower and
 * upper energy boundary, i.e.
 * \f[ n_i = \log \sqrt{{\rm LO}_i \times {\rm HI}_i} \f]
 * where
 * \f$n_i\f$ is node \f$i\f$,
 * \f${\rm LO}_i\f$ is the lower bin boundary for bin \f$i\f$, and
 * \f${\rm HI}_i\f$ is the upper bin boundary for bin \f$i\f$.
 *
 * @todo Check that none of the axis boundaries is non-positive.
 ***************************************************************************/
void GCTACubeExposure::set_eng_axis(void)
{
    // Get number of bins
    int bins = m_ebounds.size();

    // Clear node array
    m_elogmeans.clear();

    // Compute nodes
    for (int i = 0; i < m_ebounds.size(); ++i) {
     
        // Get logE/TeV
        m_elogmeans.append(m_ebounds.elogmean(i).log10TeV()); 

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
    hdu.card("DATE_OBS", date_obs,   "Observation start date");
    hdu.card("TIME_OBS", time_obs,   "Observation start time");
    hdu.card("DATE_END", date_end,   "Observation end date");
    hdu.card("TIME_END", time_end,   "Observation end time");

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
