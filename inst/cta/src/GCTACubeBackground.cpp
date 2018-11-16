/***************************************************************************
 *             GCTACubeBackground.cpp - CTA cube background class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2018 by Michael Mayer                               *
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
 * @file GCTACubeBackground.cpp
 * @brief CTA cube background class implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GLog.hpp"
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GObservations.hpp"
#include "GSkyRegionCircle.hpp"
#include "GCTAEventCube.hpp"
#include "GCTAObservation.hpp"
#include "GCTARoi.hpp"
#include "GCTAInstDir.hpp"
#include "GCTACubeBackground.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                             "GCTACubeBackground::read(GFits&)"
#define G_MC                "GCTACubeBackground::mc(GEnergy&, GTime&, GRan&)"
#define G_FILL              "GCTACubeBackground::fill(GObservations&, GLog*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_LOG_INTERPOLATION   //!< Energy interpolate log10(background rate)

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_MC_INIT
//#define G_DEBUG_CACHE

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTACubeBackground::GCTACubeBackground(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename FITS file name.
 *
 * Construct instance by loading the background information from a FITS file.
 ***************************************************************************/
GCTACubeBackground::GCTACubeBackground(const GFilename& filename)
{
    // Initialise class members
    init_members();

    // Load background from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Event cube constructor
 *
 * @param[in] cube Event cube.
 *
 * Construct background cube using the same binning and sky projection that
 * is used for the event cube.
 ***************************************************************************/
GCTACubeBackground::GCTACubeBackground(const GCTAEventCube& cube)
{
    // Initialise class members
    init_members();

    // Store energy boundaries
    m_ebounds = cube.ebounds();

    // Set GNodeArray used for interpolation
    set_eng_axis();

    // Set background cube to event cube
    m_cube = cube.counts();
    m_cube.nmaps(m_ebounds.size());

    // Set all background cube pixels to zero as we want to have a clean map
    // upon construction
    m_cube = 0.0;

    // Return
    return;

}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bgd Background.
 ***************************************************************************/
GCTACubeBackground::GCTACubeBackground(const GCTACubeBackground& bgd)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(bgd);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Background cube constructor
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
 * Constructs a background cube by specifying the sky map grid and the
 * energies.
 ***************************************************************************/
GCTACubeBackground::GCTACubeBackground(const std::string&   wcs,
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
GCTACubeBackground::~GCTACubeBackground(void)
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
 * @param[in] bgd Background cube.
 * @return Background cube.
 ***************************************************************************/
GCTACubeBackground& GCTACubeBackground::operator=(const GCTACubeBackground& bgd)
{
    // Execute only if object is not identical
    if (this != &bgd) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(bgd);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return background rate in units of events MeV\f$^{-1}\f$
 *        s\f$^{-1}\f$ sr\f$^{-1}\f$
 *
 * @param[in] dir Reconstructed event position where rate should be computed
 * @param[in] energy Energy at which the rate should be computed
 * @return Background rate (events MeV\f$^{-1}\f$ s\f$^{-1}\f$ sr\f$^{-1}\f$).
 *
 * Returns the background rate for a giben instrument direction and energy.
 * The rate is given per ontime. The method assures that the background rate
 * never becomes negative.
 *
 * The method interpolates logarithmically in the energy direction.
 ***************************************************************************/
double GCTACubeBackground::operator()(const GCTAInstDir& dir,
                                      const GEnergy&     energy) const
{
    // Initialise background rate
    double background = 0.0;

    // Set indices and weighting factors for interpolation
    update(energy.log10TeV());

    // If left weight is close to 1, use left node
    if (std::abs(m_wgt_left-1.0) < 1.0e-6) {
        background = m_cube(dir.dir(), m_inx_left);
    }

    // ... otherwise, if right weight is close to 1, use right node
    else if (std::abs(m_wgt_right-1.0) < 1.0e-6) {
        background = m_cube(dir.dir(), m_inx_right);
    }

    // ... otherwise perform interpolation
    else {
        double background_left  = m_cube(dir.dir(), m_inx_left);
        double background_right = m_cube(dir.dir(), m_inx_right);
        if (background_left > 0.0 && background_right > 0.0) {
            background = std::exp(m_wgt_left  * std::log(background_left) +
                                  m_wgt_right * std::log(background_right));
        }
    }

    // Make sure that background rate does not become negative
    if (background < 0.0) {
        background = 0.0;
    }

    // Return background rate
    return background;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear background
 *
 * This method properly resets the background to an initial state.
 ***************************************************************************/
void GCTACubeBackground::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone background
 *
 * @return Pointer to deep copy of background.
 ***************************************************************************/
GCTACubeBackground* GCTACubeBackground::clone(void) const
{
    return new GCTACubeBackground(*this);
}


/***********************************************************************//**
 * @brief Fill background cube from observation container
 *
 * @param[in] obs Observation container.
 * @param[in] log Pointer to logger (optional).
 *
 * @exception GException::invalid_value
 *            No event list found in CTA observations.
 *
 * Set the background cube by computing the ontime weighted background rate
 * for all CTA observations in an observation container. The cube pixel
 * values are computed as the sum over the background rates.
 ***************************************************************************/
void GCTACubeBackground::fill(const GObservations& obs, GLog* log)
{
    // Set energy margin
    const GEnergy energy_margin(0.01, "GeV");

    // Clear background cube
    m_cube = 0.0;

    // Set dummy GTI needed to generate an event cube. It is not important
    // what the actually value is since it will be overwritten later in any
    // case, but it's important that there is one time slice
    GGti gti(GTime(0.0), GTime(1.0));

    // Initialise event cubes for model evaluate and storage
    GCTAEventCube eventcube = GCTAEventCube(m_cube, m_ebounds, gti);

    // Initialise total ontime
    double total_ontime = 0.0;

    // Loop over all observations in container
    for (int i = 0; i < obs.size(); ++i) {

        // Get observation and continue only if it is a CTA observation
        const GCTAObservation* cta =
              dynamic_cast<const GCTAObservation*>(obs[i]);

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

        // Get observation ontime
        double ontime = cta->ontime();

        // Skip observation if ontime is zero
        if (ontime == 0.0) {
            if (log != NULL) {
                *log << "Skipping unbinned ";
                *log << cta->instrument();
                *log << " observation ";
                *log << "\"" << cta->name() << "\"";
                *log << " (id=" << cta->id() << ") due to zero ontime";
                *log << std::endl;
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

        // Make sure the observation overlaps with the cube
        GSkyRegionCircle obs_reg(roi.centre().dir(), roi.radius());
        if (!m_cube.overlaps(obs_reg)) {
            if (log != NULL) {
                *log << "Skipping unbinned ";
                *log << cta->instrument();
                *log << " observation ";
                *log << "\"" << cta->name() << "\"";
                *log << " (id=" << cta->id() << ") since it does not overlap ";
                *log << "with the background cube.";
                *log << std::endl;
            }
            continue;
        } // endif: m_cube overlaps observation

        // Announce observation usage
        if (log != NULL) {
            *log << "Including ";
            *log << cta->instrument();
            *log << " observation \"" << cta->name();
            *log << "\" (id=" << cta->id() << ")";
            *log << " in background cube computation." << std::endl;
        }

        // Set GTI of actual observations as the GTI of the event cube
        eventcube.gti(cta->gti());

        // Compute number of pixels and energy layers
        int npix   = eventcube.npix();
        int nebins = m_ebounds.size();

        // Loop over all spatial pixels
        for (int ipix = 0; ipix < npix; ++ipix) {

            // Get event bin in first layer
            GCTAEventBin* bin = eventcube[ipix];

            // Skip if bin is not contained in RoI
            if (!roi.contains(*bin)) {
                continue;
            }

            // Get CTA direction for this bin, including sky direction and
            // DETX/DETY coordinates
            GCTAInstDir dir = cta->pointing().instdir(bin->dir().dir());

            // Initialise values for integration over energy bin
            double f1 = 0.0;
            double f2 = 0.0;

            // Loop over all energy bins of the cube
            for (int iebin = 0, ibin = ipix; iebin < nebins; ++iebin, ibin += npix) {

                // If this is the first bin then compute f1 at the lower bin edge
                if (iebin == 0) {
                    GCTAEventBin* bin = eventcube[ibin];
                    bin->dir(dir); // Set instrument direction
                    bin->energy(m_ebounds.emin(iebin));
                    f1 = obs.models().eval(*bin, *cta);
                }

                // Compute f2 at the upper bin edge
                GCTAEventBin* bin = eventcube[ibin];
                bin->dir(dir); // Set instrument direction
                bin->energy(m_ebounds.emax(iebin));
                f2 = obs.models().eval(*bin, *cta);

                // Add background model value to cube if the energy bin is
                // fully container in the observation interval and if both f1
                // and f2 are positive. Add a small margin here to the energy
                // bin boundaries to take provision for rounding errors.
                if (obs_ebounds.contains(m_ebounds.emin(iebin) + energy_margin,
                                         m_ebounds.emax(iebin) - energy_margin) &&
                    (f1 > 0.0 && f2 > 0.0)) {

                    // Compute background model value by integrating over the
                    // bin and deviding the result by the bin width
                    double x1    = m_ebounds.emin(iebin).MeV();
                    double x2    = m_ebounds.emax(iebin).MeV();
                    double model = gammalib::plaw_integral(x1, f1, x2, f2) / (x2 - x1);

                    // Multiply by ontime to get the correct weighting for each
                    // observation. We divide by the total ontime later to get
                    // the background model in units of counts/(MeV s sr).
                    model *= ontime;

                    // Add existing number of counts
                    model += bin->counts();

                    // Store cumulated value (units: counts/(MeV sr))
                    bin->counts(model);

                } // endif: added background model value to cube

                // Set f1 to f2 for next energy layer
                f1 = f2;

            } // endfor: looped over energy layers

        } // endfor: looped over spatial pixels

        // Accumulate ontime
        total_ontime += ontime;

    } // endfor: looped over all observations

    // Re-normalize cube to get units of counts/(MeV s sr)
    if (total_ontime > 0.0) {

        // Loop over all bins in background cube and divide the content
        // by the total ontime.
        for (int i = 0; i < eventcube.size(); ++i) {
            GCTAEventBin* bin  = eventcube[i];
            double        rate = bin->counts() / total_ontime;
            bin->counts(rate);
        }

        // Set background cube values from event cube
        m_cube = eventcube.counts();

    } // endif: ontime was positive

    // Return
    return;

}


/***********************************************************************//**
 * @brief Return spatially integrated background rate in units of
 *        events MeV\f$^{-1}\f$ s\f$^{-1}\f$
 *
 * @param[in] logE Logarithm (base 10) of energy in TeV
 * @return Spatially integrated background rate
 *         (events MeV\f$^{-1}\f$ s\f$^{-1}\f$)
 *
 * Spatially integrates the background cube at a given energy. This method
 * performs an interpolation between the energy maps. The integration is
 * performed by summing over all bin contents multiplied by the solid angle
 * of the bin.
 ***************************************************************************/
double GCTACubeBackground::integral(const double& logE) const
{
    // Update interpolation cache
    update(logE);

    // Initialise result
    double result = 0.0;

    // Loop over all map pixels
    for (int i = 0; i < m_cube.npix(); ++i) {

        // Get bin value
        double value = m_wgt_left  * m_cube(i, m_inx_left) +
                       m_wgt_right * m_cube(i, m_inx_right);

        // Sum bin contents
        result += value * m_cube.solidangle(i);

    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Read background cube from FITS object
 *
 * @param[in] fits FITS object.
 *
 * Read the background cube from a FITS object.
 ***************************************************************************/
void GCTACubeBackground::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get HDUs
    const GFitsImage& hdu_bgdcube = *fits.image("Primary");
    const GFitsTable& hdu_ebounds = *fits.table(gammalib::extname_ebounds);

    // Read cube
    m_cube.read(hdu_bgdcube);

    // Read energy boundaries
    m_ebounds.read(hdu_ebounds);

    // Set energy node array
    set_eng_axis();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA background cube into FITS object.
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GCTACubeBackground::write(GFits& fits) const
{
    // Write cube
    m_cube.write(fits);

    // Write energy boundaries
    m_ebounds.write(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load background cube from FITS file
 *
 * @param[in] filename FITS file.
 *
 * Loads the background cube from a FITS file.
 ***************************************************************************/
void GCTACubeBackground::load(const GFilename& filename)
{
    // Put into OpenMP criticial zone
    #pragma omp critical(GCTACubeBackground_load)
    {

    // Open FITS file
    GFits fits(filename);

    // Read background cube
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
 * @brief Save background cube into FITS file
 *
 * @param[in] filename background cube FITS file name.
 * @param[in] clobber Overwrite existing file?
 *
 * Save the background cube into a FITS file.
 ***************************************************************************/
void GCTACubeBackground::save(const GFilename& filename,
                              const bool&      clobber) const
{
    // Put into OpenMP criticial zone
    #pragma omp critical(GCTACubeBackground_save)
    {

    // Open or create FITS file
    GFits fits;

    // Write background cube
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
 * @brief Print background information
 *
 * @param[in] chatter Chattiness.
 * @return String containing background information.
 ***************************************************************************/
std::string GCTACubeBackground::print(const GChatter& chatter) const
{
    // Initialise result string
     std::string result;

     // Continue only if chatter is not silent
     if (chatter != SILENT) {

         // Append header
         result.append("=== GCTACubeBackground ===");

         // Append parameters
         result.append("\n"+gammalib::parformat("Filename")+m_filename);

        // Append energies
        if (m_ebounds.size() > 0) {
            result.append("\n"+m_ebounds.print(chatter));
        }
        else {
            result.append("\n"+gammalib::parformat("Energy boundaries") +
                          "Not defined");
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
void GCTACubeBackground::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_cube.clear();
    m_ebounds.clear();
    m_elogmeans.clear();

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
 * @param[in] bgd Background.
 ***************************************************************************/
void GCTACubeBackground::copy_members(const GCTACubeBackground& bgd)
{
    // Copy members
    m_filename  = bgd.m_filename;
    m_cube      = bgd.m_cube;
    m_ebounds   = bgd.m_ebounds;
    m_elogmeans = bgd.m_elogmeans;

    // Copy cache
    m_inx_left  = bgd.m_inx_left;
    m_inx_right = bgd.m_inx_right;
    m_wgt_left  = bgd.m_wgt_left;
    m_wgt_right = bgd.m_wgt_right;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTACubeBackground::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for a logarithmic (base 10) energy axis
 *
 *
 * Set axis nodes so that each node is the logarithm of the energy values.
 ***************************************************************************/
void GCTACubeBackground::set_eng_axis(void)
{
    // Get number of bins
    int bins = m_ebounds.size();

    // Clear node array
    m_elogmeans.clear();

    // Compute nodes
    for (int i = 0; i < bins; ++i) {

        // Get logE/TeV
        m_elogmeans.append(m_ebounds.elogmean(i).log10TeV());

    }  // endfor: looped over energy bins

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
void GCTACubeBackground::update(const double& logE) const
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
