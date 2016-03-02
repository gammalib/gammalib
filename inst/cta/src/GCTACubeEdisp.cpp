/***************************************************************************
 *     GCTACubeEdisp.cpp - CTA cube analysis energy dispersion function class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Michael Mayer                                *
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
 * @file GCTACubeEdisp.cpp
 * @brief CTA cube analysis energy dispersion class implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GLog.hpp"
#include "GObservations.hpp"
#include "GCTACubeEdisp.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventCube.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SET                            "GCTACubeEdisp::set(GCTAObservation&)"
#define G_FILL                     "GCTACubeEdisp::fill(GObservations&, GLog*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_SMOOTH_EDISP                      //!< Guarantee no singularities

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
const GEnergy g_energy_margin(1.0e-12, "TeV");


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Energy dispersion.
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(const GCTACubeEdisp& cube)
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
 * @param[in] filename energy dispersion cube filename.
 *
 * Construct energy dispersion cube by loading the information from a energy dispersion cube file.
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(const GFilename& filename)
{
    // Initialise class members
    init_members();

    // Load Edisp cube from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Event cube constructor
 *
 * @param[in] cube Event cube.
 * @param[in] mmax Maximum migra.
 * @param[in] nmbins Number of migra bins.
 *
 * Construct Edisp cube using the same binning and sky projection that is
 * used for the event cube.
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(const GCTAEventCube& cube, const double& mmax,
                         const int& nmbins)
{
    // Initialise class members
    init_members();

    // Store energy boundaries
    m_ebounds = cube.ebounds();

    // Set GNodeArray used for interpolation
    set_eng_axis();

    // Set delta node array
    m_migras.clear();
    for (int i = 0; i < nmbins; ++i) {
        double binsize = mmax / double(nmbins);
        double migra   = binsize * (double(i) + 0.5); // avoid central singularity
        m_migras.append(migra);
    }

    // Set migra node array for computation
    set_migra_axis();

    // Compute number of sky maps
    int nmaps = m_ebounds.size() * m_migras.size();

    // Set Edisp cube to event cube
    m_cube = cube.map();

    // Set appropriate number of skymaps
    m_cube.nmaps(nmaps);

    // Set all Edisp cube pixels to zero as we want to have a clean map
    // upon construction
    m_cube = 0.0;

    // Return
    return;

}


/***********************************************************************//**
 * @brief Mean PSF cube constructor
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
 * @param[in] mmax    Maximum migra (deg).
 * @param[in] nmbins  Number of migra bins.
 *
 * Constructs a mean Edisp cube from cube parameters and migra binning
 *
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(const std::string&   wcs,
                         const std::string&   coords,
                         const double&        x,
                         const double&        y,
                         const double&        dx,
                         const double&        dy,
                         const int&           nx,
                         const int&           ny,
                         const GEbounds&      ebounds,
                         const double&        mmax,
                         const int&           nmbins)
{
    // Initialise class members
    init_members();

    // Store energy boundaries
    m_ebounds = ebounds;

    // Set energy node array
    set_eng_axis();

    // Set migra node array
    m_migras.clear();
    for (int i = 0; i < nmbins; ++i) {

        double binsize = mmax / double(nmbins);
        double migra   = binsize * (double(i) + 0.5); // avoid central singularity

        m_migras.append(migra);
    }

    // Set migra node array for computation
    set_migra_axis();

    // Compute number of sky maps
    int nmaps = m_ebounds.size() * m_migras.size();
    
    // Create sky map
    m_cube = GSkyMap(wcs, coords, x, y, dx, dy, nx, ny, nmaps);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTACubeEdisp::~GCTACubeEdisp(void)
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
 * @param[in] cube Edisp cube.
 * @return Mean Edisp cube.
 ***************************************************************************/
GCTACubeEdisp& GCTACubeEdisp::operator= (const GCTACubeEdisp& cube)
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
 * @brief Return edisp (in units or MeV^-1)
 *
 * @param[in] dir Coordinate of the true photon position.
 * @param[in] migra Fraction of true photon energy and reconstructed photon energy
 * @param[in] energy Energy of the true photon.
 * @return energy dispersion (in units or MeV^-1)
 *
 * Returns the energy dispersion for a given fraction of true and reconstructed
 * energy (in units of MeV^-1 for a given energy and coordinate
 ***************************************************************************/
double GCTACubeEdisp::operator()(const GSkyDir& dir,
                               const double&  migra,
                               const GEnergy& energy) const
{
    // Update indices and weighting factors for interpolation
    update(migra, energy.log10TeV());

    // Perform bi-linear interpolation
    double edisp = m_wgt1 * m_cube(dir, m_inx1) +
                 m_wgt2 * m_cube(dir, m_inx2) +
                 m_wgt3 * m_cube(dir, m_inx3) +
                 m_wgt4 * m_cube(dir, m_inx4);

    // Make sure that Edisp does not become negative
    if (edisp < 0.0) {
        edisp = 0.0;
    }

    // Return Edisp
    return edisp;
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
void GCTACubeEdisp::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of mean Edisp instance.
 ***************************************************************************/
GCTACubeEdisp* GCTACubeEdisp::clone(void) const
{
    return new GCTACubeEdisp(*this);
}


/***********************************************************************//**
 * @brief Set Edisp cube from one CTA observation
 *
 * @param[in] obs CTA observation.
 ***************************************************************************/
void GCTACubeEdisp::set(const GCTAObservation& obs)
{
    // Clear Edisp cube
    clear_cube();

    // Only continue if we have an unbinned observation
    if (obs.eventtype() == "EventList") {

        // Extract region of interest from CTA observation
        const GCTAEventList* list = dynamic_cast<const GCTAEventList*>(obs.events());
        if (list == NULL) {
            std::string msg = "CTA Observation does not contain an event "
                              "list. Event list information is needed to "
                              "retrieve the Region of Interest for each "
                              "CTA observation.";
            throw GException::invalid_value(G_SET, msg);
        }
        const GCTARoi& roi = list->roi();

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
                    double  theta = pnt.dist(dir);

                    // Loop over all exposure cube energy bins
                    for (int iebin = 0; iebin < m_ebounds.size(); ++iebin){

                        // Get logE/TeV
                    	GEnergy srcEng = m_ebounds.elogmean(iebin);

                        // Loop over migra values
                        for (int imigra = 0; imigra < m_migras.size(); ++imigra) {

                            // Compute migra
                            double migra = m_migras[imigra];

                            // Set map index
                            int imap = offset(imigra, iebin);

                            // Set Edisp cube
                            m_cube(pixel, imap) = rsp->edisp(srcEng * migra, theta, 0.0, 0.0, 0.0, srcEng.log10TeV());

                        } // endfor: looped over migra bins

                    } // endfor: looped over energy bins

                } // endif: pixel was within RoI

            } // endfor: looped over all pixels

        } // endif: response was valid

    } // endif: observation was unbinned

    // Compile option: guarantee smooth Edisp
    #if defined(G_SMOOTH_EDISP)
    set_to_smooth();
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill Edisp cube from observation container
 *
 * @param[in] obs Observation container.
 * @param[in] log Pointer towards logger (default: NULL).
 ***************************************************************************/
void GCTACubeEdisp::fill(const GObservations& obs, GLog* log)
{
    // Clear Edisp cube
    clear_cube();

    // Initialise skymap for exposure weight accumulation
    GSkyMap exposure(m_cube);

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
            *log << " in point spread function cube computation." << std::endl;
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

                // Loop over all energy bins
                for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {

                    // Skip if pixel is not within observation energy boundaries
                    if (!obs_ebounds.contains(m_ebounds.emin(iebin)+g_energy_margin,
                                              m_ebounds.emax(iebin)-g_energy_margin)) {
                        continue;
                    }

                    // Get source energy
                    GEnergy srcEng = m_ebounds.elogmean(iebin);

                    // Get LogE/Tev
                    double logE = srcEng.log10TeV();

                    // Compute exposure weight
                    double weight = rsp->aeff(theta, 0.0, 0.0, 0.0, logE) *
                                    cta->livetime();

                    // Accumulate weights
                    exposure(pixel, iebin) += weight;

                    // Loop over delta values
                    for (int imigra = 0; imigra < m_migras.size(); ++imigra) {

                        // Compute delta in radians
                        double migra = m_migras[imigra];

                        // Set map index
                        int imap = offset(imigra, iebin);

                        // Add on Edisp cube
                        m_cube(pixel, imap) +=
                            rsp->edisp(srcEng * migra, theta, 0.0, 0.0, 0.0, logE) * weight;

                    } // endfor: looped over migra bins

                } // endfor: looped over energy bins

            } // endif: pixel was within RoI

        } // endfor: looped over all pixels

    } // endfor: looped over observations

    // Compute mean Edisp cube by dividing though the weights
    for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {
        for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {
            if (exposure(pixel, iebin) > 0.0) {
                double norm = 1.0 / exposure(pixel, iebin);
                for (int imigra = 0; imigra < m_migras.size(); ++imigra) {
                    int imap = offset(imigra, iebin);
                    m_cube(pixel, imap) *= norm;
                }
            }
            else {
                for (int imigra = 0; imigra < m_migras.size(); ++imigra) {
                    int imap = offset(imigra, iebin);
                    m_cube(pixel, imap) = 0.0;
                }
            }
        }
    }

    // Compile option: guarantee smooth Edisp
    #if defined(G_SMOOTH_EDISP)
    set_to_smooth();
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Edisp cube from FITS object
 *
 * @param[in] fits FITS object.
 *
 * Read the Edisp cube from a FITS object.
 ***************************************************************************/
void GCTACubeEdisp::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get HDUs
    const GFitsImage& hdu_psfcube = *fits.image("Primary");
    const GFitsTable& hdu_ebounds = *fits.table("EBOUNDS");
    const GFitsTable& hdu_migras  = *fits.table("MIGRAS");

    // Read cube
    m_cube.read(hdu_psfcube);

    // Read energy boundaries
    m_ebounds.read(hdu_ebounds);

    // Read migra nodes
    m_migras.read(hdu_migras);

    // Set migra node array for computation
    set_migra_axis();

    // Set energy node array
    set_eng_axis();

    // Compile option: guarantee smooth Edisp
    #if defined(G_SMOOTH_EDISP)
    set_to_smooth();
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA Edisp cube into FITS object.
 *
 * @param[in] fits FITS object.
 *
 * Write the CTA Edisp cube into a FITS object.
 ***************************************************************************/
void GCTACubeEdisp::write(GFits& fits) const
{
    // Write cube
    m_cube.write(fits);

    // Write energy boundaries
    m_ebounds.write(fits);

    // Write delta nodes
    m_migras.write(fits, "MIGRAS");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Edisp cube from FITS file
 *
 * @param[in] filename Performance table file name.
 *
 * Loads the Edisp cube from a FITS file into the object.
 ***************************************************************************/
void GCTACubeEdisp::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read Edisp cube
    read(fits);

    // Close Edisp file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Edisp cube into FITS file
 *
 * @param[in] filename Edisp cube FITS file name.
 * @param[in] clobber Overwrite existing file? (default: false)
 *
 * Save the Edisp cube into a FITS file.
 ***************************************************************************/
void GCTACubeEdisp::save(const GFilename& filename, const bool& clobber) const
{
    // Create FITS file
    GFits fits;

    // Write PSF cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return true energy interval that contains the energy dispersion.
 *
 * @param[in] logEobs Log10 of the observed event energy (TeV).
 * @param[in] theta Offset angle in camera system (rad).
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 *
 * Returns the band of true photon energies outside of which the energy
 * dispersion becomes negligible for a given observed energy @p logEobs and
 * offset angle @p theta. An energy in considered negligible if inferior to
 * 1e-6.
 ***************************************************************************/
GEbounds GCTACubeEdisp::ebounds_src(const GSkyDir& dir,
                                  const GEnergy& obsEng) const
{
    // Get index of obsEng inside energy boundaries
	int index = m_ebounds.index(obsEng);

	// Get map indices containing the migration information
	int binmin = m_migras.size() * index;
	int binmax = m_migras.size() * (index + 1);

	m_cube(dir, binmin)


	// Compute only if parameters changed
    if (!m_ebounds_src_computed || theta != m_last_theta_src) {

        // Set computation flag
        m_ebounds_src_computed = true;
        m_last_theta_src       = theta;

        // Compute ebounds_src
        compute_ebounds_src(theta, phi, zenith, azimuth);

    }

    // Search index only if logEobs has changed
    if (logEobs != m_last_logEobs) {

        // Store observed energy
        m_last_logEobs = logEobs;

        // Find right index with bisection
        int low  = 0;
        int high = m_ebounds_src.size() - 1;
        while ((high-low) > 1) {
            int  mid = (low+high) / 2;
            double e = m_edisp.axis_lo(m_inx_etrue, mid);
            if (logEobs < std::log10(e)) {
                high = mid;
            }
            else {
                low = mid;
            }
        }

        // Index found
        m_index_src = low;

    }

    // Return energy boundaries
    return m_ebounds_src[m_index_src];
}

/***********************************************************************//**
 * @brief Print Edisp cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing Edisp cube information.
 ***************************************************************************/
std::string GCTACubeEdisp::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTACubeEdisp ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);

        // Append energy intervals
        if (m_ebounds.size() > 0) {
            result.append("\n"+m_ebounds.print(chatter));
        }
        else {
            result.append("\n"+gammalib::parformat("Energy intervals") +
                          "not defined");
        }

        // Append number of migra bins
        result.append("\n"+gammalib::parformat("Number of migra bins") +
                      gammalib::str(m_migras.size()));

        // Append migra range
        result.append("\n"+gammalib::parformat("Migra range"));
        if (m_migras.size() > 0) {
            result.append(gammalib::str(m_migras[0]));
            result.append(" - ");
            result.append(gammalib::str(m_migras[m_migras.size()-1]));
        }
        else {
            result.append("not defined");
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
void GCTACubeEdisp::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_cube.clear();
    m_ebounds.clear();
    m_elogmeans.clear();
    m_migras.clear();
    m_migras_cache.clear();

    // Initialise cache
    m_inx1 = 0;
    m_inx2 = 0;
    m_inx3 = 0;
    m_inx4 = 0;
    m_wgt1 = 0.0;
    m_wgt2 = 0.0;
    m_wgt3 = 0.0;
    m_wgt4 = 0.0;
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Edisp cube.
 ***************************************************************************/
void GCTACubeEdisp::copy_members(const GCTACubeEdisp& cube)
{
    // Copy members
    m_filename          = cube.m_filename;
    m_cube              = cube.m_cube;
    m_ebounds           = cube.m_ebounds;
    m_elogmeans         = cube.m_elogmeans;
    m_migras            = cube.m_migras;
    m_migras_cache      = cube.m_migras_cache;

    // Copy cache
    m_inx1 = cube.m_inx1;
    m_inx2 = cube.m_inx2;
    m_inx3 = cube.m_inx3;
    m_inx4 = cube.m_inx4;
    m_wgt1 = cube.m_wgt1;
    m_wgt2 = cube.m_wgt2;
    m_wgt3 = cube.m_wgt3;
    m_wgt4 = cube.m_wgt4;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTACubeEdisp::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clear all pixels in the Edisp cube
 ***************************************************************************/
void GCTACubeEdisp::clear_cube(void)
{
    // Loop over all maps
    for (int imap = 0; imap < m_cube.nmaps(); ++imap) {

        // Loop over all pixels in sky map
        for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

            // Reset cube value to zero
            m_cube(pixel, imap) = 0.0;

        } // endfor: looped over all pixels

    } // endfor: looped over maps

    // Return
    return;
}

/***********************************************************************//**
 * @brief Update Edisp parameter cache
 *
 * @param[in] delta Angular separation between true and measured photon
 *            directions (radians).
 * @param[in] logE Log10 true photon energy (TeV). 
 *
 * This method updates the Edisp parameter cache.
 ***************************************************************************/
void GCTACubeEdisp::update(const double& migra, const double& logE) const
{
    // Set node array for delta interpolation
	m_migras_cache.set_value(migra);

    // Set node array for energy interpolation
    m_elogmeans.set_value(logE);
   
    // Set indices for bi-linear interpolation
    m_inx1 = offset(m_migras_cache.inx_left(),  m_elogmeans.inx_left());
    m_inx2 = offset(m_migras_cache.inx_left(),  m_elogmeans.inx_right());
    m_inx3 = offset(m_migras_cache.inx_right(), m_elogmeans.inx_left());
    m_inx4 = offset(m_migras_cache.inx_right(), m_elogmeans.inx_right());

    // Set weighting factors for bi-linear interpolation
    m_wgt1 = m_migras_cache.wgt_left()  * m_elogmeans.wgt_left();
    m_wgt2 = m_migras_cache.wgt_left()  * m_elogmeans.wgt_right();
    m_wgt3 = m_migras_cache.wgt_right() * m_elogmeans.wgt_left();
    m_wgt4 = m_migras_cache.wgt_right() * m_elogmeans.wgt_right();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for migra axis
 *
 * Set the migra axis nodes.
 *
 * @todo Check that none of the axis boundaries is non-positive.
 ***************************************************************************/
void GCTACubeEdisp::set_migra_axis(void)
{
    // Initialise computation members
    m_migras_cache.clear();

    // use a linear binning
	m_migras_cache.clear();
	for (int i = 0; i < m_migras.size(); ++i) {
		m_migras_cache.append(m_migras[i]);
	}

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
void GCTACubeEdisp::set_eng_axis(void)
{
    // Get number of bins
    int bins = m_ebounds.size();

    // Clear nodes
    m_elogmeans.clear();

    // Compute nodes
    for (int i = 0; i < bins; ++i) {
     
        // Append logE/TeV
        m_elogmeans.append(m_ebounds.elogmean(i).log10TeV());

    }  // endfor: looped over energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Pad the last migra bins with zero
 *
 * Zero padding of the last delta bins assures that the Edisp goes to zero
 * without any step at the last migra value.
 ***************************************************************************/
void GCTACubeEdisp::set_to_smooth(void)
{
    // Continue only if there are migra bins
    if (m_migras.size() > 2) {

        // Set first delta bin to value of second bin (capped Psf)
        for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {
            int isrc = offset(1, iebin);
            int idst = offset(0, iebin);
            for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {
                m_cube(pixel, idst) = m_cube(pixel, isrc);
            }
        }
    
        // Get index of last delta bin
        int imigra = m_migras.size()-1;

        // Pad mean Edisp with zeros in the last delta bin
        for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {
            int imap = offset(imigra, iebin);
            for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {
                m_cube(pixel, imap) = 0.0;
            }
        }

    } // endif: there were bins to pad

    // Return
    return;
}
