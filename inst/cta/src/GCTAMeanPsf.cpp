/***************************************************************************
 *       GCTAMeanPsf.cpp - CTA mean point spread function cube class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @file GCTAMeanPsf.cpp
 * @brief CTA mean point spread function class cube implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAMeanPsf.hpp"
#include "GCTAObservation.hpp"
#include "GCTAExposure.hpp"
#include "GMath.hpp"

/* __ Method name definitions ____________________________________________ */

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
GCTAMeanPsf::GCTAMeanPsf(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Point spread function.
 ***************************************************************************/
GCTAMeanPsf::GCTAMeanPsf(const GCTAMeanPsf& cube)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Mean PSF cube constructor
 *
 * @param[in] obs     Observation container.
 * @param[in] wcs     World Coordinate System.
 * @param[in] coords  Coordinate System (CEL or GAL).
 * @param[in] x       X coordinate of sky map centre (deg).
 * @param[in] y       Y coordinate of sky map centre (deg).
 * @param[in] dx      Pixel size in x direction at centre (deg/pixel).
 * @param[in] dy      Pixel size in y direction at centre (deg/pixel).
 * @param[in] nx      Number of pixels in x direction.
 * @param[in] ny      Number of pixels in y direction.
 * @param[in] ebounds Energy boundaries.
 * @param[in] dmin    Minimum delta (deg).
 * @param[in] dmax    Maximum delta (deg.
 * @param[in] ndbins  Number of delta bins.
 *
 * Constructs a mean PSF cube by computing the mean PSF from all CTA
 * observations found in the observation container.
 ***************************************************************************/
GCTAMeanPsf::GCTAMeanPsf(const GObservations& obs,
                         const std::string&   wcs,
                         const std::string&   coords,
                         const double&        x,
                         const double&        y,
                         const double&        dx,
                         const double&        dy,
                         const int&           nx,
                         const int&           ny,
                         const GEbounds&      ebounds,
                         const double&        dmin,
                         const double&        dmax,
                         const int&           ndbins)
{
    // Initialise class members
    init_members();

    // Store energy boundaries
    m_ebounds = ebounds;

    // Set delta node array
    m_deltas.clear();
    for (int i = 0; i < nbins; ++i) {
        double binsize = (max*max - min*min)/nbins;
        double delta   = std::sqrt(binsize*0.5*i + min*min);
        m_deltas.append(delta);
    }

    // Compute number of sky maps
    int nmaps = m_ebounds.size() * m_deltas.size();
    
    // Create sky map
    m_cube = GSkymap(wcs, coords, x, y, dx, dy, nx, ny, nmaps);

    // Fill the PSF cube
    fill(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAMeanPsf::~GCTAMeanPsf(void)
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
 * @param[in] psf Mean PSF cube.
 * @return Mean PSF cube.
 ***************************************************************************/
GCTAMeanPsf& GCTAMeanPsf::operator= (const GCTAMeanPsf& cube)
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
void GCTAMeanPsf::clear(void)
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
 * @return Deep copy of mean PSF instance.
 ***************************************************************************/
GCTAMeanPsf* GCTAMeanPsf::clone(void) const
{
    return new GCTAMeanPsf(*this);
}


/***********************************************************************//**
 * @brief Set PSF cube from one CTA observation
 *
 * @param[in] obs CTA observation.
 ***************************************************************************/
void GCTAMeanPsf::set(const GCTAObservation& obs)
{
    // Clear PSF cube
    clear_cube();

    // Get references on CTA response and pointing direction
    const GCTAResponse& rsp = obs.response();
    const GSkyDir&      pnt = obs.pointing().dir();

    // Loop over all pixels in sky map
    for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

        // Compute theta angle with respect to pointing direction
        // in radians
        GSkyDir dir     = m_cube.inx2dir(pixel);
        double  theta   = pnt.dist(dir);
    
        // Loop over all exposure cube energy bins
        for (int iebin = 0; iebin < m_ebounds.size(); ++iebin){

            // Get logE/TeV
            double logE = m_ebounds.emean(iebin).log10TeV();

            // Loop over delta values
            for (int idelta = 0; idelta < m_deltas.size(); ++idelta) {

                // Compute delta in radians
                double delta = m_deltas[idelta] * gammalib::deg2rad;

                // Set map index
                int imap = offset(idelta, iebin);
                
                // Set PSF cube
                m_cube(pixel, imap) = rsp.psf(delta, theta, 0.0, 0.0, 0.0, logE);

            } // endfor: looped over delta bins
        } // endfor: looped over energy bins
    } // endfor: looped over all pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill PSF cube from observation container
 *
 * @param[in] obs Observation container.
 ***************************************************************************/
void GCTAMeanPsf::fill(const GObservations& obs)
{
    // Clear PSF cube
    clear_cube();

    // Initialise skymap for exposure weight accumulation
    GSkymap exposure(m_cube);

    // Loop over all observations in container
    for (int i = 0; i < obs.size(); ++i) {

        // Get observation and continue only if it is a CTA observation
        const GCTAObservation *cta = dynamic_cast<const GCTAObservation*>(m_obs[i]);
        if (cta != NULL) {

            // Get references on CTA response and pointing direction
            const GCTAResponse& rsp = cta.response();
            const GSkyDir&      pnt = cta.pointing().dir();

            // Loop over all pixels in sky map
            for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

                // Compute theta angle with respect to pointing direction
                // in radians
                GSkyDir dir   = m_cube.inx2dir(pixel);
                double  theta = pnt.dist(dir);
    
                // Loop over all energy bins
                for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {

                    // Get logE/TeV
                    double logE = m_ebounds.emean(iebin).log10TeV();

                    // Compute exposure weight
                    double weight = rsp.aeff(theta, 0.0, 0.0, 0.0, logE) *
                                    cta->livetime();

                    // Accumulate weights
                    exposure(pixel, iebin) += weight;

                    // Loop over delta values
                    for (int idelta = 0; idelta < m_deltas.size(); ++idelta) {

                        // Compute delta in radians
                        double delta = m_deltas[idelta] * gammalib::deg2rad;

                        // Set map index
                        int imap = offset(idelta, iebin);
                
                        // Add on PSF cube
                        m_cube(pixel, imap) +=
                           rsp.psf(delta, theta, 0.0, 0.0, 0.0, logE) * weight;

                    } // endfor: looped over delta bins

                } // endfor: looped over energy bins

            } // endfor: looped over all pixels

        } // endif: observation was a CTA observation

    } // endfor: looped over observations

    // Compute mean PSF cube by dividing though the weights
    for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {
        for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {
            if (exposure(pixel, iebin) > 0.0) {
                double norm = 1.0 / exposure(pixel, iebin);
                for (int idelta = 0; idelta < m_deltas.size(); ++idelta) {
                    int imap = offset(idelta, iebin);
                    m_cube(pixel, imap) *= norm;
                }
            }
            else {
                for (int idelta = 0; idelta < m_deltas.size(); ++idelta) {
                    m_cube(pixel, imap) = 0.0;
                }
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA PSF cube into FITS file.
 *
 * @param[in] fits FITS file.
 *
 * @todo Write also delta binning information
 ***************************************************************************/
void GCTAMeanPsf::write(GFits& fits) const
{
    // Write cube
    m_cube.write(fits);

    // Write energy boundaries
    m_ebounds.write(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load PSF cube from FITS file
 *
 * @param[in] filename Performance table file name.
 *
 * Loads the PSF cube from a FITS file into the object.
 *
 * @todo Implement method
 ***************************************************************************/
void GCTAMeanPsf::load(const std::string& filename)
{
    return;
}


/***********************************************************************//**
 * @brief Save PSF cube into FITS file
 *
 * @param[in] filename PSF cube FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * Save the PSF cube into a FITS file.
 *
 * @todo Implement method
 ***************************************************************************/
void GCTAMeanPsf::save(const std::string& filename, const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write PSF cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);
    return;
}


/***********************************************************************//**
 * @brief Print PSF cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing PSF cube information.
 *
 * @todo Add content
 ***************************************************************************/
std::string GCTAMeanPsf::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAMeanPsf ===");

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
void GCTAMeanPsf::init_members(void)
{
    // Initialise members
    m_cube.clear();
    m_ebounds.clear();
    m_deltas.clear();
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube PSF cube.
 ***************************************************************************/
void GCTAMeanPsf::copy_members(const GCTAMeanPsf& cube)
{
    // Copy members
    m_cube    = cube.m_cube;
    m_ebounds = cube.m_ebounds;
    m_deltas  = cube.m_deltas;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAMeanPsf::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clear all pixels in the PSF cube
 ***************************************************************************/
void GCTAMeanPsf::clear_cube(void)
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
