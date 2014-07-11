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
 * @param[in] wcs     World Coordinate System.
 * @param[in] coords  Coordinate System (CEL or GAL).
 * @param[in] x       X coordinate of sky map centre (deg).
 * @param[in] y       Y coordinate of sky map centre (deg).
 * @param[in] dx      Pixel size in x direction at centre (deg/pixel).
 * @param[in] dy      Pixel size in y direction at centre (deg/pixel).
 * @param[in] nx      Number of pixels in x direction.
 * @param[in] ny      Number of pixels in y direction.
 * @param[in] ebounds Energy boundaries.
 * @param[in] dmax    Maximum delta (deg).
 * @param[in] ndbins  Number of delta bins.
 *
 * Constructs a mean PSF cube by computing the mean PSF from all CTA
 * observations found in the observation container.
 *
 * @todo Think about the way how the delta node array is computed. The PSF
 * should be more finely sampled when being close to the peak. For the
 * moment a simple linear sampling is used.
 ***************************************************************************/
GCTAMeanPsf::GCTAMeanPsf(const std::string&   wcs,
                         const std::string&   coords,
                         const double&        x,
                         const double&        y,
                         const double&        dx,
                         const double&        dy,
                         const int&           nx,
                         const int&           ny,
                         const GEbounds&      ebounds,
                         const double&        dmax,
                         const int&           ndbins)
{
    // Initialise class members
    init_members();

    // Store energy boundaries
    m_ebounds = ebounds;

    // Set energy node array
    set_eng_axis();

    // Set delta node array
    m_deltas.clear();
    for (int i = 0; i < ndbins; ++i) {
        double binsize = dmax / double(ndbins);
        double delta   = binsize * (double(i)+0.5);
        m_deltas.append(delta);
    }

    // Compute number of sky maps
    int nmaps = m_ebounds.size() * m_deltas.size();
    
    // Create sky map
    m_cube = GSkymap(wcs, coords, x, y, dx, dy, nx, ny, nmaps);

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
 * @param[in] cube Mean PSF cube.
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


/***********************************************************************//**
 * @brief Return point spread function (in units of sr^-1)
 *
 * @param[in] dir Coordinate of the true photon position.
 * @param[in] delta Angular separation between true and measured photon
 *            directions (rad).
 * @param[in] energy Energy of the true photon.
 * @return point spread function (in units of sr^-1)
 *
 * Returns the point spread function for a given angular separation in units
 * of sr^-1 for a given energy and coordinate.
 ***************************************************************************/
double GCTAMeanPsf::operator()(const GSkyDir& dir, 
                               const double&  delta,
                               const GEnergy& energy) const
{
    // Update indices and weighting factors for interpolation
    update(delta, energy.log10TeV());

    // Perform bi-linear interpolation
    double psf = m_wgt1 * m_cube(dir, m_inx1) +
                 m_wgt2 * m_cube(dir, m_inx2) +
                 m_wgt3 * m_cube(dir, m_inx3) +
                 m_wgt4 * m_cube(dir, m_inx4);

    // Return PSF
    return psf;
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
            double logE = m_ebounds.elogmean(iebin).log10TeV();

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
        const GCTAObservation *cta = dynamic_cast<const GCTAObservation*>(obs[i]);
        if (cta != NULL) {

            // Get references on CTA response and pointing direction
            const GCTAResponse& rsp = cta->response();
            const GSkyDir&      pnt = cta->pointing().dir();

            // Loop over all pixels in sky map
            for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

                // Compute theta angle with respect to pointing direction
                // in radians
                GSkyDir dir   = m_cube.inx2dir(pixel);
                double  theta = pnt.dist(dir);
    
                // Loop over all energy bins
                for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {

                    // Get logE/TeV
                    double logE = m_ebounds.elogmean(iebin).log10TeV();

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
                    int imap = offset(idelta, iebin);
                    m_cube(pixel, imap) = 0.0;
                }
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read PSF cube from FITS object
 *
 * @param[in] fits FITS object.
 *
 * Read the PSF cube from a FITS object.
 ***************************************************************************/
void GCTAMeanPsf::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get HDUs
    const GFitsImage& hdu_psfcube  = *fits.image("Primary");
    const GFitsTable& hdu_ebounds = *fits.table("EBOUNDS");
    const GFitsTable& hdu_deltas  = *fits.table("DELTAS");

    // Read cube
    m_cube.read(hdu_psfcube);

    // Read energy boundaries
    m_ebounds.read(hdu_ebounds);

    // Read delta nodes
    m_deltas.read(hdu_deltas);

    // Set energy node array
    set_eng_axis();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA PSF cube into FITS object.
 *
 * @param[in] fits FITS object.
 *
 * Write the CTA PSF cube into a FITS object.
 ***************************************************************************/
void GCTAMeanPsf::write(GFits& fits) const
{
    // Write cube
    m_cube.write(fits);

    // Write energy boundaries
    m_ebounds.write(fits);

    // Write delta nodes
    m_deltas.write(fits, "DELTAS");

    // Set the nodes unit to "deg"
    (*fits.table("DELTAS"))["Value"]->unit("deg");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load PSF cube from FITS file
 *
 * @param[in] filename Performance table file name.
 *
 * Loads the PSF cube from a FITS file into the object.
 ***************************************************************************/
void GCTAMeanPsf::load(const std::string& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read PSF cube
    read(fits);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save PSF cube into FITS file
 *
 * @param[in] filename PSF cube FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * Save the PSF cube into a FITS file.
 ***************************************************************************/
void GCTAMeanPsf::save(const std::string& filename, const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write PSF cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Return
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
    m_elogmeans.clear();
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
    m_cube      = cube.m_cube;
    m_ebounds   = cube.m_ebounds;
    m_elogmeans = cube.m_elogmeans;
    m_deltas    = cube.m_deltas;

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
void GCTAMeanPsf::set_eng_axis(void)
{
    // Get number of bins
    int bins = m_ebounds.size();

    // Clear nodes
    m_elogmeans.clear();

    // Compute nodes
    for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {
     
        // Append logE/TeV
        m_elogmeans.append(m_ebounds.elogmean(iebin).log10TeV());

    }  // endfor: looped over energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update PSF parameter cache
 *
 * @param[in] delta Angular separation between true and measured photon
 *            directions (rad).
 * @param[in] logE Log10 true photon energy (TeV). 
 *
 * This method updates the PSF parameter cache.
 ***************************************************************************/
void GCTAMeanPsf::update(const double& delta, const double& logE) const
{
    // Set node array interpolation values
    m_deltas.set_value(delta*gammalib::rad2deg);
    m_elogmeans.set_value(logE);
   
    // Set indices for bi-linear interpolation
    m_inx1 = offset(m_deltas.inx_left(),  m_elogmeans.inx_left());
    m_inx2 = offset(m_deltas.inx_left(),  m_elogmeans.inx_right());
    m_inx3 = offset(m_deltas.inx_right(), m_elogmeans.inx_left());
    m_inx4 = offset(m_deltas.inx_right(), m_elogmeans.inx_right());

    // Set weighting factors for bi-linear interpolation
    m_wgt1 = m_deltas.wgt_left()  * m_elogmeans.wgt_left();
    m_wgt2 = m_deltas.wgt_left()  * m_elogmeans.wgt_right();
    m_wgt3 = m_deltas.wgt_right() * m_elogmeans.wgt_left();
    m_wgt4 = m_deltas.wgt_right() * m_elogmeans.wgt_right();

    // Return
    return;
}
