/***************************************************************************
 *                GCTAExposure.cpp - CTA exposure cube class               *
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
 * @file GCTAExposure.cpp
 * @brief CTA exposure cube class implementation
 * @author Chia-Chun Lu
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAExposure.hpp"
#include "GCTAObservation.hpp"
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
GCTAExposure::GCTAExposure(void)
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
GCTAExposure::GCTAExposure(const GCTAExposure& cube)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(cube);

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
GCTAExposure::GCTAExposure(const std::string&   wcs,
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
    m_cube = GSkymap(wcs, coords, x, y, dx, dy, nx, ny, m_ebounds.size());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAExposure::~GCTAExposure(void)
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
 * @param[in] cube exposure cube
 * @return Exposure cube.
 ***************************************************************************/
GCTAExposure& GCTAExposure::operator= (const GCTAExposure& cube)
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
double GCTAExposure::operator()(const GSkyDir& dir, const GEnergy& energy) const
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
void GCTAExposure::clear(void)
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
GCTAExposure* GCTAExposure::clone(void) const
{
    return new GCTAExposure(*this);
}


/***********************************************************************//**
 * @brief Set exposure cube from one CTA observation
 *
 * @param[in] obs CTA observation.
 *
 * Set the exposure cube from a single CTA observations. The cube pixel
 * values are computed as product of the effective area and the livetime.
 ***************************************************************************/
void GCTAExposure::set(const GCTAObservation& obs)
{
    // Clear exposure cube
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

            // Set exposure cube (effective area * lifetime)
            m_cube(pixel, iebin) = rsp.aeff(theta, 0.0, 0.0, 0.0, logE) *
                                   obs.livetime();

        } // endfor: looped over energy bins

    } // endfor: looped over all pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill exposure cube from observation container
 *
 * @param[in] obs Observation container.
 *
 * Set the exposure cube by summing the exposure for all CTA observations in
 * an observation container. The cube pixel values are computed as the sum
 * over the products of the effective area and the livetime.
 ***************************************************************************/
void GCTAExposure::fill(const GObservations& obs)
{
    // Clear exposure cube
    clear_cube();

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
                GSkyDir dir     = m_cube.inx2dir(pixel);
                double  theta   = pnt.dist(dir);
    
                // Loop over all exposure cube energy bins
                for (int iebin = 0; iebin < m_ebounds.size(); ++iebin){

                    // Get logE/TeV
                    double logE = m_ebounds.elogmean(iebin).log10TeV();

                    // Add to exposure cube (effective area * lifetime)
                    m_cube(pixel, iebin) += rsp.aeff(theta, 0.0, 0.0, 0.0, logE) *
                                            cta->livetime();

                } // endfor: looped over energy bins

            } // endfor: looped over all pixels
    
        } // endif: observation was a CTA observation

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
void GCTAExposure::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get HDUs
    const GFitsImage& hdu_expcube = *fits.image("Primary");
    const GFitsTable& hdu_ebounds = *fits.table("EBOUNDS");

    // Read cube
    m_cube.read(hdu_expcube);

    // Read energy boundaries
    m_ebounds.read(hdu_ebounds);

    // Set energy node array
    set_eng_axis();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA exposure cube into FITS object.
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GCTAExposure::write(GFits& fits) const
{
    // Write cube
    m_cube.write(fits);

    // Write energy boundaries
    m_ebounds.write(fits);

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
void GCTAExposure::load(const std::string& filename)
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
 * @brief Save exposure cube into FITS file
 *
 * @param[in] filename Exposure cube FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * Save the exposure cube into a FITS file.
 *
 * @todo Implement method
 ***************************************************************************/
void GCTAExposure::save(const std::string& filename, const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write exposure cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

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
std::string GCTAExposure::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAExposure ===");

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
void GCTAExposure::init_members(void)
{
    // Initialise members
    m_cube.clear();
    m_ebounds.clear();
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Exposure cube
 ***************************************************************************/
void GCTAExposure::copy_members(const GCTAExposure& cube)
{
    // Initialise members
    m_cube    = cube.m_cube;
    m_ebounds = cube.m_ebounds;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAExposure::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clear all pixels in the exposure cube
 ***************************************************************************/
void GCTAExposure::clear_cube(void)
{
    // Loop over all exposure cube energy bins
    for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {

        // Loop over all pixels in sky map
        for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

            // Reset cube value to zero
            m_cube(pixel, iebin) = 0.0;

        } // endfor: looped over all pixels

    } // endfor: looped over energy bins

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
void GCTAExposure::update(const double& logE) const
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
void GCTAExposure::set_eng_axis(void)
{
    // Get number of bins
    int bins = m_ebounds.size();

    // Clear node array
    m_elogmeans.clear();

    // Compute nodes
    for (int iebin = 0; iebin < m_ebounds.size(); ++iebin) {
     
        // Get logE/TeV
        m_elogmeans.append(m_ebounds.elogmean(iebin).log10TeV()); 

    }  // endfor: looped over energy bins

    // Return
    return;
}
