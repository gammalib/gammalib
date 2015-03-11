/***************************************************************************
 *             GCTACubeBackground.cpp - CTA cube background class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Michael Mayer                                    *
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
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GCTACubeBackground.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                             "GCTACubeBackground::read(GFits&)"
#define G_MC                "GCTACubeBackground::mc(GEnergy&, GTime&, GRan&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_MC_INIT
//#define G_DEBUG_CACHE
#define G_LOG_INTERPOLATION   //!< Energy interpolate log10(background rate)

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
GCTACubeBackground::GCTACubeBackground(const std::string& filename)
{
    // Initialise class members
    init_members();

    // Load background from file
    load(filename);

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
 * @brief Return background rate in units of events/s/MeV/sr
 *
 * @param[in] dir Position in the sky where rate should be computed
 * @param[in] energy Energy at which the rate should be computed
 *
 * Returns the background rate in units of events/s/MeV/sr for a given energy
 * and sky coordinate. The method assures that the background rate
 * never becomes negative.
 *
 * The method interpolates logarithmically in the energy direction.
 *
 * @todo Check whether this should be a GSkyDir or rather a GInstDir object.
 ***************************************************************************/
double GCTACubeBackground::operator()(const GSkyDir& dir,
                                      const GEnergy& energy) const
{
    // Set indices and weighting factors for interpolation
    update(energy.log10TeV());

    // Perform interpolation
    double background = m_wgt_left  * m_cube(dir, m_inx_left) +
                        m_wgt_right * m_cube(dir, m_inx_right);

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
 * @brief Set background cube from skymap and energy boundaries
 *
 * @param[in] cube Sky map.
 * @param[in] ebounds Energy boundaries.
 *
 * Set this instance from sky map, energy boundaries
 ***************************************************************************/
void GCTACubeBackground::set(const GSkymap& cube, const GEbounds& ebounds)
{
    // Bring cube to clean state
    clear();

    // Copy skymap
    m_cube = cube;

    // Copy ebounds
    m_ebounds = ebounds;

    // Set energy node array
    set_eng_axis();

    // Initialise Monte Carlo cache
    init_mc_cache();

    // Return
    return;
}



/***********************************************************************//**
 * @brief Set Monte Carlo simulation cone
 *
 * @param[in] centre Simulation cone centre.
 * @param[in] radius Simulation cone radius (degrees).
 *
 * Sets the simulation cone centre and radius that defines the directions
 * that will be simulated using the mc() method.
 *
 * @todo Checks whether this is needed
 ***************************************************************************/
void GCTACubeBackground::set_mc_cone(const GSkyDir& centre,
                                     const double&  radius)
{
    // Initialise cache
    m_mc_cache.clear();
    m_mc_spectrum.clear();

    // Determine number of cube pixels and maps
    int npix  = m_cube.npix();
    int nmaps = m_cube.nmaps();

    // Continue only if there are pixels and maps
    if (npix > 0 && nmaps > 0) {

        // Reserve space for all pixels in cache
        m_mc_cache.reserve((npix+1)*nmaps);

        // Loop over all maps
        for (int i = 0; i < nmaps; ++i) {

            // Compute pixel offset
            int offset = i * (npix+1);

            // Set first cache value to 0
            m_mc_cache.push_back(0.0);

            // Initialise cache with cumulative pixel fluxes and compute
            // total flux in skymap for normalization. Negative pixels are
            // excluded from the cumulative map.
            double total_flux = 0.0;
            for (int k = 0; k < npix; ++k) {

                // Derive effective pixel radius from half opening angle
                // that corresponds to the pixel's solid angle. For security,
                // the radius is enhanced by 50%.
                double pixel_radius =
                       std::acos(1.0 - m_cube.solidangle(k)/gammalib::twopi) *
                       gammalib::rad2deg * 1.5;

                // Add up flux with simulation cone radius + effective pixel
                // radius. The effective pixel radius is added to make sure
                // that all pixels that overlap with the simulation cone are
                // taken into account. There is no problem of having even
                // pixels outside the simulation cone taken into account as
                // long as the mc() method has an explicit test of whether a
                // simulated event is contained in the simulation cone.
                double distance = centre.dist_deg(m_cube.pix2dir(k));
                if (distance <= radius+pixel_radius) {
                    double flux = m_cube(k,i) * m_cube.solidangle(k);
                    if (flux > 0.0) {
                        total_flux += flux;
                    }
                }

                // Push back flux
                m_mc_cache.push_back(total_flux); // units: ph/cm2/s/MeV
            }

            // Normalize cumulative pixel fluxes so that the values in the
            // cache run from 0 to 1
            if (total_flux > 0.0) {
                for (int k = 0; k < npix; ++k) {
                    m_mc_cache[k+offset] /= total_flux;
                }
            }

            // Make sure that last pixel in the cache is >1
            m_mc_cache[npix+offset] = 1.0001;

            // Store centre flux in node array
            if (m_elogmeans.size() == nmaps) {
                GEnergy energy;
                energy.log10MeV(m_elogmeans[i]);

                // Only append node if flux > 0
                if (total_flux > 0.0) {
                    m_mc_spectrum.append(energy, total_flux);
                }

            }

        } // endfor: looped over all maps

        // Dump cache values for debugging
        #if defined(G_DEBUG_CACHE)
        for (int i = 0; i < m_mc_cache.size(); ++i) {
            std::cout << "i=" << i;
            std::cout << " c=" << m_mc_cache[i] << std::endl;
        }
        #endif

    } // endif: there were cube pixels and maps

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 ***************************************************************************/
GSkyDir GCTACubeBackground::mc(const GEnergy& energy,
                               const GTime&   time,
                               GRan&          ran) const
{

    // Allocate sky direction
    GSkyDir dir;

    // Determine number of skymap pixels
    int npix = m_cube.npix();

    // Continue only if there are skymap pixels
    if (npix > 0) {

        // If no energy boundaries are defined, throw an exception
        if (m_ebounds.size() < 1) {
            std::string msg = "The energy boundaries of the maps in the cube"
                              " have not been defined. Maybe the map cube file"
                              " is missing the \"ENERGIES\" extension which"
                              " defines the energy of each map in the cube.\n"
                              "Please provide the energy information.";
            throw GException::invalid_value(G_MC, msg);
        }

        // Determine the map that corresponds best to the specified energy.
        // This is not 100% clean, as ideally some map interpolation should
        // be done to the exact energy specified. However, as long as the map
        // does not change drastically with energy, taking the closest map
        // seems to be fine.
        int i = m_ebounds.index(energy);
        if (i < 0) {
            if (energy <= m_ebounds.emin()) {
                i = 0;
            }
            else if (energy >= m_ebounds.emax()) {
                i = m_ebounds.size()-1;
            }
            else {
                std::string msg = "The specified energy "+energy.print()+" does"
                                  " not fall in any of the energy boundaries of"
                                  " the map cube.\n"
                                  "Please make sure that the map cube energies"
                                  " are properly defined.";
                throw GException::invalid_value(G_MC, msg);
            }
        }

        // Get uniform random number
        double u = ran.uniform();

        // Get pixel index according to random number. We use a bi-section
        // method to find the corresponding skymap pixel
        int offset = i * (npix+1);
        int low    = offset;
        int high   = offset + npix;
        while ((high - low) > 1) {
            int mid = (low+high) / 2;
            if (u < m_mc_cache[mid]) {
                high = mid;
            }
            else if (m_mc_cache[mid] <= u) {
                low = mid;
            }
        }

        // Convert sky map index to sky map pixel
        GSkyPixel pixel = m_cube.inx2pix(low-offset);

        // Randomize pixel
        pixel.x(pixel.x() + ran.uniform() - 0.5);
        pixel.y(pixel.y() + ran.uniform() - 0.5);

        // Get sky direction
        dir = m_cube.pix2dir(pixel);

    } // endif: there were pixels in sky map

    // Return instrument direction
    return dir;
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
    const GFitsTable& hdu_ebounds = *fits.table("EBOUNDS");

    // Read cube
    m_cube.read(hdu_bgdcube);

    // Read energy boundaries
    m_ebounds.read(hdu_ebounds);

    // Set energy node array
    set_eng_axis();

    // Initialise Monte Carlo cache
    init_mc_cache();

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
void GCTACubeBackground::load(const std::string& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read background cube
    read(fits);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;

}


/***********************************************************************//**
 * @brief Save background cube into FITS file
 *
 * @param[in] filename background cube FITS file name.
 * @param[in] clobber Overwrite existing file? (default: false)
 *
 * Save the background cube into a FITS file.
 ***************************************************************************/
void GCTACubeBackground::save(const std::string& filename,
                              const bool&        clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write background cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print background information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
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
         result.append("\n"+gammalib::parformat("Cube file")+m_filename);

         // Append detailed information only if a map cube exists
         if (m_cube.npix() > 0) {

             // NORMAL: Append sky map
             if (chatter >= NORMAL) {
                 result.append("\n"+m_cube.print(chatter));
             }

             // EXPLICIT: Append energy nodes
             if (chatter >= EXPLICIT && m_elogmeans.size() > 0) {
                 result.append("\n"+gammalib::parformat("Cube energy values"));
                 if (m_elogmeans.size() > 0) {
                     for (int i = 0; i < m_elogmeans.size(); ++i) {
                         result.append("\n"+gammalib::parformat("  Map "+gammalib::str(i+1)));
                         result.append(gammalib::str(std::pow(10.0, m_elogmeans[i])));
                         result.append(" MeV (log10E=");
                         result.append(gammalib::str(m_elogmeans[i]));
                         result.append(")");
                         if (m_ebounds.size() == m_elogmeans.size()) {
                             result.append(" [");
                             result.append(m_ebounds.emin(i).print());
                             result.append(", ");
                             result.append(m_ebounds.emax(i).print());
                             result.append("]");
                         }
                     }
                 }
                 else {
                     result.append("not specified");
                 }
             }

             // VERBOSE: Append MC cache
             if (chatter >= VERBOSE) {
                 result.append("\n"+gammalib::parformat("Map flux"));
                 if (m_mc_spectrum.nodes() > 0) {
                     for (int i = 0; i < m_mc_spectrum.nodes(); ++i) {
                         result.append("\n"+gammalib::parformat("  Map "+gammalib::str(i+1)));
                         result.append(gammalib::str(m_mc_spectrum.intensity(i)));
                     }
                 }
                 else {
                     result.append("not specified");
                 }
             }

         } // endif: map cube exists

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

    // Initialise MC cache
    m_mc_cache.clear();
    m_mc_spectrum.clear();
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

    // Copy MC cache
    m_mc_cache    = bgd.m_mc_cache;
    m_mc_spectrum = bgd.m_mc_spectrum;
    m_inx_left    = bgd.m_inx_left;
    m_inx_right   = bgd.m_inx_right;
    m_wgt_left    = bgd.m_wgt_left;
    m_wgt_right   = bgd.m_wgt_right;

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
 * @brief Initialise Monte Carlo cache
 ***************************************************************************/
void GCTACubeBackground::init_mc_cache(void)
{
    // Set centre and radius to all sky
    GSkyDir centre;
    double  radius = 360.0;

    // Compute cache
    set_mc_cone(centre, radius);

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
void GCTACubeBackground::set_eng_axis(void)
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
 * @brief Compute spatially averaged background raye for a given energy
 *
 * @param[in] logE Logarithm (base 10) of energy in TeV
 * @return Spatially averaged background rate in units of events/s/MeV/sr
 *
 * Spatially integrates the background cube at a given energy. This method
 * performs an interpolation between the energy maps. The integration is
 * limited to counting the bin contents
 ***************************************************************************/
double GCTACubeBackground::integral(const double& logE) const
{
    // Update interpolation cache
    update(logE);

    // Initialise result
    double result     = 0.0;
    double solidangle = 0.0;

    // Loop over all map pixels
    for (int i = 0; i < m_cube.npix(); ++i) {

        // Sum bin contents
        result += m_wgt_left  * m_cube(i, m_inx_left) +
                  m_wgt_right * m_cube(i, m_inx_right);

        // Sum solid angles of pixels
        solidangle += m_cube.solidangle(i);

    }

    // Divide integral by solidangle of map
    if (solidangle > 0.0) {
        result /= solidangle;
    }

    // Return result
    return result;
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
