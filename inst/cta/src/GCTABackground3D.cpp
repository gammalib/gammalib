/***************************************************************************
 *              GCTABackground3D.cpp - CTA 3D background class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTABackground3D.cpp
 * @brief CTA 3D background class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GCTABackground3D.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MC                  "GCTABackground3D::mc(GEnergy&, GTime&, GRan&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

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
GCTABackground3D::GCTABackground3D(void) : GCTABackground()
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
GCTABackground3D::GCTABackground3D(const std::string& filename) :
                  GCTABackground()
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
GCTABackground3D::GCTABackground3D(const GCTABackground3D& bgd) :
                  GCTABackground(bgd)
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
GCTABackground3D::~GCTABackground3D(void)
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
 * @param[in] bgd Background.
 * @return Background.
 ***************************************************************************/
GCTABackground3D& GCTABackground3D::operator=(const GCTABackground3D& bgd)
{
    // Execute only if object is not identical
    if (this != &bgd) {

        // Copy base class members
        this->GCTABackground::operator=(bgd);

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
void GCTABackground3D::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTABackground::free_members();

    // Initialise members
    this->GCTABackground::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone background
 *
 * @return Pointer to deep copy of background.
 ***************************************************************************/
GCTABackground3D* GCTABackground3D::clone(void) const
{
    return new GCTABackground3D(*this);
}


/***********************************************************************//**
 * @brief Load background from FITS file
 *
 * @param[in] filename FITS file.
 *
 * This method loads the background information from a FITS file.
 ***************************************************************************/
void GCTABackground3D::load(const std::string& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read background from file
    read(fits);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read background from FITS file
 *
 * @param[in] fits FITS file pointer.
 *
 * Reads the background form the FITS file extension "BACKGROUND". The data
 * are stored in m_background which is of type GCTAResponseTable. The DETX
 * and DETY axes will be set to radians, the energy axis will be set to
 * log10.
 ***************************************************************************/
void GCTABackground3D::read(const GFits& fits)
{
    // Clear response table
    m_background.clear();

    // Get background table
    const GFitsTable& table = *fits.table("BACKGROUND");

    // Read background table
    m_background.read(table);

    // Set DETX and DETY axis to radians
    m_background.axis_radians(0);
    m_background.axis_radians(1);

    // Set energy axis to logarithmic scale
    m_background.axis_log10(2);

    // Convert background
    //m_background.scale(0, 1.0e4);
    //m_background.scale(1, 1.0e4);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return CTA instrument direction.
 ***************************************************************************/
GCTAInstDir GCTABackground3D::mc(const GEnergy& energy,
                                 const GTime&   time,
                                 GRan&          ran) const
{
    // Initialise Monte Carlo Cache
    if (m_mc_cache.empty()) {
        init_mc_cache();
    }

    // Allocate instrument direction
    GCTAInstDir dir;

    // Determine number of cube pixels and maps
    int n_detx = m_background.axis(0);
    int n_dety = m_background.axis(1);
    int npix   = n_detx * n_dety;
    int nmaps  = m_background.axis(2);

    // Continue only if there are cube pixels and maps
    if (npix > 0 && nmaps > 0) {

        // Determine the map that corresponds best to the specified energy.
        // This is not 100% clean, as ideally some map interpolation should
        // be done to the exact energy specified. However, as long as the map
        // does not change drastically with energy, taking the closest map
        // seems to be fine.
        int index = -1;
        for (int i = 0; i < nmaps; ++i) {
            if ((energy(m_background.axis_lo_unit(2)) >= m_background.axis_lo(2,i)) &&
                (energy(m_background.axis_hi_unit(2)) <= m_background.axis_hi(2,i))) {
                index = i;
                break;
            }
        }
        if (index == -1) {
            std::string msg = "The specified energy "+energy.print()+" does"
                              " not fall in any of the energy boundaries of"
                              " the background response cube.\n"
                              "Please select an energy range that is comprised"
                              " in the range ["+
                              gammalib::str(m_background.axis_lo(2,0))+
                              " "+m_background.axis_lo_unit(2)+" - "+
                              gammalib::str(m_background.axis_hi(2,m_background.axis(2)-1))+
                              " "+m_background.axis_hi_unit(2)+"].";
            throw GException::invalid_value(G_MC, msg);
        }

        // Get uniform random number
        double u = ran.uniform();

        // Get pixel index according to random number. We use a bi-section
        // method to find the corresponding response cube pixel
        int offset = index * (npix+1);
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

        // Get pixel indices in x and y
        int pixel = low - offset;
        int idetx = pixel % n_detx;
        int idety = pixel / n_detx;
//std::cout << pixel << " " << idetx << " " << idety << std::endl;

        // Get minimum and binsize for the pixel
        double xbin_min = m_background.axis_lo(0,idetx);
        double xbin_bin = m_background.axis_hi(0,idetx) - xbin_min;
        double ybin_min = m_background.axis_lo(1,idety);
        double ybin_bin = m_background.axis_hi(1,idety) - ybin_min;

        // Randomize detx and dety
        double detx = xbin_min + xbin_bin * ran.uniform();
        double dety = ybin_min + ybin_bin * ran.uniform();

        // Set instrument direction (in radians)
        dir.detx(detx * gammalib::deg2rad);
        dir.dety(dety * gammalib::deg2rad);

    } // endif: there were cube pixels and maps

    // Return instrument direction
    return dir;
}


/***********************************************************************//**
 * @brief Print background information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing background information.
 ***************************************************************************/
std::string GCTABackground3D::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute DETX boundaries in deg
        double detx_min = m_background.axis_lo(0,0);
        double detx_max = m_background.axis_hi(0,m_background.axis(0)-1);

        // Compute DETY boundaries in deg
        double dety_min = m_background.axis_lo(1,0);
        double dety_max = m_background.axis_hi(1,m_background.axis(1)-1);

        // Compute energy boundaries in TeV
        double emin = m_background.axis_lo(2,0);
        double emax = m_background.axis_hi(2,m_background.axis(2)-1);


        // Append header
        result.append("=== GCTABackground3D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of DETX bins") +
                      gammalib::str(m_background.axis(0)));
        result.append("\n"+gammalib::parformat("Number of DETY bins") +
                      gammalib::str(m_background.axis(1)));
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_background.axis(2)));
        result.append("\n"+gammalib::parformat("DETX range"));
        result.append(gammalib::str(detx_min)+" - "+gammalib::str(detx_max)+" deg");
        result.append("\n"+gammalib::parformat("DETX range"));
        result.append(gammalib::str(dety_min)+" - "+gammalib::str(dety_max)+" deg");
        result.append("\n"+gammalib::parformat("Log10(Energy) range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");

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
void GCTABackground3D::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_background.clear();

    // Initialise MC cache
    m_mc_cache.clear();
    m_mc_spectrum.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bgd Background.
 ***************************************************************************/
void GCTABackground3D::copy_members(const GCTABackground3D& bgd)
{
    // Copy members
    m_filename   = bgd.m_filename;
    m_background = bgd.m_background;

    // Copy MC cache
    m_mc_cache    = bgd.m_mc_cache;
    m_mc_spectrum = bgd.m_mc_spectrum;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTABackground3D::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise Monte Carlo cache
 *
 * @todo Verify assumption made about the solid angles of the response table
 *       elements.
 * @todo Add optional sampling on a finer spatial grid.
 ***************************************************************************/
void GCTABackground3D::init_mc_cache(void) const
{
    // Initialise cache
    m_mc_cache.clear();
    m_mc_spectrum.clear();

    // Determine number of cube pixels and maps
    int npix  = m_background.axis(0) * m_background.axis(1);
    int nmaps = m_background.axis(2);

    // Compute DETX boundaries and width in deg
    double detx_min = m_background.axis_lo(0,0);
    double detx_max = m_background.axis_hi(0,m_background.axis(0)-1);
    double detx_bin = (detx_max - detx_min) / m_background.axis(0);

    // Compute DETY boundaries and width in deg
    double dety_min = m_background.axis_lo(1,0);
    double dety_max = m_background.axis_hi(1,m_background.axis(1)-1);
    double dety_bin = (dety_max - dety_min) / m_background.axis(1);

    // Determine solid angle of pixel (we assume here simply that we have
    // square pixels of identical solid angle; it needs to be checked whether
    // this is a valid assumption)
    double solidangle = detx_bin * gammalib::deg2rad *
                        dety_bin * gammalib::deg2rad;

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
            // total flux in response table for normalization. Negative
            // pixels are excluded from the cumulative map.
            double total_rate = 0.0;
        	for (int k = 0, element = i*npix; k < npix; ++k, ++element) {

                // Determine background rate
                double rate = m_background(element) * solidangle;
                if (rate > 0.0) {
                    total_rate += rate;
                }

                // Push back flux
                m_mc_cache.push_back(total_rate); // units: events/s/MeV
        	}

            // Normalize cumulative pixel fluxes so that the values in the
            // cache run from 0 to 1
            if (total_rate > 0.0) {
        		for (int k = 0, element = offset; k < npix; ++k, ++element) {
        			m_mc_cache[element] /= total_rate;
        		}
        	}

            // Make sure that last pixel in the cache is >1
            m_mc_cache[npix+offset] = 1.0001;

            // Set energy value (unit independent)
            GEnergy energy;
            energy.log10(m_background.nodes(2)[i], m_background.axis_lo_unit(2));
                
            // Only append node if rate > 0
            if (total_rate > 0.0) {
                m_mc_spectrum.append(energy, total_rate);
            }

            // Dump spectrum for debugging
            #if defined(G_DEBUG_MC_INIT)
            std::cout << "Energy=" << energy;
            std::cout << " Rate=" << total_rate;
            std::cout << " events/s/MeV" << std::endl;
            #endif

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
