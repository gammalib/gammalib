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
#include "GFitsBinTable.hpp"
#include "GCTABackground3D.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MC                  "GCTABackground3D::mc(GEnergy&, GTime&, GRan&)"

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


/***********************************************************************//**
 * @brief Return background rate in units of events/s/MeV/sr
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] detx Tangential coord in nominal sys (rad).
 * @param[in] dety Tangential coord in nominal sys (rad).
 * @param[in] etrue Use true energy (true/false). Defaults to false.
 *
 * Returns the background rate in units of events/s/MeV/sr for a given energy
 * and detector coordinates. The method assures that the background rate
 * never becomes negative.
 *
 * The method supports true and reconstructed energies for logE. To access
 * the background rate as function of true energy, specify etrue=true
 * (this is the default). The obtained the background rate as function of
 * reconstructed energy, specify etrue=false.
 *
 * The method interpolates linearly in DETX and DETY, and logarithmically
 * in the energy direction.
 ***************************************************************************/
double GCTABackground3D::operator()(const double& logE, 
                                    const double& detx, 
                                    const double& dety,
                                    const bool&   etrue) const
{
    // Set parameter index
    int index = (etrue) ? 0 : 1;

    // Get background rate
    #if defined(G_LOG_INTERPOLATION)
    // Retrieve node arrays
    const GNodeArray& detx_nodes   = m_background.nodes(0);
    const GNodeArray& dety_nodes   = m_background.nodes(1);
    const GNodeArray& energy_nodes = m_background.nodes(2);

    // Set values for node arrays
    detx_nodes.set_value(detx);
    dety_nodes.set_value(dety);
    energy_nodes.set_value(logE);

    // Compute offsets of DETY in DETX-DETY plane
    int size1        = m_background.axis(0);
    int offset_left  = dety_nodes.inx_left()  * size1;
    int offset_right = dety_nodes.inx_right() * size1;

    // Set indices for bi-linear interpolation in DETX-DETY plane
    int inx_ll = detx_nodes.inx_left()  + offset_left;
    int inx_lr = detx_nodes.inx_left()  + offset_right;
    int inx_rl = detx_nodes.inx_right() + offset_left;
    int inx_rr = detx_nodes.inx_right() + offset_right;

    // Set weighting factors for bi-linear interpolation in DETX-DETY plane
    double wgt_ll = detx_nodes.wgt_left()  * dety_nodes.wgt_left();
    double wgt_lr = detx_nodes.wgt_left()  * dety_nodes.wgt_right();
    double wgt_rl = detx_nodes.wgt_right() * dety_nodes.wgt_left();
    double wgt_rr = detx_nodes.wgt_right() * dety_nodes.wgt_right();

    // Set indices for energy interpolation
    int inx_emin = energy_nodes.inx_left();
    int inx_emax = energy_nodes.inx_right();

    // Set weighting factors for energy interpolation
    double wgt_emin = energy_nodes.wgt_left();
    double wgt_emax = energy_nodes.wgt_right();

    // Compute offsets in energy dimension
    int npixels     = m_background.axis(0) * m_background.axis(1);
    int offset_emin = inx_emin * npixels;
    int offset_emax = inx_emax * npixels;

    // Bi-linear interpolation the rates in both energy layers
    double rate(0.0);
    double rate_emin = wgt_ll * m_background(index, inx_ll + offset_emin) +
                       wgt_lr * m_background(index, inx_lr + offset_emin) +
                       wgt_rl * m_background(index, inx_rl + offset_emin) +
                       wgt_rr * m_background(index, inx_rr + offset_emin);
    double rate_emax = wgt_ll * m_background(index, inx_ll + offset_emax) +
                       wgt_lr * m_background(index, inx_lr + offset_emax) +
                       wgt_rl * m_background(index, inx_rl + offset_emax) +
                       wgt_rr * m_background(index, inx_rr + offset_emax);

    // If both rates are positive then perform a logarithmic
    // interpolation in energy
    if (rate_emin > 0.0 && rate_emax > 0.0) {
        rate = std::pow(10.0, wgt_emin * std::log10(rate_emin) +
                              wgt_emax * std::log10(rate_emax));
    }

    // ... otherwise perform a linear interpolation
    else {
        rate = wgt_emin * rate_emin + wgt_emax * rate_emax;
    }
    #else
    double rate = m_background(index, detx, dety, logE);
    #endif

    // Make sure that background rate is not negative
    if (rate < 0.0) {
        rate = 0.0;
    }
    
    // Return background rate
    return rate;
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
 * @brief Save background table into FITS file
 *
 * @param[in] filename background table FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * Save the background table into a FITS file.
 * @todo Add necessary keywords.
 ***************************************************************************/
void GCTABackground3D::save(const std::string& filename, const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;
    GFitsBinTable hdu;
    hdu.extname("BACKGROUND");
    // Set additional keywords
    // hdu.card("HDUCLAS1","RESPONSE", "");
    //   hdu.card("HDUCLAS2","BGD", "");
    //   hdu.card("CDEC0001","CTA background", "");
    //  hdu.card("TDIM7","(14,14,15)", "");
    // hdu.card("TDIM8","(14,14,15)", "");
    //   hdu.card("TELESCOP", "CTA", "");
    //    hdu.card("INSTRUME", "CTA", "");
    // Write background table
    m_background.write(hdu);
   
    // Save FITS file
    fits.append(hdu);
    fits.saveto(filename, clobber);

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

    // Return
    return;
}

/***********************************************************************//**
 * @brief Write CTA background table into FITS object.
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GCTABackground3D::write(GFitsBinTable& hdu) const
{
    // Write cube
    m_background.write(hdu);

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

    // Continue only if there are cube pixels and maps
    if (m_mc_npix > 0 && m_mc_nmaps > 0) {

        // Determine the map that corresponds best to the specified energy.
        // This is not 100% clean, as ideally some map interpolation should
        // be done to the exact energy specified. However, as long as the map
        // does not change drastically with energy, taking the closest map
        // seems to be fine. And recall that we rebin in energy anyways if
        // the input binning is too coarse.
        double logE  = energy.log10TeV();
        int    index = -1;
        if (logE == m_mc_logE_min) {
            index = 0;
        } else if (logE == m_mc_logE_max) {
            index = m_mc_nmaps-1;
        }
        else {
            index = int((logE - m_mc_logE_min) / m_mc_logE_bin);
        }
        if (index < 0 || index >= m_mc_nmaps) {
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
        int offset = index * (m_mc_npix+1);
        int low    = offset;
        int high   = offset + m_mc_npix;
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
        int idetx = pixel % m_mc_nx;
        int idety = pixel / m_mc_nx;

        // Get randomized detx and dety
        double detx = m_mc_detx_min + (idetx + ran.uniform()) * m_mc_detx_bin;
        double dety = m_mc_dety_min + (idety + ran.uniform()) * m_mc_dety_bin;

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
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");
        result.append("\n"+gammalib::parformat("Maximum bin size for MC"));
        result.append(gammalib::str(m_mc_max_bin)+" deg");
        result.append("\n"+gammalib::parformat("Maximum logE step for MC"));
        result.append(gammalib::str(m_mc_max_logE)+"^10 TeV");

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
    m_mc_max_bin  = 0.05;  //!< Spatial binning not worse than 0.05 deg
    m_mc_max_logE = 0.02;  //!< Spectral binning not worse than 0.02^10 TeV

    // Initialise MC cache
    m_mc_cache.clear();
    m_mc_spectrum.clear();
    m_mc_nx       = 0;
    m_mc_ny       = 0;
    m_mc_npix     = 0;
    m_mc_nmaps    = 0;
    m_mc_detx_min = 0.0;
    m_mc_detx_max = 0.0;
    m_mc_detx_bin = 0.0;
    m_mc_dety_min = 0.0;
    m_mc_dety_max = 0.0;
    m_mc_dety_bin = 0.0;

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
    m_filename    = bgd.m_filename;
    m_background  = bgd.m_background;
    m_mc_max_bin  = bgd.m_mc_max_bin;
    m_mc_max_logE = bgd.m_mc_max_logE;

    // Copy MC cache
    m_mc_cache    = bgd.m_mc_cache;
    m_mc_spectrum = bgd.m_mc_spectrum;
    m_mc_nx       = bgd.m_mc_nx;
    m_mc_ny       = bgd.m_mc_ny;
    m_mc_npix     = bgd.m_mc_npix;
    m_mc_nmaps    = bgd.m_mc_nmaps;
    m_mc_detx_min = bgd.m_mc_detx_min;
    m_mc_detx_max = bgd.m_mc_detx_max;
    m_mc_detx_bin = bgd.m_mc_detx_bin;
    m_mc_dety_min = bgd.m_mc_dety_min;
    m_mc_dety_max = bgd.m_mc_dety_max;
    m_mc_dety_bin = bgd.m_mc_dety_bin;

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
 * @param[in] etrue Use true energy flag (true=yes, false=no)
 *
 * Initialises the cache for Monte Carlo sampling. The method uses the
 * members m_mc_max_bin and m_mc_max_logE to enforce an internal rebinning
 * in case that the provided background model information is coarsely
 * pixelised. This rebinning is needed to assure coherence between Monte
 * Carlo simulated data and the model.
 * 
 * @todo Verify assumption made about the solid angles of the response table
 *       elements.
 ***************************************************************************/
void GCTABackground3D::init_mc_cache(const bool& etrue) const
{
    // Initialise cache
    m_mc_cache.clear();
    m_mc_spectrum.clear();

    // Determine number of response cube pixels and maps
    m_mc_nx    = m_background.axis(0);
    m_mc_ny    = m_background.axis(1);
    m_mc_nmaps = m_background.axis(2);

    // Compute DETX boundaries and binsize in deg. If DETX binsize is too
    // coarse then request a finer sampling
    m_mc_detx_min = m_background.axis_lo(0,0);
    m_mc_detx_max = m_background.axis_hi(0,m_mc_nx-1);
    m_mc_detx_bin = (m_mc_detx_max - m_mc_detx_min) / m_mc_nx;
    if (m_mc_detx_bin > m_mc_max_bin) {
        m_mc_nx       = int((m_mc_detx_max - m_mc_detx_min) / m_mc_max_bin + 1.0);
        m_mc_detx_bin = (m_mc_detx_max - m_mc_detx_min) / m_mc_nx;
    }

    // Compute DETY boundaries and binsize in deg. If DETY binsize is too
    // coarse then request a finer sampling
    m_mc_dety_min = m_background.axis_lo(1,0);
    m_mc_dety_max = m_background.axis_hi(1,m_mc_ny-1);
    m_mc_dety_bin = (m_mc_dety_max - m_mc_dety_min) / m_mc_ny;
    if (m_mc_dety_bin > m_mc_max_bin) {
        m_mc_ny       = int((m_mc_dety_max - m_mc_dety_min) / m_mc_max_bin + 1.0);
        m_mc_dety_bin = (m_mc_dety_max - m_mc_dety_min) / m_mc_ny;
    }

    // Compute energy sampling. If energy sampling is too coarse then
    // request a finer sampling
    GEnergy emin(m_background.axis_lo(2,0), m_background.axis_lo_unit(2));
    GEnergy emax(m_background.axis_hi(2,m_mc_nmaps-1), m_background.axis_hi_unit(2));
    m_mc_logE_min = emin.log10TeV();
    m_mc_logE_max = emax.log10TeV();
    m_mc_logE_bin = (m_mc_logE_max - m_mc_logE_min) / m_mc_nmaps;
    if (m_mc_logE_bin > m_mc_max_logE) {
        m_mc_nmaps    = int((m_mc_logE_max - m_mc_logE_min) / m_mc_max_logE + 1.0);
        m_mc_logE_bin = (m_mc_logE_max - m_mc_logE_min) / m_mc_nmaps;
    }

    // Determine number of MC cube pixels
    m_mc_npix = m_mc_nx * m_mc_ny;

    // Determine solid angle of pixel (we assume here simply that we have
    // square pixels of identical solid angle; it needs to be checked whether
    // this is a valid assumption)
    double solidangle = m_mc_detx_bin * gammalib::deg2rad *
                        m_mc_dety_bin * gammalib::deg2rad;

    // Continue only if there are pixels and maps
    if (m_mc_npix > 0 && m_mc_nmaps > 0) {

        // Reserve space for all pixels in cache
        m_mc_cache.reserve((m_mc_npix+1)*m_mc_nmaps);

        // Loop over all maps
        for (int i = 0; i < m_mc_nmaps; ++i) {

            // Compute pixel offset
            int offset = i * (m_mc_npix+1);

            // Get logE
            //double logE = m_background.nodes(2)[i];
            double logE = m_mc_logE_min + (i+0.5) * m_mc_logE_bin;

            // Set first cache value to 0
            m_mc_cache.push_back(0.0);

            // Initialise cache with cumulative pixel fluxes and compute
            // total flux in response table for normalization. Negative
            // pixels are excluded from the cumulative map.
            double total_rate = 0.0;
            double xmin       = m_mc_detx_min * gammalib::deg2rad;
            double xbin       = m_mc_detx_bin * gammalib::deg2rad;
            double ymin       = m_mc_dety_min * gammalib::deg2rad;
            double ybin       = m_mc_dety_bin * gammalib::deg2rad;
            double detx       = xmin + 0.5 * xbin;
            for (int ix = 0; ix < m_mc_nx; ++ix, detx += xbin) {
                double dety = ymin + 0.5 * ybin;
                for (int iy = 0; iy < m_mc_ny; ++iy, dety += ybin) {
                    double rate = (*this)(logE, detx, dety, etrue) * solidangle;
                    if (rate > 0.0) {
                        total_rate += rate;
                    }
                    m_mc_cache.push_back(total_rate); // units: events/s/MeV
                }
            }

            // Normalize cumulative pixel fluxes so that the values in the
            // cache run from 0 to 1
            if (total_rate > 0.0) {
        		for (int k = 0, element = offset; k < m_mc_npix; ++k, ++element) {
        			m_mc_cache[element] /= total_rate;
        		}
        	}

            // Make sure that last pixel in the cache is >1
            m_mc_cache[m_mc_npix+offset] = 1.0001;

            // Set energy value (unit independent)
            GEnergy energy;
            energy.log10(logE, m_background.axis_lo_unit(2));
                
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
