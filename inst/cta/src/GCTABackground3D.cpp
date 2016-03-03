/***************************************************************************
 *              GCTABackground3D.cpp - CTA 3D background class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2016 by Juergen Knoedlseder                         *
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
#include "GException.hpp"
#include "GFilename.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GCTABackground3D.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                          "GCTABackground3D::read(GFitsTable&)"

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
 *
 * Constructs empty background.
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
 * Constructs background from a FITS file.
 ***************************************************************************/
GCTABackground3D::GCTABackground3D(const GFilename& filename) :
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
 *
 * Constructs background by copying from another background.
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
 *
 * Destructs background.
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
 *
 * Assigns background.
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
 * @param[in] detx Tangential coord in nominal sys (radians).
 * @param[in] dety Tangential coord in nominal sys (radians).
 *
 * Returns the background rate in units of events/s/MeV/sr for a given energy
 * and detector coordinates. The method assures that the background rate
 * never becomes negative.
 *
 * The method interpolates linearly in DETX and DETY, and logarithmically
 * in energy.
 *
 * If the background model is not valid the operator returns 0.
 ***************************************************************************/
double GCTABackground3D::operator()(const double& logE, 
                                    const double& detx, 
                                    const double& dety) const
{
    // Initialise rate
    double rate(0.0);

    // Continue only if background model is valid
    if (is_valid()) {

        // Retrieve references to node arrays
        const GNodeArray& detx_nodes   = m_background.axis_nodes(m_inx_detx);
        const GNodeArray& dety_nodes   = m_background.axis_nodes(m_inx_dety);
        const GNodeArray& energy_nodes = m_background.axis_nodes(m_inx_energy);

        // Set values for node arrays
        detx_nodes.set_value(detx);
        dety_nodes.set_value(dety);
        energy_nodes.set_value(logE);

        // Compute offsets of DETY in DETX-DETY plane
        int size1        = m_background.axis_bins(m_inx_detx);
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
        int npixels     = m_background.axis_bins(m_inx_detx) *
                          m_background.axis_bins(m_inx_dety);
        int offset_emin = inx_emin * npixels;
        int offset_emax = inx_emax * npixels;

        // Bi-linear interpolation the rates in both energy layers
        double rate_emin = wgt_ll * m_background(m_inx_bgd, inx_ll + offset_emin) +
                           wgt_lr * m_background(m_inx_bgd, inx_lr + offset_emin) +
                           wgt_rl * m_background(m_inx_bgd, inx_rl + offset_emin) +
                           wgt_rr * m_background(m_inx_bgd, inx_rr + offset_emin);
        double rate_emax = wgt_ll * m_background(m_inx_bgd, inx_ll + offset_emax) +
                           wgt_lr * m_background(m_inx_bgd, inx_lr + offset_emax) +
                           wgt_rl * m_background(m_inx_bgd, inx_rl + offset_emax) +
                           wgt_rr * m_background(m_inx_bgd, inx_rr + offset_emax);

        // If both rates are positive then perform a logarithmic
        // interpolation in energy
        if (rate_emin > 0.0 && rate_emax > 0.0) {
            rate = std::exp(wgt_emin * std::log(rate_emin) +
                            wgt_emax * std::log(rate_emax));
        }
        else if (rate_emin > 0.0) {
            rate = rate_emin;
        }
        else if (rate_emax > 0.0) {
            rate = rate_emax;
        }

        // Make sure that background rate is not negative
        if (rate < 0.0) {
            rate = 0.0;
        }

    } // endif: background model was valid

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
 * Clears background.
 ***************************************************************************/
void GCTABackground3D::clear(void)
{
    // Free class members
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
 * @return Deep copy of background.
 *
 * Returns a pointer to a deep copy of the background.
 ***************************************************************************/
GCTABackground3D* GCTABackground3D::clone(void) const
{
    return new GCTABackground3D(*this);
}


/***********************************************************************//**
 * @brief Read background from FITS table
 *
 * @param[in] table FITS table.
 *
 * @exception GException::invalid_value
 *            Response table is not three-dimensional.
 *
 * Reads the background form the FITS @p table. The following column names
 * are mandatory:
 *
 *     DETX_LO  - DETX lower bin boundaries
 *     DETX_HI  - DETX upper bin boundaries
 *     DETY_LO  - DETY lower bin boundaries
 *     DETY_HI  - DETY upper bin boundaries
 *     ENERG_LO - Energy lower bin boundaries
 *     ENERG_HI - Energy upper bin boundaries
 *     BGD      - Background template
 *
 * The data are stored in the m_background member. The DETX and DETY axes
 * will be set to radians, the energy axis will be set to log10.
 ***************************************************************************/
void GCTABackground3D::read(const GFitsTable& table)
{
    // Clear response table
    m_background.clear();

    // Read background table
    m_background.read(table);

    // Get mandatory indices (throw exception if not found)
    m_inx_detx   = m_background.axis("DETX");
    m_inx_dety   = m_background.axis("DETY");
    m_inx_energy = m_background.axis("ENERG");
    m_inx_bgd    = m_background.table("BGD");

    // Throw an exception if the table is not three-dimensional
    if (m_background.axes() != 3) {
        std::string msg = "Expected three-dimensional background "
                          "response table but found "+
                          gammalib::str(m_background.axes())+
                          " dimensions. Please specify a three-dimensional "
                          "background.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Set DETX and DETY axis to radians
    m_background.axis_radians(m_inx_detx);
    m_background.axis_radians(m_inx_dety);

    // Set energy axis to logarithmic scale
    m_background.axis_log10(m_inx_energy);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Write background into FITS table
 *
 * @param[in] table FITS binary table.
 *
 * Writes background into a FITS binary @p table.
 *
 * @todo Add necessary keywords.
 ***************************************************************************/
void GCTABackground3D::write(GFitsBinTable& table) const
{
    // Write background table
    m_background.write(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load background from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the background from a FITS file.
 *
 * If no extension name is provided, the background will be loaded from the
 * "BACKGROUND" extension.
 ***************************************************************************/
void GCTABackground3D::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get background table
    const GFitsTable& table = *fits.table(filename.extname("BACKGROUND"));

    // Read effective area from table
    read(table);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save background into FITS file
 *
 * @param[in] filename Background table FITS file name.
 * @param[in] clobber Overwrite existing file? (default: false)
 *
 * Saves background into a FITS file. If a file with the given @p filename
 * does not yet exist it will be created, otherwise the method opens the
 * existing file. The method will create a (or replace an existing)
 * background extension. The extension name can be specified as part
 * of the @p filename, or if no extension name is given, is assumed to be
 * "BACKGROUND".
 *
 * An existing file will only be modified if the @p clobber flag is set to
 * true.
 ***************************************************************************/
void GCTABackground3D::save(const GFilename& filename, const bool& clobber) const
{
    // Get extension name
    std::string extname = filename.extname("BACKGROUND");

    // Open or create FITS file (without extension name since the requested
    // extension may not yet exist in the file)
    GFits fits(filename.url(), true);

    // Remove extension if it exists already
    if (fits.contains(extname)) {
        fits.remove(extname);
    }

    // Create binary table
    GFitsBinTable table;

    // Write the background table
    write(table);

    // Set binary table extension name
    table.extname(extname);

    // Append table to FITS file
    fits.append(table);

    // Save to file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in] energy Event energy.
 * @param[in] time Event trigger time.
 * @param[in,out] ran Random number generator.
 * @return CTA instrument direction.
 *
 * Returns a Monte Carlo simulated instrument direction for a given @p energy
 * and event trigger @p time. The simulation is done using a rejection
 * method that makes use of the background rate access operator. This assures
 * that the simulation is consistent with the background rate interpolation
 * that is done by the access operator.
 ***************************************************************************/
GCTAInstDir GCTABackground3D::mc(const GEnergy& energy,
                                 const GTime&   time,
                                 GRan&          ran) const
{
    // Initialise Monte Carlo Cache
    if (m_mc_max.empty()) {
        init_mc_cache();
    }

    // Allocate instrument direction
    GCTAInstDir dir;

    // Continue only if are maps
    if (!m_mc_max.empty()) {

        // Get reference to node array for the energy axis
        const GNodeArray& energy_nodes = m_background.axis_nodes(m_inx_energy);

        // Get log10(energy) in TeV
        double logE = energy.log10TeV();

        // Set values for node arrays
        energy_nodes.set_value(logE);

        // Get indices for energy interpolation
        int inx_left  = energy_nodes.inx_left();
        int inx_right = energy_nodes.inx_right();

        // Get maximum background rates
        double max_rate_left  = m_mc_max[inx_left];
        double max_rate_right = m_mc_max[inx_right];

        // Compute maximum background rate
        double max_rate = 0.0;
        if (max_rate_left > 0.0 && max_rate_right > 0.0) {
            double wgt_left  = energy_nodes.wgt_left();
            double wgt_right = energy_nodes.wgt_right();
            max_rate         = std::exp(wgt_left  * std::log(max_rate_left) +
                                        wgt_right * std::log(max_rate_right));
        }
        else if (max_rate_left > 0.0) {
            max_rate = max_rate_left;
        }
        else if (max_rate_right > 0.0) {
            max_rate = max_rate_right;
        }

        // Get instrument direction
        while (true) {

            // Get randomized detx and dety (radians)
            double detx = (m_mc_detx_min +
                           ran.uniform() * (m_mc_detx_max - m_mc_detx_min)) *
                          gammalib::deg2rad;
            double dety = (m_mc_dety_min +
                           ran.uniform() * (m_mc_dety_max - m_mc_dety_min)) *
                          gammalib::deg2rad;

            // Get background rate for these coordinates
            double value = this->operator()(logE, detx, dety);

            // Get uniform random number
            double uniform = ran.uniform() * max_rate;

            // Exit loop if we're not larger than the background rate value
            if (uniform <= value) {
                dir.detx(detx);
                dir.dety(dety);
                break;
            }

        } // endwhile: loop until instrument direction was accepted

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

        // Append header
        result.append("=== GCTABackground3D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);

        // If background model is valid then print information
        if (is_valid()) {

            // Compute DETX boundaries in deg
            double detx_min = m_background.axis_lo(m_inx_detx,0);
            double detx_max = m_background.axis_hi(m_inx_detx,
                              m_background.axis_bins(m_inx_detx)-1);

            // Compute DETY boundaries in deg
            double dety_min = m_background.axis_lo(m_inx_dety,0);
            double dety_max = m_background.axis_hi(m_inx_dety,
                              m_background.axis_bins(m_inx_dety)-1);

            // Compute energy boundaries in TeV
            double emin = m_background.axis_lo(m_inx_energy,0);
            double emax = m_background.axis_hi(m_inx_energy,
                          m_background.axis_bins(m_inx_energy)-1);

            // Append information
            result.append("\n"+gammalib::parformat("Number of DETX bins") +
                          gammalib::str(m_background.axis_bins(m_inx_detx)));
            result.append("\n"+gammalib::parformat("Number of DETY bins") +
                          gammalib::str(m_background.axis_bins(m_inx_dety)));
            result.append("\n"+gammalib::parformat("Number of energy bins") +
                          gammalib::str(m_background.axis_bins(m_inx_energy)));
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

        } // endif: there were 3 axis

        // ... otherwise show empty array
        else {
            result.append("\n"+gammalib::parformat("Number of DETX bins") +
                          gammalib::str(0));
            result.append("\n"+gammalib::parformat("Number of DETY bins") +
                          gammalib::str(0));
            result.append("\n"+gammalib::parformat("Number of energy bins") +
                          gammalib::str(0));
        }

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
    m_inx_detx    = 0;
    m_inx_dety    = 1;
    m_inx_energy  = 2;
    m_inx_bgd     = 0;

    // Initialise MC cache
    m_mc_max.clear();
    m_mc_spectrum.clear();
    m_mc_detx_min = 0.0;
    m_mc_detx_max = 0.0;
    m_mc_detx_bin = 0.0;
    m_mc_dety_min = 0.0;
    m_mc_dety_max = 0.0;
    m_mc_dety_bin = 0.0;
    m_mc_logE_min = 0.0;
    m_mc_logE_max = 0.0;
    m_mc_logE_bin = 0.0;

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
    m_inx_detx    = bgd.m_inx_detx;
    m_inx_dety    = bgd.m_inx_dety;
    m_inx_energy  = bgd.m_inx_energy;
    m_inx_bgd     = bgd.m_inx_bgd;

    // Copy MC cache
    m_mc_max      = bgd.m_mc_max;
    m_mc_spectrum = bgd.m_mc_spectrum;
    m_mc_detx_min = bgd.m_mc_detx_min;
    m_mc_detx_max = bgd.m_mc_detx_max;
    m_mc_detx_bin = bgd.m_mc_detx_bin;
    m_mc_dety_min = bgd.m_mc_dety_min;
    m_mc_dety_max = bgd.m_mc_dety_max;
    m_mc_dety_bin = bgd.m_mc_dety_bin;
    m_mc_logE_min = bgd.m_mc_logE_min;
    m_mc_logE_max = bgd.m_mc_logE_max;
    m_mc_logE_bin = bgd.m_mc_logE_bin;

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
 * @exception GException::invalid_value
 *            No valid background model defined.
 *
 * Initialises the cache for Monte Carlo sampling.
 * 
 * @todo Verify assumption made about the solid angles of the response table
 *       elements.
 * @todo Implement method that computes the "flux" in a background map
 *       pixel (inspired from GSkyMap::flux())
 ***************************************************************************/
void GCTABackground3D::init_mc_cache(void) const
{
    // Initialise cache
    m_mc_max.clear();
    m_mc_spectrum.clear();

    // Determine number of response cube pixels and maps
    int nx    = m_background.axis_bins(m_inx_detx);
    int ny    = m_background.axis_bins(m_inx_dety);
    int npix  = nx * ny;
    int nmaps = m_background.axis_bins(m_inx_energy);

    // Compute DETX boundaries and binsize in degress
    m_mc_detx_min = m_background.axis_lo(m_inx_detx,0);
    m_mc_detx_max = m_background.axis_hi(m_inx_detx,nx-1);
    m_mc_detx_bin = (m_mc_detx_max - m_mc_detx_min) / nx;

    // Compute DETY boundaries and binsize in degress
    m_mc_dety_min = m_background.axis_lo(m_inx_dety,0);
    m_mc_dety_max = m_background.axis_hi(m_inx_dety,ny-1);
    m_mc_dety_bin = (m_mc_dety_max - m_mc_dety_min) / ny;

    // Compute energy boundaries and binsize in log10(TeV)
    GEnergy emin(m_background.axis_lo(m_inx_energy,0),
                 m_background.axis_lo_unit(m_inx_energy));
    GEnergy emax(m_background.axis_hi(m_inx_energy,nmaps-1),
                 m_background.axis_hi_unit(m_inx_energy));
    m_mc_logE_min = emin.log10TeV();
    m_mc_logE_max = emax.log10TeV();
    m_mc_logE_bin = (m_mc_logE_max - m_mc_logE_min) / nmaps;

    // Determine solid angle of the pixel. We assume here simply that we
    // have square pixels of identical solid angle. It needs to be checked
    // whether this is a valid assumption.
    double solidangle = m_mc_detx_bin * gammalib::deg2rad *
                        m_mc_dety_bin * gammalib::deg2rad;

    // Compute flux normalisation factor. As we add-up 3 intensity we have
    // do divide by 3, and as we add up 8 wedges we have to divide the
    // solid angle by 8.
    double flux_norm  = gammalib::onethird * solidangle / 8.0;

    // Continue only if there are pixels and maps
    if (npix > 0 && nmaps > 0) {

        // Loop over all maps
        for (int i = 0; i < nmaps; ++i) {

            // Get logE
            double logE = m_mc_logE_min + (i+0.5) * m_mc_logE_bin;

            // Initialise cache with cumulative pixel fluxes and compute
            // total flux in response table for normalization. Negative
            // pixels are excluded from the cumulative map.
            double total_flux = 0.0; // units: events/s/MeV
            double max_rate   = 0.0; // units: events/s/sr/MeV
            double xmin       = m_mc_detx_min * gammalib::deg2rad;
            double xbin       = m_mc_detx_bin * gammalib::deg2rad;
            double ymin       = m_mc_dety_min * gammalib::deg2rad;
            double ybin       = m_mc_dety_bin * gammalib::deg2rad;
            double dx         = 0.5 * xbin;
            double dy         = 0.5 * ybin;
            double detx       = xmin + 0.5 * xbin;
            for (int ix = 0; ix < nx; ++ix, detx += xbin) {
                double dety = ymin + 0.5 * ybin;
                for (int iy = 0; iy < ny; ++iy, dety += ybin) {

                    // Compute flux by summing the flux in 8 pixel wedges
                    double i0    = (*this)(logE, detx,    dety);
                    double i1    = (*this)(logE, detx-dx, dety-dy);
                    double i2    = (*this)(logE, detx,    dety-dy);
                    double i3    = (*this)(logE, detx+dx, dety-dy);
                    double i4    = (*this)(logE, detx+dx, dety);
                    double i5    = (*this)(logE, detx+dx, dety+dy);
                    double i6    = (*this)(logE, detx,    dety+dy);
                    double i7    = (*this)(logE, detx-dx, dety+dy);
                    double i8    = (*this)(logE, detx-dx, dety);
                    double flux1 = flux_norm * (i1 + i2 + i0);
                    double flux2 = flux_norm * (i2 + i3 + i0);
                    double flux3 = flux_norm * (i3 + i4 + i0);
                    double flux4 = flux_norm * (i4 + i5 + i0);
                    double flux5 = flux_norm * (i5 + i6 + i0);
                    double flux6 = flux_norm * (i6 + i7 + i0);
                    double flux7 = flux_norm * (i7 + i8 + i0);
                    double flux8 = flux_norm * (i8 + i1 + i0);
                    double flux  = (flux1 + flux2  + flux3  + flux4 +
                                    flux5 + flux6  + flux7  + flux8);

                    // Sum flux
                    if (flux > 0.0) {
                        total_flux += flux;
                    }

                    // Get maximum rate
                    if (i0 > max_rate) {
                        max_rate = i0;
                    }
                }
            }

            // Append maximum rate
            m_mc_max.push_back(max_rate);

            // Set energy value (unit independent)
            GEnergy energy;
            energy.log10(logE, m_background.axis_lo_unit(m_inx_energy));

            // Only append node if flux is positive
            if (total_flux > 0.0) {
                m_mc_spectrum.append(energy, total_flux);
            }

            // Dump spectrum for debugging
            #if defined(G_DEBUG_MC_INIT)
            std::cout << "Energy=" << energy;
            std::cout << " Rate_max=" << max_rate;
            std::cout << " events/s/sr/MeV";
            std::cout << " Rate_total=" << total_flux;
            std::cout << " events/s/MeV" << std::endl;
            #endif

        } // endfor: looped over all maps

    } // endif: there were cube pixels and maps

    // Return
    return;
}
