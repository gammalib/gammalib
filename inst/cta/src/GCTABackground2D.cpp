/***************************************************************************
 *              GCTABackground2D.cpp - CTA 2D background class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GCTABackground2D.cpp
 * @brief CTA 2D background class implementation
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
#include "GCTAInstDir.hpp"
#include "GCTABackground2D.hpp"
#include "GCTASupport.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MC                  "GCTABackground2D::mc(GEnergy&, GTime&, GRan&)"
#define G_SET_MEMBERS                       "GCTABackground2D::set_members()"
#define G_INIT_MC_CACHE                   "GCTABackground2D::init_mc_cache()"

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
GCTABackground2D::GCTABackground2D(void) : GCTABackground()
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
GCTABackground2D::GCTABackground2D(const GFilename& filename) :
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
GCTABackground2D::GCTABackground2D(const GCTABackground2D& bgd) :
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
GCTABackground2D::~GCTABackground2D(void)
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
GCTABackground2D& GCTABackground2D::operator=(const GCTABackground2D& bgd)
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
 * @brief Return background rate in units of events MeV\f$^{-1}\f$
 *        s\f$^{-1}\f$ sr\f$^{-1}\f$
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] detx Tangential X coordinate in nominal system (radians).
 * @param[in] dety Tangential Y coordinate in nominal system (radians).
 * @return Background rate (events MeV\f$^{-1}\f$ s\f$^{-1}\f$ sr\f$^{-1}\f$).
 *
 * Returns the background rate for a given energy and detector coordinates.
 * The rate is given per ontime. The operator computes the offset angle
 *
 * \f[\theta = \sqrt{ {\rm DETX}^2 + {\rm DETY}^2 }\f]
 *
 * and interpolates linearly in \f$\theta\f$ and logarithmically in energy.
 *
 * The operator assures that the background rate never becomes negative. For
 * invalid background models or detector coordinates outside the range
 * covered by the response cube, the operator returns a rate of zero.
 ***************************************************************************/
double GCTABackground2D::operator()(const double& logE,
                                    const double& detx, 
                                    const double& dety) const
{
    // Initialise background rate
    double rate = 0.0;

    // Continue only if background model is valid
    if (is_valid()) {

        // Compute theta in radians
        double theta = std::sqrt(detx * detx + dety * dety);

        // Continue only if theta is within validity range
        if ((theta >= m_theta_min) && (theta <= m_theta_max)) {

            // Retrieve references to energy node array
            const GNodeArray& energy_nodes = m_background.axis_nodes(m_inx_energy);

            // Set energy for node array interpolation
            energy_nodes.set_value(logE);

            // Bi-linearly interpolate the rates in both energy layers
            double rate_emin = this->rate(energy_nodes.inx_left(),  theta);
            double rate_emax = this->rate(energy_nodes.inx_right(), theta);

            // If both rates are positive then perform a logarithmic interpolation
            // in energy
            if (rate_emin > 0.0 && rate_emax > 0.0) {

                // Perform energy interpolation
                rate = std::exp(energy_nodes.wgt_left()  * std::log(rate_emin) +
                                energy_nodes.wgt_right() * std::log(rate_emax));

                // Make sure that background rate is not negative
                if (rate < 0.0) {
                    rate = 0.0;
                }

            } // endif: both rates were positive

        } // endif: theta was valid

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
void GCTABackground2D::clear(void)
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
GCTABackground2D* GCTABackground2D::clone(void) const
{
    return new GCTABackground2D(*this);
}


/***********************************************************************//**
 * @brief Read background from FITS table
 *
 * @param[in] table FITS table.
 *
 * Reads the background form the FITS @p table. The following column names
 * are mandatory:
 *
 *     THETA_LO - THETA lower bin boundaries
 *     THETA_HI - THETA upper bin boundaries
 *     ENERG_LO - Energy lower bin boundaries
 *     ENERG_HI - Energy upper bin boundaries
 *     BKG      - Background template
 *
 * The data are stored in the m_background member. The ``THETA`` axis will
 * be set to radians, the energy axis will be set to log10. The data are
 * expected to be given in unit of events per MeV, second and steradian,
 * where the time is the real ontime.
 *
 * This method ignores all column names that are not the mandatory column
 * names in the FITS @p table.
 ***************************************************************************/
void GCTABackground2D::read(const GFitsTable& table)
{
    // Clear response table
    m_background.clear();

    // Copy response table and skip all columns that are not mandatory
    GFitsTable *ptr = table.clone();  // Make sure that table is of same type
    for (int i = 0; i < table.ncols(); ++i) {
        std::string colname = table[i]->name();
        if ((colname != "THETA_LO") && (colname != "THETA_HI") &&
            (colname != "ENERG_LO") && (colname != "ENERG_HI") &&
            (colname != "BKG")) {
            ptr->remove(colname);
        }
    }

    // Read background table from reduced FITS table. The background table
    // is expected to provide the number of events per MeV, second and
    // steradian, where the time is the real ontime.
    m_background.read(*ptr);

    // Delete reduced FITS table
    delete ptr;

    // Set class members
    set_members();

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
void GCTABackground2D::write(GFitsBinTable& table) const
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
 * If no extension name is given the method scans the `HDUCLASS` keywords
 * of all extensions and loads the background from the first extension
 * for which `HDUCLAS4=BKG_2D`.
 *
 * Otherwise, the background will be loaded from the `BKG` extension.
 ***************************************************************************/
void GCTABackground2D::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get the default extension name. If no GADF compliant name was found
    // then set the default extension name to "BKG".
    std::string extname = gammalib::gadf_hduclas4(fits, "BKG_2D");
    if (extname.empty()) {
        extname = gammalib::extname_cta_background2d;
    }

    // Get background table
    const GFitsTable& table = *fits.table(filename.extname(extname));

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
 * @param[in] clobber Overwrite existing file?
 *
 * Saves background into a FITS file. If a file with the given @p filename
 * does not yet exist it will be created, otherwise the method opens the
 * existing file. The method will create a (or replace an existing)
 * background extension. The extension name can be specified as part
 * of the @p filename, or if no extension name is given, is assumed to be
 * `BKG`.
 *
 * An existing file will only be modified if the @p clobber flag is set to
 * true.
 ***************************************************************************/
void GCTABackground2D::save(const GFilename& filename, const bool& clobber) const
{
    // Get extension name
    std::string extname = filename.extname(gammalib::extname_cta_background2d);

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
GCTAInstDir GCTABackground2D::mc(const GEnergy& energy,
                                 const GTime&   time,
                                 GRan&          ran) const
{
    // If Monte Carlo has not been initialised then initialise the cache now
    if (m_mc_max.empty()) {
        init_mc_cache();
    }

    // Determine energy range of response table
    GEnergy emin = m_mc_spectrum.energy(0);
    GEnergy emax = m_mc_spectrum.energy(m_mc_spectrum.nodes()-1);
    if (energy < emin || energy > emax) {
        std::string msg = "Event energy "+energy.print()+" is outside the "
                          "energy range ["+emin.print()+", "+emax.print()+"] "
                          "covered by the background response table. Please "
                          "restrict the energy range of the simulation to the "
                          "validity range of the background response table.";
        throw GException::invalid_value(G_MC, msg);
    }

    // Allocate instrument direction
    GCTAInstDir dir;

    // Continue only if there is a background model
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

        // Get maximum background rate as the maximum of the left and the right
        // energy node. Add some margin.
        double max_rate = m_mc_max[inx_left];
        if (m_mc_max[inx_right] > max_rate) {
            max_rate = m_mc_max[inx_right];
        }
        max_rate *= 1.5;

        // Get theta width
        double theta_width = 2.0 * m_theta_max;

        // Get instrument direction using a rejection method. This assures
        // that the Monte Carlo sample follows the model distribution
        while (true) {

            // Get randomized detx and dety in radians
            double detx = ran.uniform() * theta_width - m_theta_max;
            double dety = ran.uniform() * theta_width - m_theta_max;

            // Get background rate for these coordinates. If the background
            // rate is zero then get other coordinates.
            double value = this->operator()(logE, detx, dety);
            if (value <= 0.0) {
                continue;
            }

            // Dump a warning if value is larger than the maximum rate. This
            // should never happen!
            if (value > max_rate) {
                std::string msg = "Background rate "+gammalib::str(value)+" "
                                  "for "+energy.print()+" at "
                                  "DETXY=["+
                                  gammalib::str(detx*gammalib::rad2deg)+","+
                                  gammalib::str(dety*gammalib::rad2deg)+"] "
                                  "is larger than the maximum expected rate "+
                                  gammalib::str(max_rate)+". Something is "
                                  "wrong with the code.";
                gammalib::warning(G_MC, msg);
            }

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
 * @brief Returns background rate integrated over energy interval in units
 *        of events s\f$^{-1}\f$ sr\f$^{-1}\f$
 *
 * @param[in] dir Instrument direction.
 * @param[in] emin Minimum energy of energy interval.
 * @param[in] emax Maximum energy of energy interval.
 * @return Integrated background count rate (events s\f$^{-1}\f$ sr\f$^{-1}\f$).
 *
 * Returns the background count rate for a given instrument direction that
 * is integrated over a specified energy interval. The rate is given per
 * ontime.
 *
 * If the energy interval is not positive, a zero background rate is
 * returned.
 ***************************************************************************/
double GCTABackground2D::rate_ebin(const GCTAInstDir& dir,
                                   const GEnergy&     emin,
                                   const GEnergy&     emax) const
{
    // Initialise rate
    double rate = 0.0;

    // Continue only if background model is valid and the energy interval is
    // positive
    if (is_valid() && (emax > emin)) {

        // Compute theta in radians
        double theta = std::sqrt(dir.detx() * dir.detx() + dir.dety() * dir.dety());

        // Continue only if theta is within validity range
        if ((theta >= m_theta_min) && (theta <= m_theta_max)) {

            // Initialise first and second node
            double x1 = emin.MeV();
            double x2 = 0.0;
            double f1 = (*this)(emin.log10TeV(), dir.detx(), dir.dety());
            double f2 = 0.0;

            // Loop over all nodes
            for (int i = 0; i < m_energy.size(); ++i) {

                // If node energy is below x1 then skip node
                if (m_energy[i].MeV() <= x1) {
                    continue;
                }

                // If node energy is above emax then use emax as energy
                if (m_energy[i] > emax) {
                    x2 = emax.MeV();
                    f2 = (*this)(emax.log10TeV(), dir.detx(), dir.dety());
                }

                // ... otherwise use node energy
                else {
                    x2 = m_energy[i].MeV();
                    f2 = this->rate(i, theta);
                }

                // Compute integral
                if (f1 > 0.0 && f2 > 0.0) {
                    rate += gammalib::plaw_integral(x1, f1, x2, f2);
                }

                // Set second node as first node
                x1 = x2;
                f1 = f2;

                // If node energy is above emax then break now
                if (m_energy[i] > emax) {
                    break;
                }

            } // endfor: looped over all nodes

            // If last node energy is below emax then compute last part of
            // integral up to emax
            if (x1 < emax.MeV()) {
                x2    = emax.MeV();
                f2    = (*this)(emax.log10TeV(), dir.detx(), dir.dety());
                if (f1 > 0.0 && f2 > 0.0) {
                    rate += gammalib::plaw_integral(x1, f1, x2, f2);
                }
            }

        } // endif: theta was in valid range

    } // endif: background model was valid and energy interval was positive

    // Return background rate
    return rate;
}


/***********************************************************************//**
 * @brief Print background information
 *
 * @param[in] chatter Chattiness.
 * @return String containing background information.
 ***************************************************************************/
std::string GCTABackground2D::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTABackground2D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);

        // If background model is valid then print information
        if (is_valid()) {

            // Compute THETA boundaries in deg
            double theta_min = m_background.axis_lo(m_inx_theta, 0);
            double theta_max = m_background.axis_hi(m_inx_theta, m_num_theta-1);

            // Compute energy boundaries in TeV
            double emin = m_background.axis_lo(m_inx_energy, 0);
            double emax = m_background.axis_hi(m_inx_energy, m_num_energy-1);

            // Append information
            result.append("\n"+gammalib::parformat("Number of THETA bins") +
                          gammalib::str(m_num_theta));
            result.append("\n"+gammalib::parformat("Number of energy bins") +
                          gammalib::str(m_num_energy));
            result.append("\n"+gammalib::parformat("Energy range"));
            result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");
            result.append("\n"+gammalib::parformat("Offset angle range"));
            result.append(gammalib::str(theta_min)+" - "+gammalib::str(theta_max)+" deg");

        } // endif: there were 2 axis

        // ... otherwise show empty array
        else {
            result.append("\n"+gammalib::parformat("Number of THETA bins") +
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
void GCTABackground2D::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_background.clear();
    m_energy.clear();
    m_inx_theta  = 0;
    m_inx_energy = 1;
    m_inx_bgd    = 0;
    m_num_theta  = 0;
    m_num_energy = 0;
    m_num[0]     = 0;
    m_num[1]     = 0;
    m_theta_min  = 0.0;
    m_theta_max  = 0.0;
    m_logE_min   = 0.0;
    m_logE_max   = 0.0;

    // Initialise MC cache
    m_mc_max.clear();
    m_mc_spectrum.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bgd Background.
 ***************************************************************************/
void GCTABackground2D::copy_members(const GCTABackground2D& bgd)
{
    // Copy members
    m_filename   = bgd.m_filename;
    m_background = bgd.m_background;
    m_energy     = bgd.m_energy;
    m_inx_theta  = bgd.m_inx_theta;
    m_inx_energy = bgd.m_inx_energy;
    m_inx_bgd    = bgd.m_inx_bgd;
    m_num_theta  = bgd.m_num_theta;
    m_num_energy = bgd.m_num_energy;
    m_num[0]     = bgd.m_num[0];
    m_num[1]     = bgd.m_num[1];
    m_theta_min  = bgd.m_theta_min;
    m_theta_max  = bgd.m_theta_max;
    m_logE_min   = bgd.m_logE_min;
    m_logE_max   = bgd.m_logE_max;

    // Copy MC cache
    m_mc_max       = bgd.m_mc_max;
    m_mc_spectrum  = bgd.m_mc_spectrum;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTABackground2D::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set members from background table
 *
 * @exception GException::invalid_value
 *            Response table is not three-dimensional.
 *
 * Set class members based on the background table. The THETA axis
 * will be set to radians, the energy axis will be set to log10.
 ***************************************************************************/
void GCTABackground2D::set_members(void)
{
    // Throw an exception if the table is not two-dimensional
    if (m_background.axes() != 2) {
        std::string msg = "Expected two-dimensional background "
                          "response table but found "+
                          gammalib::str(m_background.axes())+
                          " dimensions. Please specify a two-dimensional "
                          "background.";
        throw GException::invalid_value(G_SET_MEMBERS, msg);
    }

    // Get mandatory indices (throw exception if not found)
    m_inx_theta  = m_background.axis("THETA");
    m_inx_energy = m_background.axis("ENERG");
    m_inx_bgd    = m_background.table("BKG");

    // Get axes dimensions
    m_num_theta  = m_background.axis_bins(m_inx_theta);
    m_num_energy = m_background.axis_bins(m_inx_energy);

    // Set dimension array
    m_num[m_inx_theta]  = m_num_theta;
    m_num[m_inx_energy] = m_num_energy;

    // Set THETA axis to radians
    m_background.axis_radians(m_inx_theta);

    // Set energy axis to logarithmic scale
    m_background.axis_log10(m_inx_energy);

    // Compute THETA boundaries in radians
    m_theta_min = m_background.axis_lo(m_inx_theta, 0) *
                  gammalib::deg2rad;
    m_theta_max = m_background.axis_hi(m_inx_theta, m_num_theta-1) *
                  gammalib::deg2rad;

    // Compute energy boundaries in log10(TeV)
    GEnergy emin(m_background.axis_lo(m_inx_energy, 0),
                 m_background.axis_lo_unit(m_inx_energy));
    GEnergy emax(m_background.axis_hi(m_inx_energy, m_num_energy-1),
                 m_background.axis_hi_unit(m_inx_energy));
    m_logE_min = emin.log10TeV();
    m_logE_max = emax.log10TeV();

    // Setup energy vector
    m_energy.clear();
    for (int i = 0; i < m_num_energy; ++i) {
        GEnergy elo(m_background.axis_lo(m_inx_energy, i),
                    m_background.axis_lo_unit(m_inx_energy));
        GEnergy ehi(m_background.axis_hi(m_inx_energy, i),
                    m_background.axis_hi_unit(m_inx_energy));
        GEnergy energy(std::sqrt(elo.TeV()*ehi.TeV()), "TeV");
        m_energy.append(energy);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return background rate bin index
 *
 * @param[in] itheta THETA index.
 * @param[in] iebin Energy index.
 * @return Background rate bin index.
 *
 * Returns the background rate bin index independent of the ordering.
 ***************************************************************************/
int GCTABackground2D::index(const int& itheta,
                            const int& iebin) const
{
    // Set index arrays
    int inx[2];
    inx[m_inx_theta]  = itheta;
    inx[m_inx_energy] = iebin;

    // Compute index
    int index = inx[0] + inx[1] * m_num[0];

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Initialise Monte Carlo cache
 *
 * @exception GException::invalid_value
 *            Background response table was empty.
 *
 * Initialises the cache for Monte Carlo sampling.
 ***************************************************************************/
void GCTABackground2D::init_mc_cache(void) const
{
    // Set number of subbins per energy bin of response cube
    const int nsubbins = 10;

    // Initialise cache
    m_mc_max.clear();
    m_mc_spectrum.clear();

    // Initialise maximum rate array
    init_mc_max_rate();

    // Compute DETX and DETY binsize in radians
    double detx_bin = 2.0 * m_theta_max / m_num_theta;
    double dety_bin = 2.0 * m_theta_max / m_num_theta;

    // Set energy bins for Monte Carlo cache spectrum. Make sure that the
    // actual energy nodes are part of the energy bins to avoid that sharp
    // structures in the background rate variations are missed.
    GEnergies energies;
    for (int i = 0; i < m_num_energy; ++i) {

        // If this is the first node then add energies below the first node
        // down to the lower energy limit of the background response
        if (i == 0) {
            GEnergy   emin(m_background.axis_lo(m_inx_energy, i),
                           m_background.axis_lo_unit(m_inx_energy));
            GEnergies enodes(nsubbins/2, emin, m_energy[i]);
            enodes.remove(enodes.size()-1); // avoids duplication of energy
            energies.extend(enodes);
        }

        // If this is not the last node then append nodes between current
        // node and next node
        if (i < m_num_energy-1) {
            GEnergies enodes(nsubbins, m_energy[i], m_energy[i+1]);
            enodes.remove(enodes.size()-1); // avoids duplication of energy
            energies.extend(enodes);
        }

        // ... otherwise this is the last node. Add energies above the last
        // node up to the upper energy limit of the background response
        else {
            GEnergy   emax(m_background.axis_hi(m_inx_energy, i),
                           m_background.axis_hi_unit(m_inx_energy));
            GEnergies enodes(nsubbins/2, m_energy[i], emax);
            energies.extend(enodes);
        }

    } // endfor: looped over energy bins of background response

    // Loop over all energy bins
    for (int i = 0; i < energies.size(); ++i) {

        // Get logE of energy bin
        double logE = energies[i].log10TeV();

        // Initialise cache with cumulative pixel fluxes and compute total flux
        // in response table for normalization. Negative pixels are excluded
        // from the cumulative map.
        double total_flux = 0.0; // units: events/s/MeV
        double dx         = 0.5 * detx_bin;
        double dy         = 0.5 * dety_bin;
        double detx       = -m_theta_max + dx;
        for (int ix = 0; ix < 2*m_num_theta; ++ix, detx += detx_bin) {
            double dety = -m_theta_max + dy;
            for (int iy = 0; iy < 2*m_num_theta; ++iy, dety += dety_bin) {

                // Compute intensities
                double i0 = (*this)(logE, detx,    dety);
                double i1 = (*this)(logE, detx-dx, dety-dy);
                double i2 = (*this)(logE, detx,    dety-dy);
                double i3 = (*this)(logE, detx+dx, dety-dy);
                double i4 = (*this)(logE, detx+dx, dety);
                double i5 = (*this)(logE, detx+dx, dety+dy);
                double i6 = (*this)(logE, detx,    dety+dy);
                double i7 = (*this)(logE, detx-dx, dety+dy);
                double i8 = (*this)(logE, detx-dx, dety);

                // Compute solid angle of the 8 pixel wedges
                double s1 = solid_angle(detx, dety, detx-dx, dety-dy, detx,    dety-dy);
                double s2 = solid_angle(detx, dety, detx,    dety-dy, detx+dx, dety-dy);
                double s3 = solid_angle(detx, dety, detx+dx, dety-dy, detx+dx, dety);
                double s4 = solid_angle(detx, dety, detx+dx, dety,    detx+dx, dety+dy);
                double s5 = solid_angle(detx, dety, detx+dx, dety+dy, detx,    dety+dy);
                double s6 = solid_angle(detx, dety, detx,    dety+dy, detx-dx, dety+dy);
                double s7 = solid_angle(detx, dety, detx-dx, dety+dy, detx-dx, dety);
                double s8 = solid_angle(detx, dety, detx-dx, dety,    detx-dx, dety-dy);

                // Compute flux by summing the flux in 8 pixel wedges
                double flux1 = s1 * (i1 + i2 + i0);
                double flux2 = s2 * (i2 + i3 + i0);
                double flux3 = s3 * (i3 + i4 + i0);
                double flux4 = s4 * (i4 + i5 + i0);
                double flux5 = s5 * (i5 + i6 + i0);
                double flux6 = s6 * (i6 + i7 + i0);
                double flux7 = s7 * (i7 + i8 + i0);
                double flux8 = s8 * (i8 + i1 + i0);
                double flux  = (flux1 + flux2  + flux3  + flux4 +
                                flux5 + flux6  + flux7  + flux8) / 3.0;

                // Sum flux
                if (flux > 0.0) {
                    total_flux += flux;
                }

            } // endfor: loop over DETY

        } // endfor: loop over DETX

        // Append node if flux is positive
        if (total_flux > 0.0) {
            m_mc_spectrum.append(energies[i], total_flux);
        }

        // Dump spectrum for debugging
        #if defined(G_DEBUG_MC_INIT)
        std::cout << "Energy=" << energies[i];
        std::cout << " Rate_total=" << total_flux;
        std::cout << " events/(s MeV)" << std::endl;
        #endif

    } // endfor: looped over energy bins

    // Make sure that spectrum is not empty (if this is the case the entire
    // background response table was empty)
    if (m_mc_spectrum.nodes() == 0) {
        std::string msg = "Background response table is empty. Please provide "
                          "a valid three-dimensional background template.";
        throw GException::invalid_value(G_INIT_MC_CACHE, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise array of maximum background rate
 *
 * Initialises the array m_mc_max that holds the maximum background rate for
 * each energy of the background model.
 ***************************************************************************/
void GCTABackground2D::init_mc_max_rate(void) const
{
    // Initialise cache
    m_mc_max.clear();

    // Loop over all maps to determine the maximum map intensity in units
    // of events/s/sr/MeV
    for (int iebin = 0; iebin < m_num_energy; ++iebin) {

        // Initialise maximum rate
        double max_rate = 0.0;

        // Loop over THETA angles
        for (int itheta = 0; itheta < m_num_theta; ++itheta) {

            // Get bin index
            int inx = index(itheta, iebin);

            // Get background rate
            double rate = m_background(m_inx_bgd, inx);

            // If background rate is larger than maximum then store the
            // background rate as new maximum
            if (rate > max_rate) {
                max_rate = rate;
            }

        } // endfor: looped over THETA

        // Append maximum rate
        m_mc_max.push_back(max_rate);

    } // endif: looped over all maps

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute solid angle of pixel wedge
 *
 * @param[in] detx1 DETX of first edge.
 * @param[in] dety1 DETY of first edge.
 * @param[in] detx2 DETX of second edge.
 * @param[in] dety2 DETY of second edge.
 * @param[in] detx3 DETX of third edge.
 * @param[in] dety3 DETY of third edge.
 * @return Solid angle (steradians).
 *
 * Estimate the solid angle subtended by 3 coordinates using Huilier's
 * theorem.
 *
 * Below, the definiton of the pixel cornes and sides are shown as used
 * within the code.
 *
 *             a12
 *         1---------2
 *         |        / 
 *         |       /  
 *         |      /   
 *         |     /    
 *      a13|    /a23
 *         |   / 
 *         |  / 
 *         | /
 *         |/
 *         3
 *
 ***************************************************************************/
double GCTABackground2D::solid_angle(const double& detx1, const double& dety1,
                                     const double& detx2, const double& dety2,
                                     const double& detx3, const double& dety3) const
{
    // Set sky directions
    double  theta1 = std::sqrt(detx1 * detx1 + dety1 * dety1);
    double  phi1   = std::atan2(dety1, detx1);
    double  theta2 = std::sqrt(detx2 * detx2 + dety2 * dety2);
    double  phi2   = std::atan2(dety2, detx2);
    double  theta3 = std::sqrt(detx3 * detx3 + dety3 * dety3);
    double  phi3   = std::atan2(dety3, detx3);
    GSkyDir dir1;
    GSkyDir dir2;
    GSkyDir dir3;
    dir1.radec(phi1, gammalib::pihalf - theta1);
    dir2.radec(phi2, gammalib::pihalf - theta2);
    dir3.radec(phi3, gammalib::pihalf - theta3);

    // Compute angular distances between pixel corners
    double a12 = dir1.dist(dir2);
    double a13 = dir1.dist(dir3);
    double a23 = dir2.dist(dir3);

    // Compute solid angle
    double s          = 0.5 * (a12 + a23 + a13);
    double solidangle = 4.0 * std::atan(std::sqrt(std::tan(0.5*s) *
                                                  std::tan(0.5*(s-a12)) *
                                                  std::tan(0.5*(s-a23)) *
                                                  std::tan(0.5*(s-a13))));

    // Return solid angle
    return solidangle;
}


/***********************************************************************//**
 * @brief Return background rate for a given energy bin and offset angle
 *        value (events/s/MeV/sr)
 *
 * @param[in] iebin Energy index.
 * @param[in] theta Offset angle in nominal system (radians).
 * @return Background rate (events/s/MeV/sr)
 ***************************************************************************/
double GCTABackground2D::rate(const int& iebin, const double& theta) const
{
    // Retrieve reference to node array
    const GNodeArray& thetas = m_background.axis_nodes(m_inx_theta);

    // Set offset angle value for node array interpolation
    thetas.set_value(theta);

    // Get bin indices
    int inx_left  = index(thetas.inx_left(),  iebin);
    int inx_right = index(thetas.inx_right(), iebin);

    // Bi-linearly interpolate the rates
    double rate = thetas.wgt_left()  * m_background(m_inx_bgd, inx_left) +
                  thetas.wgt_right() * m_background(m_inx_bgd, inx_right);

    // Make sure that background rate is not negative
    if (rate < 0.0) {
        rate = 0.0;
    }

    // Return
    return rate;
}
