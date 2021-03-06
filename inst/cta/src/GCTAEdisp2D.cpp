/***************************************************************************
 *            GCTAEdisp2D.cpp - CTA 2D energy dispersion class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2021 by Florent Forest                              *
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
 * @file GCTAEdisp2D.cpp
 * @brief CTA 2D energy dispersion class implementation
 * @author Florent Forest
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <vector>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GNdarray.hpp"
#include "GFft.hpp"
#include "GException.hpp"
#include "GIntegral.hpp"
#include "GFilename.hpp"
#include "GRan.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GCTAEdisp2D.hpp"
#include "GCTASupport.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                               "GCTAEdisp2D::read(GFitsTable&)"
#define G_MC   "GCTAEdisp2D::mc(GRan&, GEnergy&, double&, double&, double&, "\
                                                                   "double&)"
#define G_FETCH                                        "GCTAEdisp2D::fetch()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_SMOOTH_EDISP_KLUDGE    //!< Get rid of noise in energy dispersion
//#define G_SMOOTH_EDISP_SAVE_TEST //!< Save test matrix

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_HESS_RENORMALIZATION

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs empty energy dispersion.
 ***************************************************************************/
GCTAEdisp2D::GCTAEdisp2D(void) : GCTAEdisp()
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
 * Constructs energy dispersion from a FITS file.
 ***************************************************************************/
GCTAEdisp2D::GCTAEdisp2D(const GFilename& filename) : GCTAEdisp()
{

    // Initialise class members
    init_members();

    // Load energy dispersion from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] edisp Energy dispersion.
 *
 * Constructs energy dispersion by copying from another energy dispersion.
 ***************************************************************************/
GCTAEdisp2D::GCTAEdisp2D(const GCTAEdisp2D& edisp) : GCTAEdisp(edisp)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(edisp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destructs energy dispersion.
 ***************************************************************************/
GCTAEdisp2D::~GCTAEdisp2D(void)
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
 * @param[in] edisp Energy dispersion.
 * @return Energy dispersion.
 *
 * Assigns energy dispersion.
 ***************************************************************************/
GCTAEdisp2D& GCTAEdisp2D::operator=(const GCTAEdisp2D& edisp)
{
    // Execute only if object is not identical
    if (this != &edisp) {

        // Copy base class members
        this->GCTAEdisp::operator=(edisp);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(edisp);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return energy dispersion in units of MeV\f$^{-1}\f$
 *
 * @param[in] ereco Reconstructed photon energy.
 * @param[in] etrue True photon energy.
 * @param[in] theta Offset angle in camera system (radians).
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Energy dispersion (MeV\f$^{-1}\f$)
 *
 * Returns the energy dispersion
 *
 * \f[
 *    E_{\rm disp}(E_{\rm reco} | E_{\rm true}, \theta) =
 *    \frac{E_{\rm disp}(m | E_{\rm true}, \theta)}{E_{\rm true}}
 * \f]
 *
 * in units of MeV\f$^{-1}\f$ where
 * \f$E_{\rm reco}\f$ is the reconstructed energy,
 * \f$E_{\rm true}\f$ is the true energy, and
 * \f$\theta\f$ is the offset angle.
 ***************************************************************************/
double GCTAEdisp2D::operator()(const GEnergy& ereco,
                               const GEnergy& etrue,
                               const double&  theta,
                               const double&  phi,
                               const double&  zenith,
                               const double&  azimuth) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Initalize edisp
    double edisp = 0.0;

    // Get log10 of true photon energy
    double logEsrc = etrue.log10TeV();

    // Continue only if logEsrc and theta are in validity range
    if ((logEsrc >= m_logEsrc_min)  && (logEsrc <= m_logEsrc_max) &&
        (theta   >= m_theta_min)    && (theta   <= m_theta_max)) {

        // Compute migration (migra=Ereco/Etrue)
        double migra = ereco.MeV() / etrue.MeV();

        // Continue only if migration is in validity range
        if ((migra >= m_migra_min)  && (migra <= m_migra_max)) {

            // Setup argument vector
            double arg[3];
            arg[m_inx_etrue] = logEsrc;
            arg[m_inx_migra] = migra;
            arg[m_inx_theta] = theta;

            // Compute edisp
            edisp = m_edisp(m_inx_matrix, arg[0], arg[1], arg[2]);

            // Make sure that energy dispersion is non-negative
            if (edisp < 0.0) {
                edisp = 0.0;
            }

            // If energy dispersion is positive then divide the energy
            // dispersion value by true photon energy to get the energy
            // dispersion per MeV
            else {
                edisp /= etrue.MeV();
            }

        } // endif: migration was in valid range

    } // endif: true energy and offset angle was in valid range

    // Return
    return edisp;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear energy dispersion
 *
 * Clears energy dispersion.
 ***************************************************************************/
void GCTAEdisp2D::clear(void)
{
    // Free class members
    free_members();
    this->GCTAEdisp::free_members();

    // Initialise members
    this->GCTAEdisp::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone energy dispersion
 *
 * @return Deep copy of energy dispersion.
 *
 * Returns a pointer to a deep copy of the point spread function.
 ***************************************************************************/
GCTAEdisp2D* GCTAEdisp2D::clone(void) const
{
    return new GCTAEdisp2D(*this);
}


/***********************************************************************//**
 * @brief Assign response table
 *
 * @param[in] table Response table.
 *
 * Assigns the response table for an energy dispersion.
 ***************************************************************************/
void GCTAEdisp2D::table(const GCTAResponseTable& table)
{
    // Assign response table
    m_edisp = table;

    // Set table
    set_table();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy dispersion from FITS table
 *
 * @param[in] table FITS table.
 *
 * @exception GException::invalid_value
 *            Response table is not three-dimensional.
 *
 * Reads the energy dispersion form the FITS @p table. The following column
 * names are mandatory:
 *
 *     ENERG_LO - True energy lower bin boundaries (alternative name: ETRUE_LO)
 *     ENERG_HI - True energy upper bin boundaries (alternative name: ETRUE_HI)
 *     MIGRA_LO - Migration lower bin boundaries
 *     MIGRA_HI - Migration upper bin boundaries
 *     THETA_LO - Offset angle lower bin boundaries
 *     THETA_HI - Offset angle upper bin boundaries
 *     MATRIX   - Migration matrix
 *
 * The data are stored in the m_edisp member. The energy axis will be set to
 * \f$log_{10} E_{\rm true}\f$, the offset angle axis to radians.
 *
 * The method assures that the energy dispersion is properly normalised.
 ***************************************************************************/
void GCTAEdisp2D::read(const GFitsTable& table)
{
    // Clear response table
    m_edisp.clear();

    // Read energy dispersion table
    m_edisp.read(table);

    // Throw an exception if the table is not two-dimensional
    if (m_edisp.axes() != 3) {
        std::string msg = "Expected three-dimensional energy dispersion "
                          "response table but found "+
                          gammalib::str(m_edisp.axes())+
                          " dimensions. Please specify a three-dimensional "
                          "energy dispersion.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Set table
    set_table();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energy dispersion into FITS binary table
 *
 * @param[in] table FITS binary table.
 *
 * Writes the energy dispersion into the FITS binary @p table.
 *
 * @todo Add keywords.
 ***************************************************************************/
void GCTAEdisp2D::write(GFitsBinTable& table) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Write background table
    m_edisp.write(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy dispersion from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the energy dispersion from a FITS file.
 *
 * The method does not actually load the FITS file, which will only be loaded
 * on request. Only the filename is stored. See the fetch() method for the
 * actually loading of the energy dispersion.
 ***************************************************************************/
void GCTAEdisp2D::load(const GFilename& filename)
{
    // Clear object
    clear();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy dispersion table into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file?
 *
 * Saves energy dispersion into a FITS file. If a file with the given
 * @p filename does not yet exist it will be created, otherwise the method
 * opens the existing file. The method will create a (or replace an existing)
 * energy dispersion extension. The extension name can be specified as part
 * of the @p filename, or if no extension name is given, is assumed to be
 * `ENERGY DISPERSION`.
 *
 * An existing file will only be modified if the @p clobber flag is set to
 * true.
 ***************************************************************************/
void GCTAEdisp2D::save(const GFilename& filename, const bool& clobber) const
{
    // Get extension name
    std::string extname = filename.extname(gammalib::extname_cta_edisp2d);

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
 * @brief Simulate energy dispersion
 *
 * @param[in] ran Random number generator.
 * @param[in] etrue True photon energy.
 * @param[in] theta Offset angle in camera system (radians).
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Reconstructed energy.
 *
 * @exception GException::invalid_return_value
 *            No energy dispersion information found for parameters or
 *            energy dispersion matrix is empty.
 *
 * Draws reconstructed energy value given a true energy @p etrue and offset
 * angle @p theta. If no energy dispersion information is available the
 * method will return the true photon energy.
 ***************************************************************************/
GEnergy GCTAEdisp2D::mc(GRan&          ran,
                        const GEnergy& etrue,
                        const double&  theta,
                        const double&  phi,
                        const double&  zenith,
                        const double&  azimuth) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Initialise reconstructed event energy with true photon energy
    GEnergy ereco = etrue;

    // Get boundaries for reconstructed energy
    GEbounds ebounds = ereco_bounds(etrue, theta, phi, zenith, azimuth);
    double   emin    = ebounds.emin().TeV();
    double   emax    = ebounds.emax().TeV();

    // Throw an exception if minimum energy is equal or larger than maximum
    // energy (this means that no energy dispersion information is available
    // for that energy)
    if (emin >= emax) {
        std::string msg = "No valid energy dispersion information available "
                          "for true photon energy "+etrue.print()+", "
                          "offset angle "+
                          gammalib::str(theta*gammalib::rad2deg)+" deg "
                          "and azimuth angle "+
                          gammalib::str(phi*gammalib::rad2deg)+" deg. Either "
                          "provide an energy dispersion matrix covering these "
                          "parameters or restrict the parameter space.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Get maximum energy dispersion value (including some margin)
    double max_edisp = 1.5 * get_max_edisp(etrue, theta);

    // Throw an exception if maximum energy dispersion is zero
    if (max_edisp <= 0.0) {
        std::string msg = "Energy dispersion matrix is empty. Please provide "
                          "a valid energy dispersion matrix.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Initialise rejection method
    double    ewidth               = emax - emin;
    double    f                    = 0.0;
    double    ftest                = 1.0;
    const int max_subsequent_zeros = 10;

    // Find energy by rejection method
    while (ftest > f) {

        // Draw random observed energy and evaluate energy dispersion matrix
        // until a non-zero energy dispersion value is found. Subsequent zero
        // energy dispersion values may indicate that there is no valid energy
        // dispersion information available. After a limited number of
        // subsequent zeros the loop is terminated with an exception
        int zeros = 0;
        do {
            ereco.TeV(emin + ewidth * ran.uniform());
            f = operator()(ereco, etrue, theta, phi, zenith, azimuth);
            if (f == 0.0) {
                zeros++;
                if (zeros > max_subsequent_zeros) {
                    std::string msg = "No valid energy dispersion information "
                                      "available for true photon energy "+
                                      etrue.print()+", offset angle "+
                                      gammalib::str(theta*gammalib::rad2deg)+
                                      " deg and azimuth angle "+
                                      gammalib::str(phi*gammalib::rad2deg)+
                                      " deg. Either provide an energy "
                                      "dispersion matrix covering these "
                                      "parameters or restrict the parameter "
                                      "space.";
                    throw GException::invalid_return_value(G_MC, msg);
                }
            }
        } while (f == 0);

        // Get uniform random value between zero and the maximum energy
        // dispersion value. The loop is quit if the random number is smaller
        // than the energy dispersion value
        ftest = ran.uniform() * max_edisp;

    } // endwhile:

    // Return energy
    return ereco;
}


/***********************************************************************//**
 * @brief Return observed energy interval that contains the energy dispersion.
 *
 * @param[in] etrue True photon energy.
 * @param[in] theta Offset angle in camera system (radians).
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Observed energy boundaries.
 *
 * Returns the band of observed energies outside of which the energy
 * dispersion becomes negligible for a given true energy @p etrue and
 * offset angle @p theta.
 ***************************************************************************/
GEbounds GCTAEdisp2D::ereco_bounds(const GEnergy& etrue,
                                   const double&  theta,
                                   const double&  phi,
                                   const double&  zenith,
                                   const double&  azimuth) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Compute only if parameters changed
    if (!m_ereco_bounds_computed || theta != m_last_theta_ereco) {

        // Set computation flag
        m_ereco_bounds_computed = true;
        m_last_theta_ereco      = theta;
        m_last_etrue.TeV(0.0); // force update

        // Compute reconstructed energy bounds
        compute_ereco_bounds(theta, phi, zenith, azimuth);

    }

    // Search index only if true photon energy has changed
    if (etrue != m_last_etrue) {

        // Store last true photon energy
        m_last_etrue = etrue;

        // Get true energy in TeV
        double etrue_TeV = etrue.TeV();

        // Find right index with bisection
        int low  = 0;
        int high = m_ereco_bounds.size() - 1;
        while ((high-low) > 1) {
            int  mid = (low+high) / 2;
            double e = m_edisp.axis_lo(m_inx_etrue, mid);
            if (etrue_TeV < e) {
                high = mid;
            }
            else {
                low = mid;
            }
        }

        // Index found
        m_index_ereco = low;

    }

    // Return energy boundaries
    return m_ereco_bounds[m_index_ereco];
}


/***********************************************************************//**
 * @brief Return true energy interval that contains the energy dispersion.
 *
 * @param[in] ereco Reconstructed event energy.
 * @param[in] theta Offset angle in camera system (radians).
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return True energy boundaries.
 *
 * Returns the band of true photon energies outside of which the energy
 * dispersion becomes negligible for a given observed energy @p ereco and
 * offset angle @p theta.
 ***************************************************************************/
GEbounds GCTAEdisp2D::etrue_bounds(const GEnergy& ereco,
                                   const double&  theta,
                                   const double&  phi,
                                   const double&  zenith,
                                   const double&  azimuth) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Compute true energy boundaries only if parameters changed
    if (!m_etrue_bounds_computed || theta != m_last_theta_etrue) {

        // Set computation flag
        m_etrue_bounds_computed = true;
        m_last_theta_etrue        = theta;
        m_last_ereco.TeV(0.0); // force update

        // Compute true energy boundaries
        compute_etrue_bounds(theta, phi, zenith, azimuth);

    }

    // Search index only if reconstructed event energy has changed
    if (ereco != m_last_ereco) {

        // Store reconstructed event energy
        m_last_ereco = ereco;

        // Get reconstructed event energy in TeV
        double ereco_TeV = ereco.TeV();

        // Find right index with bisection
        int low  = 0;
        int high = m_etrue_bounds.size() - 1;
        while ((high-low) > 1) {
            int  mid = (low+high) / 2;
            double e = m_edisp.axis_lo(m_inx_etrue, mid);
            if (ereco_TeV < e) {
                high = mid;
            }
            else {
                low = mid;
            }
        }

        // Index found
        m_index_etrue = low;

    }

    // Return energy boundaries
    return m_etrue_bounds[m_index_etrue];
}


/***********************************************************************//**
 * @brief Fetch energy dispersion
 *
 * @exception GException::file_error
 *            File not found.
 *            Unable to load energy dispersion.
 *
 * Fetches the energy dispersion by reading it from a FITS file.
 *
 * If the filename contains no extension name the method scans the `HDUCLASS`
 * keywords of all extensions and loads the energy dispersion from the first
 * extension for which `HDUCLAS4=EDISP_2D`.
 *
 * Otherwise, the background will be loaded from the `ENERGY DISPERSION`
 * extension.
 *
 * This method does nothing if the energy dispersion is already loaded, if
 * there is nothing to fetch, or if the m_filename member is empty.
 *
 * The method is thread save. The method checks whether the file from which
 * the energy dispersion should be loaded actually exists.
 ***************************************************************************/
void GCTAEdisp2D::fetch(void) const
{
    // Continue only if energy dispersion has not yet been fetched
    if (!m_fetched) {

        // Continue only if the file name is not empty
        if (!m_filename.is_empty()) {

            // Throw an exception if the file does not exist
            if (!m_filename.exists()) {
                std::string msg = "File \""+m_filename+"\" not found. Cannot "
                                  "fetch energy dispersion. Maybe the file has "
                                  "been deleted in the meantime.";
                GException::file_error(G_FETCH, msg);
            }

            // Signal that energy dispersion will be fetched (has to come
            // before reading since the ebounds_obs() and ebounds_src() methods
            // will otherwise call the fetch() method recursively.
            m_fetched = true;

            // Initialise exception flag
            bool has_exception = false;

            // Load energy dispersion. Catch any exception. Put the code into
            // a critical zone as it might be called from within a parallelized
            // thread.
            #pragma omp critical(GCTAEdisp2D_fetch)
            {
            try {


                // Open FITS file
                GFits fits(m_filename);

                // Get the default extension name. If no GADF compliant name
                // was found then set the default extension name to
                // "ENERGY DISPERSION".
                std::string extname = gammalib::gadf_hduclas4(fits, "EDISP_2D");
                if (extname.empty()) {
                    extname = gammalib::extname_cta_edisp2d;
                }

                // Get energy dispersion table
                const GFitsTable& table = *fits.table(m_filename.extname(extname));

                // Read energy dispersion from table
                const_cast<GCTAEdisp2D*>(this)->read(table);

                // Close FITS file
                fits.close();

            }
            catch (...) {
                has_exception = true;
            }
            }

            // Throw an exception if an exception has occured
            if (has_exception) {
                std::string msg = "Unable to load energy dispersion from "
                                  "file \""+m_filename+"\".";
                throw GException::file_error(G_FETCH, msg);
            }

        } // endif: filename was not empty

    } // endif: energy dispersion had not yet been fetched

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return energy dispersion probability for reconstructed energy
 *        interval
 *
 * @param[in] ereco_min Minimum of reconstructed energy interval.
 * @param[in] ereco_max Maximum of reconstructed energy interval.
 * @param[in] etrue True energy.
 * @param[in] theta Offset angle (radians).
 * @return Integrated energy dispersion probability.
 *
 * Computes
 *
 * \f[
 *    \int_{E_{\rm reco}^{\rm min}}^{E_{\rm reco}^{\rm max}}
 *    E_{\rm disp}(E_{\rm reco} | E_{\rm true}, \theta) \, dE_{\rm reco}
 * \f]
 *
 * where
 * \f$E_{\rm reco}\f$ is the reconstructed energy,
 * \f$E_{\rm true}\f$ is the true energy and
 * \f$\theta\f$ is the offset angle.
 *
 * The method takes into account that the energy dispersion data are stored
 * in a matrix and uses the trapezoidal rule for integration of the tabulated
 * data. This assures the fastest possible integration.
 ***************************************************************************/
double GCTAEdisp2D::prob_erecobin(const GEnergy& ereco_min,
                                  const GEnergy& ereco_max,
                                  const GEnergy& etrue,
                                  const double&  theta) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Initalize probability
    double prob = 0.0;

    // Get log10 of true energy
    double logEsrc = etrue.log10TeV();

    // Get migration limits
    double migra_min = ereco_min / etrue;
    double migra_max = ereco_max / etrue;

    // Continue only if logEsrc and theta are in validity range and migration
    // interval overlaps with response table
    if ((logEsrc   >= m_logEsrc_min) && (logEsrc   <= m_logEsrc_max) &&
        (theta     >= m_theta_min)   && (theta     <= m_theta_max)   &&
        (migra_max >  m_migra_min)   && (migra_min <  m_migra_max)) {

        // Constrain migration limits to response table
        if (migra_min < m_migra_min) {
            migra_min = m_migra_min;
        }
        if (migra_max > m_migra_max) {
            migra_max = m_migra_max;
        }

        // Retrieve references to node arrays
        const GNodeArray& etrue_nodes = m_edisp.axis_nodes(m_inx_etrue);
        const GNodeArray& theta_nodes = m_edisp.axis_nodes(m_inx_theta);
        const GNodeArray& migra_nodes = m_edisp.axis_nodes(m_inx_migra);

        // Set logEsrc and theta for node array interpolation
        etrue_nodes.set_value(logEsrc);
        theta_nodes.set_value(theta);

        // Compute base indices
        int base_ll = table_index(etrue_nodes.inx_left(),  0, theta_nodes.inx_left());
        int base_lr = table_index(etrue_nodes.inx_left(),  0, theta_nodes.inx_right());
        int base_rl = table_index(etrue_nodes.inx_right(), 0, theta_nodes.inx_left());
        int base_rr = table_index(etrue_nodes.inx_right(), 0, theta_nodes.inx_right());

        // Get migration stride
        int stride = table_stride(m_inx_migra);

        // Initialise first and second node
        double x1 = migra_min;
        double f1 = table_value(base_ll, base_lr, base_rl, base_rr,
                                etrue_nodes.wgt_left(),
                                etrue_nodes.wgt_right(),
                                theta_nodes.wgt_left(),
                                theta_nodes.wgt_right(), x1);
        double x2 = 0.0;
        double f2 = 0.0;

        // Loop over all migration nodes
        for (int i = 0, offset = 0; i < migra_nodes.size(); ++i, offset += stride) {

            // If migration value is below migra_min then skip node
            if (migra_nodes[i] <= migra_min) {
                continue;
            }

            // If migration value is above maximum migration value then use
            // migra_max as migration value
            if (migra_nodes[i] > migra_max) {
                x2 = migra_max;
                f2 = table_value(base_ll, base_lr, base_rl, base_rr,
                                 etrue_nodes.wgt_left(),
                                 etrue_nodes.wgt_right(),
                                 theta_nodes.wgt_left(),
                                 theta_nodes.wgt_right(), x2);
            }

            // ... otherwise use migration value
            else {
                x2 = migra_nodes[i];
                f2 = table_value(base_ll, base_lr, base_rl, base_rr,
                                 etrue_nodes.wgt_left(),
                                 etrue_nodes.wgt_right(),
                                 theta_nodes.wgt_left(),
                                 theta_nodes.wgt_right(), offset);
            }

            // Compute integral
            prob += 0.5 * (f1 + f2) * (x2 - x1);

            // Set second node as first node
            x1 = x2;
            f1 = f2;

            // If node energy is above migra_max then break now
            if (migra_nodes[i] > migra_max) {
                break;
            }

        } // endfor: looped over all nodes

        // If last node energy is below migra_max then compute last part of
        // integral up to emax
        if (x1 < migra_max) {
            x2    = migra_max;
            f2    = table_value(base_ll, base_lr, base_rl, base_rr,
                                etrue_nodes.wgt_left(),
                                etrue_nodes.wgt_right(),
                                theta_nodes.wgt_left(),
                                theta_nodes.wgt_right(), x2);
            prob += 0.5 * (f1 + f2) * (x2 - x1);
        }

    } // endif: logEsrc, theta and migration range were valid

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Print energy dispersion information
 *
 * @param[in] chatter Chattiness.
 * @return String containing energy dispersion information.
 ***************************************************************************/
std::string GCTAEdisp2D::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAEdisp2D ===");
        result.append("\n"+gammalib::parformat("Filename")+m_filename);

        // Make sure that energy dispersion is online
        fetch();

        // Initialise information
        int    nebins     = 0;
        int    nmigrabins = 0;
        int    nthetabins = 0;
        double emin       = 0.0;
        double emax       = 0.0;
        double mmin       = 0.0;
        double mmax       = 0.0;
        double omin       = 0.0;
        double omax       = 0.0;

        // Extract information if there are axes in the response table
        if (m_edisp.axes() > 0) {
            nebins     = m_edisp.axis_bins(m_inx_etrue);
            nmigrabins = m_edisp.axis_bins(m_inx_migra);
            nthetabins = m_edisp.axis_bins(m_inx_theta);
            emin       = m_edisp.axis_lo(m_inx_etrue,0);
            emax       = m_edisp.axis_hi(m_inx_etrue,nebins-1);
            mmin       = m_edisp.axis_lo(m_inx_migra,0);
            mmax       = m_edisp.axis_hi(m_inx_migra,nmigrabins-1);
            omin       = m_edisp.axis_lo(m_inx_theta,0);
            omax       = m_edisp.axis_hi(m_inx_theta,nthetabins-1);
        }

        // Append information
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(nebins));
        result.append("\n"+gammalib::parformat("Number of migration bins") +
                      gammalib::str(nmigrabins));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(nthetabins));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");
        result.append("\n"+gammalib::parformat("Migration range"));
        result.append(gammalib::str(mmin)+" - "+gammalib::str(mmax));
        result.append("\n"+gammalib::parformat("Offset angle range"));
        result.append(gammalib::str(omin)+" - "+gammalib::str(omax)+" deg");

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
void GCTAEdisp2D::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_edisp.clear();
    m_fetched     = false;
    m_inx_etrue   = 0;
    m_inx_migra   = 1;
    m_inx_theta   = 2;
    m_inx_matrix  = 0;
    m_logEsrc_min = 0.0;
    m_logEsrc_max = 0.0;
    m_migra_min   = 0.0;
    m_migra_max   = 0.0;
    m_theta_min   = 0.0;
    m_theta_max   = 0.0;
    m_max_edisp.clear();

    // Initialise computation cache
    m_ereco_bounds_computed = false;
    m_etrue_bounds_computed = false;
    m_index_ereco           =     0;
    m_index_etrue           =     0;
    m_last_theta_ereco      =  -1.0;
    m_last_theta_etrue      =  -1.0;
    m_last_etrue.TeV(0.0);
    m_last_ereco.TeV(0.0);
    m_ereco_bounds.clear();
    m_etrue_bounds.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] edisp Energy dispersion.
 ***************************************************************************/
void GCTAEdisp2D::copy_members(const GCTAEdisp2D& edisp)
{
    // Copy members
    m_filename    = edisp.m_filename;
    m_edisp       = edisp.m_edisp;
    m_fetched     = edisp.m_fetched;
    m_inx_etrue   = edisp.m_inx_etrue;
    m_inx_migra   = edisp.m_inx_migra;
    m_inx_theta   = edisp.m_inx_theta;
    m_inx_matrix  = edisp.m_inx_matrix;
    m_logEsrc_min = edisp.m_logEsrc_min;
    m_logEsrc_max = edisp.m_logEsrc_max;
    m_migra_min   = edisp.m_migra_min;
    m_migra_max   = edisp.m_migra_max;
    m_theta_min   = edisp.m_theta_min;
    m_theta_max   = edisp.m_theta_max;
    m_max_edisp   = edisp.m_max_edisp;

    // Copy computation cache
    m_ereco_bounds_computed = edisp.m_ereco_bounds_computed;
    m_etrue_bounds_computed = edisp.m_etrue_bounds_computed;
    m_index_ereco           = edisp.m_index_ereco;
    m_index_etrue           = edisp.m_index_etrue;
    m_last_theta_ereco      = edisp.m_last_theta_ereco;
    m_last_theta_etrue      = edisp.m_last_theta_etrue;
    m_last_etrue            = edisp.m_last_etrue;
    m_last_ereco            = edisp.m_last_ereco;
    m_ereco_bounds          = edisp.m_ereco_bounds;
    m_etrue_bounds          = edisp.m_etrue_bounds;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEdisp2D::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get true energy
 *
 * @param[in] ietrue True energy index.
 * @return True energy.
 ***************************************************************************/
GEnergy GCTAEdisp2D::etrue(const int& ietrue) const
{
    // Compute true energy in TeV
    double etrue_TeV = std::sqrt(m_edisp.axis_lo(m_inx_etrue, ietrue) *
                                 m_edisp.axis_hi(m_inx_etrue, ietrue));

    // Set true energy
    GEnergy etrue;
    etrue.TeV(etrue_TeV);

    // Return true energy
    return etrue;
}


/***********************************************************************//**
 * @brief Get migration
 *
 * @param[in] imigra Migration index.
 * @return Migration.
 ***************************************************************************/
double GCTAEdisp2D::migra(const int& imigra) const
{
    // Compute migration
    double migra = 0.5 * (m_edisp.axis_hi(m_inx_migra, imigra) +
                          m_edisp.axis_lo(m_inx_migra, imigra));

    // Return migration
    return migra;
}


/***********************************************************************//**
 * @brief Get offset angle in radiaus
 *
 * @param[in] itheta Offset angle index.
 * @return Offset angle (radians).
 ***************************************************************************/
double GCTAEdisp2D::theta(const int& itheta) const
{
    // Compute offset angle in radians
    double theta = 0.5 * (m_edisp.axis_lo(m_inx_theta, itheta) +
                          m_edisp.axis_hi(m_inx_theta, itheta)) *
                         gammalib::deg2rad;

    // Return offset angle
    return theta;
}


/***********************************************************************//**
 * @brief Compute vector of reconstructed energy bounds
 *
 * @param[in] theta Offset angle (radians).
 * @param[in] phi Azimuth angle (radians).
 * @param[in] zenith Zenith angle (radians).
 * @param[in] azimuth Azimuth angle (radians).
 *
 * Computes for all true energies the energy boundaries of the reconstructed
 * energies covered by valid migration matrix elements. Only matrix elements
 * with values >= 1.0e-12 are considered as valid elements. In case that no
 * matrix elements are found of a given true energy, the interval of observed
 * energies will be set to [0,0] (i.e. an empty interval).
 ***************************************************************************/
void GCTAEdisp2D::compute_ereco_bounds(const double& theta,
                                       const double& phi,
                                       const double& zenith,
                                       const double& azimuth) const
{
    // Clear ebounds_obs vector
    m_ereco_bounds.clear();

    // Set epsilon
    const double eps = 1.0e-12;

    // Determine number of true energy and migration bins
    int netrue = m_edisp.axis_bins(m_inx_etrue);
    int nmigra = m_edisp.axis_bins(m_inx_migra);

    // Loop over true photon energy
    for (int ietrue = 0; ietrue < netrue; ++ietrue) {

        // Compute true photon energy
        GEnergy etrue = this->etrue(ietrue);

        // Initialise results
        double ereco_min = 0.0;
        double ereco_max = 0.0;
        bool   found_min = false;
        bool   found_max = false;

        // Find minimum boundary
        for (int imigra = 0; imigra < nmigra; ++imigra) {

            // Compute migra value
            double migra = this->migra(imigra);

            // Get matrix term
            double arg[3];
            arg[m_inx_etrue] = etrue.log10TeV();
            arg[m_inx_migra] = migra;
            arg[m_inx_theta] = theta;
            double edisp     = m_edisp(m_inx_matrix, arg[0], arg[1], arg[2]);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                found_min = true;
                ereco_min = migra * etrue.TeV();
                break;
            }

        } // endfor: find minimum boundary

        // Find maximum boundary
        for (int imigra = nmigra-1; imigra >= 0; imigra--) {

            // Compute migra value
            double migra = this->migra(imigra);

            // Get matrix term
            double arg[3];
            arg[m_inx_etrue] = etrue.log10TeV();
            arg[m_inx_migra] = migra;
            arg[m_inx_theta] = theta;
            double edisp     = m_edisp(m_inx_matrix, arg[0], arg[1], arg[2]);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                found_max = true;
                ereco_max = migra * etrue.TeV();
                break;
            }

        } // endfor: find maximum boundary

        // If we did not find boundaries then set the interval to
        // a zero interval for safety
        if (!found_min || !found_max) {
            ereco_min = 0.0;
            ereco_max = 0.0;
        }

        // Set energy boundaries
        GEnergy emin;
        GEnergy emax;
        emin.TeV(ereco_min);
        emax.TeV(ereco_max);

        // Add energy boundaries to vector
        m_ereco_bounds.push_back(GEbounds(emin, emax));

    } // endfor: looped over logEsrc

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute vector of true energy boundaries
 *
 * @param[in] theta Offset angle (radians).
 * @param[in] phi Azimuth angle (radians).
 * @param[in] zenith Zenith angle (radians).
 * @param[in] azimuth Azimuth angle (radians).
 *
 * Computes for all observed energies the energy boundaries of the true
 * energies covered by valid migration matrix elements. Only matrix elements
 * with values >= 1.0e-12 are considered as valid elements. In case that no
 * matrix elements are found of a given observed energy, the interval of true
 * energies will be set to [0,0] (i.e. an empty interval).
 ***************************************************************************/
void GCTAEdisp2D::compute_etrue_bounds(const double& theta,
                                       const double& phi,
                                       const double& zenith,
                                       const double& azimuth) const
{
    // Initialise energies
    GEnergy etrue;
    GEnergy ereco;

    // Clear ebounds_src vector
    m_etrue_bounds.clear();

    // Set epsilon
    const double eps = 1.0e-12;

    // Determine number of true energy bins
    int nereco  = m_edisp.axis_bins(m_inx_etrue);
    int nmigras = m_edisp.axis_bins(m_inx_migra);

    // Loop over reconstructed energy
    for (int iereco = 0; iereco < nereco; ++iereco) {

        // Set ereco (we use the true energy axis for the computation)
        GEnergy ereco = this->etrue(iereco);

        // Initialise results
        GEnergy etrue_min;
        GEnergy etrue_max;
        bool    found_min   = false;
        bool    found_max   = false;

        // Find minimum boundary
        for (int imigra = 0; imigra < nmigras; ++imigra) {

            // Compute migra value
            double migra = this->migra(imigra);

            // Fall through if migration is zero
            if (migra == 0.0) {
                continue;
            }

            // Compute true energy
            etrue = ereco / migra;

            // Get energy dispersion
            double edisp = operator()(ereco, etrue, theta);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                found_max = true;
                etrue_max = etrue;
                break;
            }

        } // endfor: find minimum boundary

        // Find maximum boundary
        for (int imigra = nmigras-1; imigra >= 0; imigra--) {

            // Compute migra value
            double migra = this->migra(imigra);

            // Fall through if migration is zero
            if (migra == 0.0) {
                continue;
            }

            // Compute true energy
            etrue = ereco / migra;

            // Get energy dispersion
            double edisp = operator()(ereco, etrue, theta);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                found_min = true;
                etrue_min = etrue;
                break;
            }

        } // endfor: find maximum boundary

        // If we did not find boundaries then set the interval to a zero
        // interval for safety
        if (!found_min || !found_max) {
            etrue_min.clear();
            etrue_max.clear();
        }

        // Add energy boundaries to vector
        m_etrue_bounds.push_back(GEbounds(etrue_min, etrue_max));

    } // endfor : loopend over observed energy

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set table
 *
 * After assigning or loading a response table, performs all necessary steps
 * to set the table.
 *
 * For CTA, the method applies a kludge that smoothes the energy dispersion
 * table to get rid of numerical fluctuations.
 *
 * Sets the data members
 *      m_inx_etrue
 *      m_inx_migra
 *      m_inx_theta
 *      m_inx_matrix
 *      m_logEsrc_min
 *      m_logEsrc_max
 *      m_migra_min
 *      m_migra_max
 *      m_theta_min
 *      m_theta_max
 *      m_max_edisp
 ***************************************************************************/
void GCTAEdisp2D::set_table(void)
{
    // Set table indices
    if (m_edisp.has_axis("ENERG")) {
        m_inx_etrue = m_edisp.axis("ENERG");
    }
    else {
        m_inx_etrue = m_edisp.axis("ETRUE"); // Old name, should not be used
    }
    m_inx_migra  = m_edisp.axis("MIGRA");
    m_inx_theta  = m_edisp.axis("THETA");
    m_inx_matrix = m_edisp.table("MATRIX");

    // Set true energy axis to logarithmic scale
    m_edisp.axis_log10(m_inx_etrue);

    // Set offset angle axis to radians
    m_edisp.axis_radians(m_inx_theta);

    // Set table boundaries
    set_boundaries();

    // Smooth energy dispersion table (kludge only for CTA)
    #if defined(G_SMOOTH_EDISP_KLUDGE)
    if (gammalib::strip_whitespace(m_edisp.telescope()) == "CTA") {
        smooth_table();
    }
    #endif

    // Normalize energy dispersion table
    normalize_table();

    // Set maximum energy dispersion value array
    set_max_edisp();

    // Save test
    #if defined(G_SMOOTH_EDISP_SAVE_TEST)
    this->save("test_edisp.fits", true);
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy dispersion boundaries
 *
 * Sets the data members m_logEsrc_min, m_logEsrc_max, m_migra_min,
 * m_migra_max, m_theta_min and m_theta_max that define the validity range
 * of the energy dispersion.
 ***************************************************************************/
void GCTAEdisp2D::set_boundaries(void)
{
    // Get number of bins
    int netrue = m_edisp.axis_bins(m_inx_etrue);
    int nmigra = m_edisp.axis_bins(m_inx_migra);
    int ntheta = m_edisp.axis_bins(m_inx_theta);

    // Compute minimum and maximum logEsrc boundaries
    m_logEsrc_min = std::log10(m_edisp.axis_lo(m_inx_etrue, 0));
    m_logEsrc_max = std::log10(m_edisp.axis_hi(m_inx_etrue, netrue-1));

    // Compute minimum and maximum migration boundaries
    m_migra_min = m_edisp.axis_lo(m_inx_migra, 0);
    m_migra_max = m_edisp.axis_hi(m_inx_migra, nmigra-1);

    // Compute minimum and maximum theta boundaries
    m_theta_min = m_edisp.axis_lo(m_inx_theta, 0)        * gammalib::deg2rad;
    m_theta_max = m_edisp.axis_hi(m_inx_theta, ntheta-1) * gammalib::deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set array of maximum energy dispersion values
 *
 * Sets an internal array of maximum energy dispersion values for each given
 * true photon energy and offset angle in the response table. This array will
 * be used for Monte Carlo simulations.
 ***************************************************************************/
void GCTAEdisp2D::set_max_edisp(void)
{
    // Get number of bins
    int netrue = m_edisp.axis_bins(m_inx_etrue);
    int nmigra = m_edisp.axis_bins(m_inx_migra);
    int ntheta = m_edisp.axis_bins(m_inx_theta);

    // Initialise maximum energy dispersion
    m_max_edisp.assign(netrue*ntheta, 0.0);

    // Loop over offset angle
    for (int itheta = 0, index = 0; itheta < ntheta; ++itheta) {

        // Loop over true photon energy
        for (int ietrue = 0; ietrue < netrue; ++ietrue, ++index) {

            // Get true photon energy
            GEnergy etrue = this->etrue(ietrue);

            // Get offset and stride
            int offset = table_index(ietrue, 0, itheta);
            int stride = table_stride(m_inx_migra);

            // Find maximum
            double maximum = 0.0;
            for (int imigra = 0, i = offset; imigra < nmigra; ++imigra, i += stride) {
                double value = m_edisp(m_inx_matrix, i);
                if (value > maximum) {
                    maximum = value;
                }
            }

            // Set maximum
            m_max_edisp[index] = maximum / etrue.MeV();

        } // endfor: looped over true photon energy

    } // endfor: looped over offset angle

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get maximum energy dispersion value
 *
 * @param[in] etrue True photon energy.
 * @param[in] theta Offset angle (radians).
 * @return Maximum energy dispersion value.
 *
 * For a given true energy @p etrue and offset angle @p theta, returns the
 * maximum energy dispersion value that may occur. This method makes use
 * of the internal array pre-computed by set_max_edisp().
 *
 * If set_max_edisp() was not called, the method returns 0.
 ***************************************************************************/
double GCTAEdisp2D::get_max_edisp(const GEnergy& etrue,
                                  const double&  theta) const
{
    // Initialise maximum energy dispersion
    double max_edisp = 0.0;

    // Continue only if maximum energy dispersion array is set
    if (m_max_edisp.size() > 0) {

        // Get number of true energy bins
        int netrue = m_edisp.axis_bins(m_inx_etrue);

        // Get references to axis nodes
        const GNodeArray& etrue_nodes = m_edisp.axis_nodes(m_inx_etrue);
        const GNodeArray& theta_nodes = m_edisp.axis_nodes(m_inx_theta);

        // Set values
        etrue_nodes.set_value(etrue.log10TeV());
        theta_nodes.set_value(theta);

        // Get indices that bound the true energy and offset angle
        int inx_left  = theta_nodes.inx_left()  * netrue;
        int inx_right = theta_nodes.inx_right() * netrue;
        int index_1   = inx_left  + etrue_nodes.inx_left();
        int index_2   = inx_left  + etrue_nodes.inx_right();
        int index_3   = inx_right + etrue_nodes.inx_left();
        int index_4   = inx_right + etrue_nodes.inx_right();

        // Get maximum energy dispersion
        max_edisp          = m_max_edisp[index_1];
        double max_edisp_2 = m_max_edisp[index_2];
        double max_edisp_3 = m_max_edisp[index_3];
        double max_edisp_4 = m_max_edisp[index_4];
        if (max_edisp_2 > max_edisp) {
            max_edisp = max_edisp_2;
        }
        if (max_edisp_3 > max_edisp) {
            max_edisp = max_edisp_3;
        }
        if (max_edisp_4 > max_edisp) {
            max_edisp = max_edisp_4;
        }

    } // endif: Maximum energy dispersion array was set

    // Return maximum energy dispersion
    return max_edisp;
}


/***********************************************************************//**
 * @brief Normalize energy dispersion table
 *
 * Normalize the energy dispersion table using
 *
 * \f[
 *    \int_{E_{\rm reco}^{\rm min}}^{E_{\rm reco}^{\rm max}}
 *    E_{\rm disp}(E_{\rm reco} | E_{\rm true}, \theta) \, dE_{\rm reco} = 1
 * \f]
 *
 * where
 * \f$E_{\rm reco}_{\rm min}\f$ and \f$E_{\rm reco}_{\rm max}\f$ is the
 * minimum and maximum migration value for a given true energy as returned
 * by the ebounds_obs() method, and
 * \f$E_{\rm disp}(E_{\rm true}, E_{\rm reco}, \theta)\f$ is the energy
 * dispersion returned by the operator(), given in units of \f$MeV^{-1}\f$.
 *
 * The normalisation is performed for each \f$E_{\rm true}\f$ and
 * \f$\theta\f$ bin using a Romberg integration method. Since the
 * normalisation may affect the energy boundaries
 * \f$E_{\rm reco}^{\rm min}\f$ and \f$E_{\rm reco}^{\rm max}\f$, two
 * normalisation passes are performed to stabilize the result.
 ***************************************************************************/
void GCTAEdisp2D::normalize_table(void)
{
    // Get axes dimensions
    int etrue_size = m_edisp.axis_bins(m_inx_etrue);
    int migra_size = m_edisp.axis_bins(m_inx_migra);
    int theta_size = m_edisp.axis_bins(m_inx_theta);

    // Now normalize the migration matrix vectors. We do this by integrating
    // for each true energy and offset angle the measured energy. We perform
    // here two passes as the ebounds_obs() method depends on the content of
    // the matrix, and by having two passes we are sure that we get a
    // consistent energy interval.
    for (int pass = 0; pass < 2; ++pass) {

        // Initialise sums
        std::vector<double> sums;
        sums.reserve(theta_size*etrue_size);

        // Loop over offset angle.
        for (int i_theta = 0; i_theta < theta_size; ++i_theta) {

            // Get offset angle (in radians)
            double theta = this->theta(i_theta);

            // Loop over true photon energy
            for (int i_etrue = 0; i_etrue < etrue_size; ++i_etrue) {

                // Get true photon energy
                GEnergy etrue = this->etrue(i_etrue);

                // Initialise integration
                double sum = 0.0;

                // Get integration boundaries
                GEbounds ebounds = ereco_bounds(etrue, theta);

                // Loop over all energy intervals
                for (int i = 0; i < ebounds.size(); ++i) {

                    // Get energy boundaries in log10 of energy in MeV
                    double emin = ebounds.emin(i).log10MeV();
                    double emax = ebounds.emax(i).log10MeV();

                    // Setup integration function
                    edisp_ereco_kern integrand(this, etrue, theta);
                    GIntegral        integral(&integrand);

                    // Set integration precision
                    integral.eps(1.0e-6);

                    // Do Romberg integration
                    sum += integral.romberg(emin, emax);

                } // endfor: looped over energy intervals

                // Store sum
                sums.push_back(sum);

            } // endfor: looped over Etrue

        } // endfor: looped over theta

        // Normalize vectors of migration matrix for all etrue and theta.
        for (int i_theta = 0, inx = 0; i_theta < theta_size; ++i_theta) {
            for (int i_etrue = 0; i_etrue < etrue_size; ++i_etrue, ++inx) {
                double sum = sums[inx];
                if (sum > 0.0) {
                    int offset = table_index(i_etrue, 0, i_theta);
                    int stride = table_stride(m_inx_migra);
                    for (int k = 0, i = offset; k < migra_size; ++k, i += stride) {
                        m_edisp(m_inx_matrix,i) /= sum;
                    }
                }
            }
        }

    } // endfor: made 2 passes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return index of response table element
 *
 * @param[in] ietrue True energy index.
 * @param[in] imigra Migration index.
 * @param[in] itheta Offset index.
 * @return Index of response table element.
 ***************************************************************************/
int GCTAEdisp2D::table_index(const int& ietrue,
                             const int& imigra,
                             const int& itheta) const
{
    // Set index vector
    int inx[3];
    inx[m_inx_etrue] = ietrue;
    inx[m_inx_migra] = imigra;
    inx[m_inx_theta] = itheta;

    // Compute index
    int index = inx[0] + (inx[1] + inx[2] * m_edisp.axis_bins(1)) *
                         m_edisp.axis_bins(0);

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Return stride of response table axis
 *
 * @param[in] axis Response table axis.
 * @return Stride of response table axis.
 ***************************************************************************/
int GCTAEdisp2D::table_stride(const int& axis) const
{
    // Initialise stride
    int stride = 1;

    // Multiply in higher
    if (axis > 0) {
        stride *= m_edisp.axis_bins(0);
    }
    if (axis > 1) {
        stride *= m_edisp.axis_bins(1);
    }

    // Return stride
    return stride;
}


/***********************************************************************//**
 * @brief Return bi-linearly interpolate table value for given migration bin
 *
 * @param[in] base_ll Base index for left true energy and left offset angle.
 * @param[in] base_lr Base index for left true energy and right offset angle.
 * @param[in] base_rl Base index for right true energy and left offset angle.
 * @param[in] base_rr Base index for right true energy and right offset angle.
 * @param[in] wgt_el Weighting for left true energy.
 * @param[in] wgt_er Weighting for right true energy.
 * @param[in] wgt_tl Weighting for left offset angle.
 * @param[in] wgt_tr Weighting for right offset angle.
 * @param[in] offset Offset of migration bin with respect to base indices.
 * @return Bi-linearly interpolate table value.
 ***************************************************************************/
double GCTAEdisp2D::table_value(const int&    base_ll,
                                const int&    base_lr,
                                const int&    base_rl,
                                const int&    base_rr,
                                const double& wgt_el,
                                const double& wgt_er,
                                const double& wgt_tl,
                                const double& wgt_tr,
                                const int&    offset) const
{
    // Compute table element indices
    int inx_ll = base_ll + offset;
    int inx_lr = base_lr + offset;
    int inx_rl = base_rl + offset;
    int inx_rr = base_rr + offset;

    // Get
    double value = wgt_el * (m_edisp(m_inx_matrix, inx_ll) * wgt_tl  +
                             m_edisp(m_inx_matrix, inx_lr) * wgt_tr) +
                   wgt_er * (m_edisp(m_inx_matrix, inx_rl) * wgt_tl +
                             m_edisp(m_inx_matrix, inx_rr) * wgt_tr);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return bi-linearly interpolate table value for given migration value
 *
 * @param[in] base_ll Base index for left true energy and left offset angle.
 * @param[in] base_lr Base index for left true energy and right offset angle.
 * @param[in] base_rl Base index for right true energy and left offset angle.
 * @param[in] base_rr Base index for right true energy and right offset angle.
 * @param[in] wgt_el Weighting for left true energy.
 * @param[in] wgt_er Weighting for right true energy.
 * @param[in] wgt_tl Weighting for left offset angle.
 * @param[in] wgt_tr Weighting for right offset angle.
 * @param[in] migra Migration value.
 * @return Bi-linearly interpolate table value.
 ***************************************************************************/
double GCTAEdisp2D::table_value(const int&    base_ll,
                                const int&    base_lr,
                                const int&    base_rl,
                                const int&    base_rr,
                                const double& wgt_el,
                                const double& wgt_er,
                                const double& wgt_tl,
                                const double& wgt_tr,
                                const double& migra) const
{
    // Retrieve references to migration node array
    const GNodeArray& migra_nodes = m_edisp.axis_nodes(m_inx_migra);

    // Set migration value for node array interpolation
    migra_nodes.set_value(migra);

    // Get migration stride
    int stride = table_stride(m_inx_migra);

    // Get left value
    double value_left = table_value(base_ll, base_lr, base_rl, base_rr,
                                    wgt_el, wgt_er, wgt_tl, wgt_tr,
                                    migra_nodes.inx_left() * stride);

    // Get right value
    double value_right = table_value(base_ll, base_lr, base_rl, base_rr,
                                     wgt_el, wgt_er, wgt_tl, wgt_tr,
                                     migra_nodes.inx_right() * stride);

    // Interpolate result
    double value = migra_nodes.wgt_left()  * value_left +
                   migra_nodes.wgt_right() * value_right;

    // Make sure that interpolation result is not negative
    if (value < 0.0) {
        value = 0.0;
    }

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Smoothed energy dispersion table
 *
 * This method implements the kludge for CTA that reduces the statistical
 * noise in the energy dispersion matrix.
 ***************************************************************************/
void GCTAEdisp2D::smooth_table(void)
{
    // Get axes dimensions
    int netrue = m_edisp.axis_bins(m_inx_etrue);
    int nmigra = m_edisp.axis_bins(m_inx_migra);
    int ntheta = m_edisp.axis_bins(m_inx_theta);
    int npix   = netrue * nmigra;

    // Loop over all offset angles
    for (int itheta = 0; itheta < ntheta; ++itheta) {

        // Compute means, rms and total as function of etrue
        GNdarray mean(netrue);
        GNdarray rms(netrue);
        GNdarray total(netrue);
        get_moments(itheta, &mean, &rms, &total);

        // Smooth all three
        mean  = smooth_array(mean,  30, 30, 0.5);
        rms   = smooth_array(rms,   30, 30, 0.03);
        total = smooth_array(total, 30, 30, 0.0);

        // Make sure that total is not negative
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {
            if (total(ietrue) < 0.0) {
                total(ietrue) = 0.0;
            }
        }

        // Replace matrix by Gaussians
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {

            // Continue only if there is information
            if (total(ietrue) > 0.0) {

                // Get Gaussian
                GNdarray gauss = gaussian_array(mean(ietrue), rms(ietrue), total(ietrue));

                // Compute base index
                int ibase = itheta * npix + ietrue;
                if ((mean(ietrue) > 0.0) && (rms(ietrue))) {
                    for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
                        m_edisp(m_inx_matrix, i) = gauss(imigra);
                    }
                }
                else {
                    for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
                        m_edisp(m_inx_matrix, i) = 0.0;
                    }
                }

            } // endif: there was information

            // ... otherwise reset the matrix to zero
            else {
                int ibase = itheta * npix + ietrue;
                for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
                    m_edisp(m_inx_matrix, i) = 0.0;
                }
            }

        } // endfor: looped over true energies

    } // endfor: looped over offset angles

    // Return
    return;
}


/***********************************************************************//**
 * @brief Smoothed array
 *
 * @param[in] array Array.
 * @param[in] nstart Number of nodes used for regression at first array value.
 * @param[in] nstop Number of nodes used for regression at last array value.
 * @param[in] limit Limit for array element exclusion.
 * @return Smoothed array.
 *
 * Returns a smoothed array that is computed by locally adjusting a linear
 * regression law to the array elements.
 ***************************************************************************/
GNdarray GCTAEdisp2D::smooth_array(const GNdarray& array,
                                   const int&      nstart,
                                   const int&      nstop,
                                   const double&   limit) const
{
    // Initialise smoothed array
    GNdarray smoothed_array(array.size());

    // Compute node step size
    double nstep = double(nstop - nstart)/double(array.size());

    // Loop over all array elements and computed the smoothed value
    for (int i = 0; i < array.size(); ++i) {
        int nodes         = nstart + int(i*nstep);
        smoothed_array(i) = smoothed_array_value(i, array, nodes, limit);
    }

    // Return smoothed array
    return smoothed_array;
}


/***********************************************************************//**
 * @brief Get smoothed array value
 *
 * @param[in] inx Index of array element.
 * @param[in] array Array.
 * @param[in] nodes Number of nodes used for regression.
 * @param[in] limit Limit for array element exclusion.
 * @return Smoothed array value.
 *
 * Returns a smoothed array value that is derived from adjusting a linear
 * law using regression to @p nodes adjacent array elements. Array elements
 * with values below @p limit are excluded from the regression.
 ***************************************************************************/
double GCTAEdisp2D::smoothed_array_value(const int&      inx,
                                         const GNdarray& array,
                                         const int&      nodes,
                                         const double&   limit) const
{
    // Initialise variables
    double mean_x = 0.0;
    double mean_y = 0.0;
    int    ileft  = inx - 1;
    int    iright = inx + 1;
    int    nleft  = 0;
    int    nright = 0;

    // Allocate vector of elements used for regression
    std::vector<double> x_vals;
    std::vector<double> y_vals;

    // Add nodes on the left of the element of interest
    while (nleft < nodes) {
        if (ileft < 0) {
            break;
        }
        if (array(ileft) > limit) {
            double x = double(ileft);
            double y = array(ileft);
            x_vals.push_back(x);
            y_vals.push_back(y);
            mean_x += x;
            mean_y += y;
            nleft++;
        }
        ileft--;
    }

    // Add remaining nodes on the right of the element of interest
    while (nright < 2*nodes-nleft) {
        if (iright >= array.size()) {
            break;
        }
        if (array(iright) > limit) {
            double x = double(iright);
            double y = array(iright);
            x_vals.push_back(x);
            y_vals.push_back(y);
            mean_x += x;
            mean_y += y;
            nright++;
        }
        iright++;
    }

    // Add remaining nodes on the left if all right nodes were not exhausted
    while (nleft < 2*nodes-nright) {
        if (ileft < 0) {
            break;
        }
        if (array(ileft) > limit) {
            double x = double(ileft);
            double y = array(ileft);
            x_vals.push_back(x);
            y_vals.push_back(y);
            mean_x += x;
            mean_y += y;
            nleft++;
        }
        ileft--;
    }

    // Compute mean x and y values
    double total = double(nleft+nright);
    if (total > 0.0) {
        mean_x /= total;
        mean_y /= total;
    }

    // Compute regression line slope
    double beta_nom   = 0.0;
    double beta_denom = 0.0;
    for (int i = 0; i < x_vals.size(); ++i) {
        double x    = x_vals[i];
        double y    = y_vals[i];
        beta_nom   += (x - mean_x) * (y - mean_y);
        beta_denom += (x - mean_x) * (x - mean_x);
    }
    double beta = beta_nom / beta_denom;

    // Compute regression line offset
    double alpha = mean_y - beta * mean_x;

    // Compute value
    double value = alpha + beta * double(inx);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Compute moments
 *
 * @param[in] itheta Offset angle index.
 * @param[out] mean Pointer to mean array.
 * @param[out] rms Pointer to rms array.
 * @param[out] total Pointer to total array.
 *
 * Computes the mean and root mean square values as function of true energy
 * for a given offset angle.
 ***************************************************************************/
void GCTAEdisp2D::get_moments(const int& itheta,
                              GNdarray*  mean,
                              GNdarray*  rms,
                              GNdarray*  total) const
{
    // Get axes dimensions
    int netrue = m_edisp.axis_bins(m_inx_etrue);
    int nmigra = m_edisp.axis_bins(m_inx_migra);
    int npix   = netrue * nmigra;

    // Loop over all true energies
    for (int ietrue = 0; ietrue < netrue; ++ietrue) {

        // Compute base index
        int ibase = itheta * npix + ietrue;

        // Extract array
        GNdarray array(nmigra);
        for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
            array(imigra) = m_edisp(m_inx_matrix, i);
        }

        // Store sum
        (*total)(ietrue) = sum(array);

        // Get Gaussian mean and rms
        double mean_value = 0.0;
        double rms_value  = 0.0;
        get_mean_rms(array, &mean_value, &rms_value);

        // Store mean and rms
        (*mean)(ietrue) = mean_value;
        (*rms)(ietrue)  = rms_value;

    } // endfor: looped over all true energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute mean and root mean square of migration array
 *
 * @param[in] array Energy dispersion array.
 * @param[out] mean Pointer to mean migration value.
 * @param[out] rms Pointer to root mean square value.
 *
 * Computes the mean and the root mean square of the migration array. If the
 * migration array is empty the mean and root mean square values are set to
 * zero.
 *
 * The method does not check whether the @p mean and @p rms pointers are
 * valid.
 ***************************************************************************/
void GCTAEdisp2D::get_mean_rms(const GNdarray& array,
                               double*         mean,
                               double*         rms) const
{
    // Initialise mean and rms
    *mean = 0.0;
    *rms  = 0.0;

    // Get reference to migration values
    const GNodeArray& migras = m_edisp.axis_nodes(m_inx_migra);
    int               nmigra = migras.size();

    // Pre-compute values for mean and rms computation
    double sum = 0.0;
    for (int imigra = 0; imigra < nmigra; ++imigra) {
        double migra   = migras[imigra];
        double weight  = migra * array(imigra);
        *mean         += migra * array(imigra);
        *rms          += weight * migra;
        sum           += array(imigra);
    }

    // If the array is not empty then compute now mean and standard deviation
    if (sum > 0.0) {
        *mean /= sum;
        *rms  /= sum;
        *rms  -= (*mean) * (*mean);
        *rms   = std::sqrt(*rms);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return Gaussian approximation of energy dispersion array
 *
 * @param[in] mean Gaussian mean.
 * @param[in] rms Gaussian rms.
 * @param[in] total Gaussian total.
 * @return Gaussian approximation of energy dispersion array.
 *
 * Returns a Gaussian approximation of the energy dispersion array by
 * computing the mean migration value and its root mean square and by using
 * these values as the centre and the width of a Gaussian function. The
 * Gaussian function is normalized so that the sum of the output array is
 * unity.
 ***************************************************************************/
GNdarray GCTAEdisp2D::gaussian_array(const double& mean,
                                     const double& rms,
                                     const double& total) const
{
    // Initialise empty Gaussian array
    int      nmigra = m_edisp.axis_bins(m_inx_migra);
    GNdarray gaussian(nmigra);

    // If the rms is valid then compute Gaussian
    if (rms > 0.0) {

        // Get reference to migration values
        const GNodeArray& migras = m_edisp.axis_nodes(m_inx_migra);
        int               nmigra = migras.size();

        // Compute Gaussian
        double total_gauss = 0.0;
        for (int imigra = 0; imigra < nmigra; ++imigra) {
            double arg       = (migras[imigra] - mean) / rms;
            double value     = std::exp(-0.5*arg*arg);
            gaussian(imigra) = value;
            total_gauss     += value;
        }

        // Normalise Gaussian
        if (total_gauss > 0.0) {
            gaussian *= total / total_gauss;
        }

    } // endif: computed Gaussian

    // Return Gaussian array
    return gaussian;
}


/***********************************************************************//**
 * @brief Integration kernel for GCTAEdisp2D::edisp_ereco_kern class
 *
 * @param[in] log10Ereco Log10 of reconstructed energy (\f$\log_{10}\f$ MeV).
 * @return Energy dispersion (\f$(\log_{10}\f$ MeV\f$)^{-1}\f$).
 *
 * This method implements the function
 *
 * \f[
 *    E_{\rm disp}(\log_{10} E_{\rm reco} | E_{\rm true}, \theta) =
 *    E_{\rm disp}(m | E_{\rm true}, \theta) \times
 *    \frac{\log 10 \times E_{\rm reco}}{E_{\rm true}}
 * \f]
 *
 * which is the integration kernel needed for the
 * GCTAEdisp2D::edisp_ereco_kern class. The class is used by
 * GCTAEdisp2D::normalize_table() to integrate the energy dispersion
 * information using
 *
 * \f[
 *    \int_{\log_{10} E_{\rm reco}^{\rm min}}^{\log_{10} E_{\rm reco}^{\rm max}}
 *    E_{\rm disp}(\log_{10} E_{\rm reco} | E_{\rm true}, \theta) \,
 *    d(\log_{10} E_{\rm reco})
 * \f]
 ***************************************************************************/
double GCTAEdisp2D::edisp_ereco_kern::eval(const double& log10Ereco)
{
    // Set reconstructued event energy
    GEnergy ereco;
    ereco.log10MeV(log10Ereco);

    // Get function value
    double value = m_parent->operator()(ereco, m_etrue, m_theta);

    // Normalize energy dispersion so that the units are per log10 MeV
    value *= gammalib::ln10 * ereco.MeV();

    // Compile option: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTAEdisp2D::edisp_ereco_kern::eval";
        std::cout << "(log10Ereco=" << log10Ereco << "): ";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}
