/***************************************************************************
 *            GCTAEdisp2D.cpp - CTA 2D energy dispersion class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2018 by Florent Forest                              *
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

/* __ Method name definitions ____________________________________________ */
#define G_READ                               "GCTAEdisp2D::read(GFitsTable&)"
#define G_MC    "GCTAEdisp2D::mc(GRan&, double&, double&, double&, double&, "\
                                                                   "double&)"
#define G_FETCH                                        "GCTAEdisp2D::fetch()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_SMOOTH_EDISP_KLUDGE    //!< Get rid of noise in energy dispersion
#define G_SMOOTH_EDISP_SAVE_TEST //!< Save test matrix

/* __ Debug definitions __________________________________________________ */

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
 * @brief Return energy dispersion in units of s^-1 MeV^-1
 *
 * @param[in] logEobs Log10 of the measured energy (TeV).
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Defaults to 0.0.
 * @param[in] phi Azimuth angle in camera system (rad). Not used in this method.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used in this method.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used in this method.
 *
 * Returns the energy dispersion in units of s^-1 MeV^-1 for a given observed
 * energy @p logEobs, true photon energy @p logEsrc and offset angle
 * @p theta.
 ***************************************************************************/
double GCTAEdisp2D::operator()(const double& logEobs, 
                               const double& logEsrc, 
                               const double& theta, 
                               const double& phi,
                               const double& zenith,
                               const double& azimuth) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Initalize edisp
    double edisp = 0.0;

    // Continue only if logEsrc and theta are in validity range
    if ((logEsrc >= m_logEsrc_min)  && (logEsrc <= m_logEsrc_max) &&
        (theta   >= m_theta_min)    && (theta   <= m_theta_max)) {

        // Compute Eobs/Esrc
        double EobsOverEsrc = std::exp((logEobs-logEsrc) * gammalib::ln10);

        // Continue only if migration is in validity range
        if ((EobsOverEsrc >= m_migra_min)  && (EobsOverEsrc <= m_migra_max)) {

            // Setup argument vector
            double arg[3];
            arg[m_inx_etrue] = logEsrc;
            arg[m_inx_migra] = EobsOverEsrc;
            arg[m_inx_theta] = theta;

            // Compute edisp
            edisp = m_edisp(m_inx_matrix, arg[0], arg[1], arg[2]);

            // Make sure that energy dispersion is non-negative
            if (edisp < 0.0) {
                edisp = 0.0;
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
 *     ETRUE_LO - True energy lower bin boundaries
 *     ETRUE_HI - True energy upper bin boundaries
 *     MIGRA_LO - Migration lower bin boundaries
 *     MIGRA_HI - Migration upper bin boundaries
 *     THETA_LO - Offset angle lower bin boundaries
 *     THETA_HI - Offset angle upper bin boundaries
 *     MATRIX   - Migration matrix
 *
 * The data are stored in the m_edisp member. The energy axis will be set to
 * log10, the offset angle axis to radians.
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
 * If no extension name is provided, the energy dispersion will be loaded
 * from the `ENERGY DISPERSION` extension.
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
 * @param[in] clobber Overwrite existing file? (default: false)
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
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad).
 * @param[in] phi Azimuth angle in camera system (rad).
 * @param[in] zenith Zenith angle in Earth system (rad).
 * @param[in] azimuth Azimuth angle in Earth system (rad).
 * @return Observed energy.
 *
 * @exception GException::invalid_return_value
 *            No energy dispersion information found for parameters or
 *            energy dispersion matrix is invalid.
 *
 * Draws observed energy value given a true energy @p logEsrc and offset
 * angle @p theta. If no energy dispersion information is available the
 * method will return the true photon energy.
 ***************************************************************************/
GEnergy GCTAEdisp2D::mc(GRan&         ran,
                        const double& logEsrc,
                        const double& theta,
                        const double& phi,
                        const double& zenith,
                        const double& azimuth) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Get boundaries for observed energy
    GEbounds ebounds = ebounds_obs(logEsrc, theta, phi, zenith, azimuth);
    double   emin    = ebounds.emin().log10TeV();
    double   emax    = ebounds.emax().log10TeV();

    // Throw an exception if minimum energy is equal or larger than maximum
    // energy (this means that no energy dispersion information is available
    // for that energy)
    if (emin >= emax) {
        double      eng = std::pow(10.0, logEsrc);
        std::string msg = "No valid energy dispersion information available "
                          "for true photon energy "+gammalib::str(eng)+" TeV,"
                          "offset angle "+
                          gammalib::str(theta*gammalib::rad2deg)+" deg "
                          "and azimuth angle "+
                          gammalib::str(phi*gammalib::rad2deg)+" deg. Either "
                          "provide an energy dispersion matrix covering these "
                          "parameters or restrict the parameter space.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Throw an exception if maximum energy dispersion is zero
    if (m_max_edisp <= 0.0) {
        std::string msg = "Energy dispersion matrix is empty. Please provide "
                          "a valid energy dispersion matrix.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Initialise rejection method
    double    ewidth               = emax - emin;
    double    logEobs              = logEsrc;
    double    f                    = 0.0;
    double    ftest                = 1.0;
    int       zeros                = 0;
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
            logEobs = emin + ewidth * ran.uniform();
            f       = operator()(logEobs, logEsrc, theta, phi, zenith, azimuth);
            if (f == 0.0) {
                zeros++;
                if (zeros > max_subsequent_zeros) {
                    double      eng = std::pow(10.0, logEsrc);
                    std::string msg = "No valid energy dispersion information "
                                      "available  for true photon energy "+
                                      gammalib::str(eng)+" TeV, offset angle "+
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
        // dispersion value. The loop is quite if the random number is smaller
        // than the energy dispersion value
        ftest = ran.uniform() * m_max_edisp;

    } // endwhile:

    // Set observed energy
    GEnergy energy;
    energy.log10TeV(logEobs);

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Return observed energy interval that contains the energy dispersion.
 *
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad).
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @return Observed energy boundaries.
 *
 * Returns the band of observed energies outside of which the energy
 * dispersion becomes negligible for a given true energy @p logEsrc and
 * offset angle @p theta. An energy is considered negligible if inferior
 * to 1e-6.
 ***************************************************************************/
GEbounds GCTAEdisp2D::ebounds_obs(const double& logEsrc,
                                  const double& theta,
                                  const double& phi,
                                  const double& zenith,
                                  const double& azimuth) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Compute only if parameters changed
    if (!m_ebounds_obs_computed || theta != m_last_theta_obs) {

        // Set computation flag
        m_ebounds_obs_computed = true;
        m_last_theta_obs       = theta;
        m_last_logEsrc         = -30.0; // force update

        // Compute ebounds_obs
        compute_ebounds_obs(theta, phi, zenith, azimuth);

    }

    // Search index only if logEsrc has changed
    if (logEsrc != m_last_logEsrc) {

        // Store last log(Esrc)
        m_last_logEsrc = logEsrc;

        // Find right index with bisection
        int low  = 0;
        int high = m_ebounds_obs.size() - 1;
        while ((high-low) > 1) {
            int  mid = (low+high) / 2;
            double e = m_edisp.axis_lo(m_inx_etrue, mid);
            if (logEsrc < std::log10(e)) {
                high = mid;
            }
            else {
                low = mid;
            }
        }

        // Index found
        m_index_obs = low;

    }

    // Return energy boundaries
    return m_ebounds_obs[m_index_obs];
}


/***********************************************************************//**
 * @brief Return true energy interval that contains the energy dispersion.
 *
 * @param[in] logEobs Log10 of the observed event energy (TeV).
 * @param[in] theta Offset angle in camera system (rad).
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @return True energy boundaries.
 *
 * Returns the band of true photon energies outside of which the energy
 * dispersion becomes negligible for a given observed energy @p logEobs and
 * offset angle @p theta. An energy in considered negligible if inferior to
 * 1e-6.
 ***************************************************************************/
GEbounds GCTAEdisp2D::ebounds_src(const double& logEobs,
                                  const double& theta,
                                  const double& phi,
                                  const double& zenith,
                                  const double& azimuth) const
{
    // Make sure that energy dispersion is online
    fetch();

    // Compute true energy boundaries only if parameters changed
    if (!m_ebounds_src_computed || theta != m_last_theta_src) {

        // Set computation flag
        m_ebounds_src_computed = true;
        m_last_theta_src       = theta;
        m_last_logEobs         = -30.0; // force update

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
 * @brief Fetch energy dispersion
 *
 * @exception GException::file_error
 *            File not found.
 *            Unable to load energy dispersion.
 * @exception GException::invalid_value
 *            No file name has been specified.
 *
 * Fetches the energy dispersion by reading it from a FITS file. This method
 * does nothing if the energy dispersion is already loaded, if there is
 * nothing to fetch, or if the m_filename member is empty.
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

                // Initialise energy dispersion extension name
                std::string extname =
                            m_filename.extname(gammalib::extname_cta_edisp2d);

                // Get energy dispersion table
                const GFitsTable& table = *fits.table(extname);

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

        // Throw an exception if the FITS file name is not known
        else {
            std::string msg = "Unable to fetch energy dispersion since no "
                              "filename is specified.";
            throw GException::invalid_value(G_FETCH, msg);
        }
    
    } // endif: energy dispersion had not yet been fetched
    
    // Return
    return;
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

        // Make sure that energy dispersion is online
        fetch();

        // Compute energy boundaries in TeV
        double emin = m_edisp.axis_lo(m_inx_etrue,0);
        double emax = m_edisp.axis_hi(m_inx_etrue,
                                      m_edisp.axis_bins(m_inx_etrue)-1);

        // Compute mingration
        double mmin = m_edisp.axis_lo(m_inx_migra,0);
        double mmax = m_edisp.axis_hi(m_inx_migra,
                                      m_edisp.axis_bins(m_inx_migra)-1);

        // Compute offset angle boundaries in deg
        double omin = m_edisp.axis_lo(m_inx_theta,0);
        double omax = m_edisp.axis_hi(m_inx_theta,
                                      m_edisp.axis_bins(m_inx_theta)-1);

        // Append header
        result.append("=== GCTAEdisp2D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_edisp.axis_bins(m_inx_etrue)));
        result.append("\n"+gammalib::parformat("Number of migration bins") +
                      gammalib::str(m_edisp.axis_bins(m_inx_migra)));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(m_edisp.axis_bins(m_inx_theta)));
        result.append("\n"+gammalib::parformat("Log10(Energy) range"));
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

    // Initialise computation cache
    m_ebounds_obs_computed = false;
    m_ebounds_src_computed = false;
    m_last_theta_obs       =  -1.0;
    m_last_theta_src       =  -1.0;
    m_last_logEsrc         = -30.0;
    m_last_logEobs         = -30.0;
    m_index_obs            =     0;
    m_index_src            =     0;
    m_max_edisp            =   0.0;
    m_ebounds_obs.clear();
    m_ebounds_src.clear();

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

    // Copy computation cache
    m_ebounds_obs_computed = edisp.m_ebounds_obs_computed;
    m_ebounds_src_computed = edisp.m_ebounds_src_computed;
    m_last_theta_obs       = edisp.m_last_theta_obs;
    m_last_theta_src       = edisp.m_last_theta_src;
    m_last_logEsrc         = edisp.m_last_logEsrc;
    m_last_logEobs         = edisp.m_last_logEobs;
    m_index_obs            = edisp.m_index_obs;
    m_index_src            = edisp.m_index_src;
    m_max_edisp            = edisp.m_max_edisp;
    m_ebounds_obs          = edisp.m_ebounds_obs;
    m_ebounds_src          = edisp.m_ebounds_src;

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
 * @brief Compute ebounds_obs vector
 *
 * @param[in] theta Offset angle (rad).
 * @param[in] phi Azimuth angle (rad).
 * @param[in] zenith Zenith angle (rad).
 * @param[in] azimuth Azimuth angle (rad).
 *
 * Computes for all true energies the energy boundaries of the observed
 * energies covered by valid migration matrix elements. Only matrix elements
 * with values >= 1.0e-12 are considered as valid elements. In case that no
 * matrix elements are found of a given true energy, the interval of observed
 * energies will be set to [1 TeV, 1 TeV] (i.e. an empty interval).
 ***************************************************************************/
void GCTAEdisp2D::compute_ebounds_obs(const double& theta,
                                      const double& phi,
                                      const double& zenith,
                                      const double& azimuth) const
{
    // Clear ebounds_obs vector
    m_ebounds_obs.clear();

    // Set epsilon
    const double eps = 1.0e-12;

    // Determine number of true energy and migration bins
    int n_bins  = m_edisp.axis_bins(m_inx_etrue);
    int n_migra = m_edisp.axis_bins(m_inx_migra);

    // Loop over Esrc
    for (int isrc = 0; isrc < n_bins; ++isrc) {

        // Set Esrc
        double Esrc    = std::sqrt(m_edisp.axis_hi(m_inx_etrue,isrc) *
                                   m_edisp.axis_lo(m_inx_etrue,isrc));
        double logEsrc = std::log10(Esrc);

        // Initialise results
        double logEobsMin = 0.0;
        double logEobsMax = 0.0;
        bool   minFound   = false;
        bool   maxFound   = false;

        // Find minimum boundary
        for (int i = 0; i < n_migra; ++i) {

            // Compute EobsOverEsrc value
            double EobsOverEsrc = 0.5 * (m_edisp.axis_hi(m_inx_migra,i) +
                                         m_edisp.axis_lo(m_inx_migra,i));

            // Get matrix term
            double arg[3];
            arg[m_inx_etrue] = logEsrc;
            arg[m_inx_migra] = EobsOverEsrc;
            arg[m_inx_theta] = theta;
            double edisp     = m_edisp(m_inx_matrix, arg[0], arg[1], arg[2]);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                minFound   = true;
                logEobsMin = std::log10(EobsOverEsrc * Esrc);
                break;
            }

        } // endfor: find minimum boundary

        // Find maximum boundary
        for (int i = n_migra-1; i >= 0; i--) {

            // Compute EobsOverEsrc value
            double EobsOverEsrc = 0.5 * (m_edisp.axis_hi(m_inx_migra,i) +
                                         m_edisp.axis_lo(m_inx_migra,i));

            // Get matrix term
            double arg[3];
            arg[m_inx_etrue] = logEsrc;
            arg[m_inx_migra] = EobsOverEsrc;
            arg[m_inx_theta] = theta;
            double edisp     = m_edisp(m_inx_matrix, arg[0], arg[1], arg[2]);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                maxFound   = true;
                logEobsMax = std::log10(EobsOverEsrc * Esrc);
                break;
            }

        } // endfor: find maximum boundary

        // If we did not find boundaries then set the interval to
        // a zero interval for safety
        if (!minFound || !maxFound) {
            logEobsMin = 0.0;
            logEobsMax = 0.0;
        }

        // Set energy boundaries
        GEnergy emin;
        GEnergy emax;
        emin.log10TeV(logEobsMin);
        emax.log10TeV(logEobsMax);

        // Add energy boundaries to vector
        m_ebounds_obs.push_back(GEbounds(emin, emax));

    } // endfor: looped over logEsrc

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute ebounds_src vector
 *
 * @param[in] theta Offset angle (rad).
 * @param[in] phi Azimuth angle (rad).
 * @param[in] zenith Zenith angle (rad).
 * @param[in] azimuth Azimuth angle (rad).
 *
 * Computes for all observed energies the energy boundaries of the true
 * energies covered by valid migration matrix elements. Only matrix elements
 * with values >= 1.0e-12 are considered as valid elements. In case that no
 * matrix elements are found of a given observed energy, the interval of true
 * energies will be set to [1 TeV, 1 TeV] (i.e. an empty interval).
 ***************************************************************************/
void GCTAEdisp2D::compute_ebounds_src(const double& theta,
                                      const double& phi,
                                      const double& zenith,
                                      const double& azimuth) const
{
    // Clear ebounds_src vector
    m_ebounds_src.clear();

    // Set epsilon
    const double eps = 1.0e-12;

    // Determine number of true energy bins
    int n_bins   = m_edisp.axis_bins(m_inx_etrue);
    int n_migras = m_edisp.axis_bins(m_inx_migra);

    // Loop over Eobs
    for (int iobs = 0; iobs < n_bins; ++iobs) {

        // Set Eobs
        double Eobs    = std::sqrt(m_edisp.axis_hi(m_inx_etrue,iobs) *
                                   m_edisp.axis_lo(m_inx_etrue,iobs));
        double logEobs = std::log10(Eobs);

        // Initialise results
        double logEsrcMin = 0.0;
        double logEsrcMax = 0.0;
        bool   minFound   = false;
        bool   maxFound   = false;

        // Find minimum boundary
        for (int imigra = 0; imigra < n_migras; ++imigra) {

            // Get mean migra
            double migra = 0.5 * (m_edisp.axis_lo(m_inx_migra,imigra) +
                                  m_edisp.axis_lo(m_inx_migra,imigra));

            // Fall through if migration is zero
            if (migra == 0.0) {
                continue;
            }

            // Compute log10 of true energy
            double logEsrc = logEobs - std::log10(migra);

            // Get matrix term
            double edisp = operator()(logEobs, logEsrc, theta);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                maxFound   = true;
                logEsrcMax = logEsrc;
                break;
            }

        } // endfor: find minimum boundary

        // Find maximum boundary
        for (int imigra = n_migras-1; imigra >= 0; imigra--) {

            // Get mean migra
            double migra = 0.5 * (m_edisp.axis_lo(m_inx_migra,imigra) +
                                  m_edisp.axis_lo(m_inx_migra,imigra));

            // Fall through if migration is zero
            if (migra == 0.0) {
                continue;
            }

            // Compute log10 of true energy
            double logEsrc = logEobs - std::log10(migra);

            // Get matrix term
            double edisp = operator()(logEobs, logEsrc, theta);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                minFound   = true;
                logEsrcMin = logEsrc;
                break;
            }

        } // endfor: find maximum boundary

        // If we did not find boundaries then set the interval to
        // a zero interval for safety
        if (!minFound || !maxFound) {
            logEsrcMin = 0.0;
            logEsrcMax = 0.0;
        }

        // Set energy boundaries
        GEnergy emin;
        GEnergy emax;
        emin.log10TeV(logEsrcMin);
        emax.log10TeV(logEsrcMax);

        // Add energy boundaries to vector
        m_ebounds_src.push_back(GEbounds(emin, emax));

    } // endfor : loopend over observed energy

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set table
 *
 * After assigning or loading a response table, performs all necessary steps
 * to set the table. The method also smoothes the table get rid of numerical
 * fluctuations and normalises the table properly.
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
    m_inx_etrue  = m_edisp.axis("ETRUE");
    m_inx_migra  = m_edisp.axis("MIGRA");
    m_inx_theta  = m_edisp.axis("THETA");
    m_inx_matrix = m_edisp.table("MATRIX");

    // Set true energy axis to logarithmic scale
    m_edisp.axis_log10(m_inx_etrue);

    // Set offset angle axis to radians
    m_edisp.axis_radians(m_inx_theta);

    // Set table boundaries
    set_boundaries();

    // Smooth energy dispersion table
    #if defined(G_SMOOTH_EDISP_KLUDGE)
    denoise_table();
    clean_table(15);
    clip_table(0.001);
    /*
    denoise_table();
    clean_table(20);
    smooth_table(5.0);
    clip_table(0.01);
    */
    #endif

    // Normalize energy dispersion table
    normalize_table();

    // Set maximum energy dispersion value
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
 * @brief Set maximum energy dispersion value
 ***************************************************************************/
void GCTAEdisp2D::set_max_edisp(void) const
{
    // Initialise maximum
    m_max_edisp = 0.0;

    // Loop over all response table elements
    for (int i = 0; i < m_edisp.elements(); ++i) {
        double value = m_edisp(m_inx_matrix,i);
        if (value > m_max_edisp) {
            m_max_edisp = value;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Normalize energy dispersion table
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
            double theta = 0.5 * (m_edisp.axis_lo(m_inx_theta,i_theta) +
                                  m_edisp.axis_hi(m_inx_theta,i_theta)) *
                           gammalib::deg2rad;

            // Loop over true photon energy
            for (int i_etrue = 0; i_etrue < etrue_size; ++i_etrue) {

                // Get energy
                double emin    = std::log10(m_edisp.axis_lo(m_inx_etrue,i_etrue));
                double emax    = std::log10(m_edisp.axis_hi(m_inx_etrue,i_etrue));
                double logEsrc = 0.5*(emin+emax);

                // Initialise integration
                double sum = 0.0;

                // Get integration boundaries
                GEbounds ebounds = ebounds_obs(logEsrc, theta);

                // Loop over all energy intervals
                for (int i = 0; i < ebounds.size(); ++i) {

                    // Get energy boundaries
                    emin = ebounds.emin(i).log10TeV();
                    emax = ebounds.emax(i).log10TeV();

                    // Setup integration function
                    edisp_kern integrand(this, logEsrc, theta);
                    GIntegral  integral(&integrand);

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
                    int offset = i_etrue + (i_theta*migra_size) * etrue_size;
                    for (int k = 0, i = offset; k < migra_size; ++k, i += etrue_size) {
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
 * @brief Clip energy dispersion table
 *
 * @param[in] threshold Minimum number of consecutive non-zeros
 *
 * Clips noise in the energy dispersion table.
 *
 * For each offset angle and true energy, the method determines the maximum
 * value of the energy dispersion table and then walks towards smaller and
 * larger migration value to determine the first pixel for which the energy
 * dispersion becomes zero. All pixels beyond these pixels are then clipped
 * to zero. This removes noise that is at small or large migration values.
 *
 * The method then computes the number of remaining contiguos non-zero pixels
 * for each offset angle and true energy, and if this number is smaller than
 * a threshold, all pixels are set to zero. This removes true energy bins
 * that are sparesely filled.
 ***************************************************************************/
void GCTAEdisp2D::clean_table(const int& threshold)
{
    // Get axes dimensions
    int netrue = m_edisp.axis_bins(m_inx_etrue);
    int nmigra = m_edisp.axis_bins(m_inx_migra);
    int ntheta = m_edisp.axis_bins(m_inx_theta);

    // Loop over all offset angles
    for (int itheta = 0; itheta < ntheta; ++itheta) {

        // Loop ober all true energies
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {

            // Compute base index
            int ibase = ietrue + itheta * netrue * nmigra;

            // Determine index of maximum value
            int    imax = 0;
            double max  = m_edisp(m_inx_matrix,ibase+ietrue);
            for (int imigra = 1; imigra < nmigra; ++imigra) {
                int    inx   = ibase + imigra*netrue;
                double value = m_edisp(m_inx_matrix, inx);
                if (value > max) {
                    imax = imigra;
                    max  = value;
                }
            }

            // Do nothing if all pixels are zero
            if (max == 0.0) {
                continue;
            }

            // Walk from maximum to the first pixel
            bool clip = false;
            for (int imigra = imax; imigra >= 0; imigra--) {
                int inx = ibase + imigra*netrue;
                if (clip) {
                    m_edisp(m_inx_matrix, inx) = 0.0;
                }
                else {
                    double value = m_edisp(m_inx_matrix, inx);
                    if (value == 0.0) {
                        clip = true;
                    }
                }
            }

            // Walk from maximum to the last pixel
            clip = false;
            for (int imigra = imax; imigra < nmigra; ++imigra) {
                int inx = ibase + imigra*netrue;
                if (clip) {
                    m_edisp(m_inx_matrix, inx) = 0.0;
                }
                else {
                    double value = m_edisp(m_inx_matrix, inx);
                    if (value == 0.0) {
                        clip = true;
                    }
                }
            }

            // Determine number of contiguous non-zero pixels
            int nonzero            = 0;
            int contiguous_nonzero = 0;
            for (int imigra = 0; imigra < nmigra; ++imigra) {
                int    inx   = ibase + imigra*netrue;
                double value = m_edisp(m_inx_matrix, inx);
                if (value != 0.0) {
                    nonzero++;
                    if (nonzero > contiguous_nonzero) {
                        contiguous_nonzero = nonzero;
                    }
                }
                else {
                    nonzero = 0;
                }
            }

            // If there are less than threshold contiguous non-zero pixels we
            // only have noise, hence all pixels are set to zero
            if (contiguous_nonzero < threshold) {
                for (int imigra = 0; imigra < nmigra; ++imigra) {
                    int inx = ibase + imigra*netrue;
                    m_edisp(m_inx_matrix, inx) = 0.0;
                }
            }

        } // endfor: looped over all true energies

    } // endfor: looped over all offset angles

    // Return
    return;
}


/***********************************************************************//**
 * @brief Smooth energy dispersion table
 *
 * @param[in] sigma Smoothing width (number of bins in true energy)
 *
 * Smooth the energy dispersion table in true energy. The smoothing is done
 * using a fast-fourrier transform. For this purpose the energy disperison
 * table for each offset angle is copied into a 2D array, padded with some
 * zero pixels on the left and the right to avoid wrap around.
 ***************************************************************************/
void GCTAEdisp2D::smooth_table(const double& sigma)
{
    // Get axes dimensions
    int netrue        = m_edisp.axis_bins(m_inx_etrue);
    int nmigra        = m_edisp.axis_bins(m_inx_migra);
    int ntheta        = m_edisp.axis_bins(m_inx_theta);
    int netrue_pad    = int(3.0*sigma);
    int netrue_padded = netrue + 2*netrue_pad; // Array padded with zeros
    int npix          = netrue * nmigra;

    // Get smoothing kernel
    GFft fft_kernel = fft_smooth_kernel(netrue_padded, nmigra, sigma);

    // Loop over all offset angles
    for (int itheta = 0; itheta < ntheta; ++itheta) {

        // Compute base index
        int ibase = itheta * npix;

        // Extract matrix in GNdarray
        GNdarray array(netrue_padded, nmigra);
        for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra) {
            double* ptr = array.data() + imigra * netrue_padded + netrue_pad;
            for (int ietrue = 0; ietrue < netrue; ++ietrue, ++i) {
                *ptr++ = m_edisp(m_inx_matrix, i);
            }
        }

        // FFT of array
        GFft fft_array(array);

        // Smooth array
        GFft fft_smooth = fft_array * fft_kernel;

        // Backward transform array
        GNdarray smooth = fft_smooth.backward();

        // Put back array values in matrix
        for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra) {
            double* ptr = smooth.data() + imigra * netrue_padded + netrue_pad;
            for (int ietrue = 0; ietrue < netrue; ++ietrue, ++i) {
                m_edisp(m_inx_matrix, i) = *ptr++;
            }
        }

    } // endfor: looped over all offset angles

    // Return
    return;
}


/***********************************************************************//**
 * @brief Smooth energy dispersion table
 *
 * @param[in] sigma Smoothing width (number of bins in true energy)
 *
 * Smooth the energy dispersion table in true energy. The smoothing is done
 * using a fast-fourrier transform. For this purpose the energy disperison
 * table for each offset angle is copied into a 2D array, padded with some
 * zero pixels on the left and the right to avoid wrap around.
 ***************************************************************************/
void GCTAEdisp2D::smooth_table2(const double& sigma_etrue,
                                const double& sigma_migra)
{
    // Get axes dimensions
    int netrue        = m_edisp.axis_bins(m_inx_etrue);
    int nmigra        = m_edisp.axis_bins(m_inx_migra);
    int ntheta        = m_edisp.axis_bins(m_inx_theta);
    int netrue_pad    = int(3.0*sigma_etrue);
    int netrue_padded = netrue + 2*netrue_pad; // Array padded with zeros
    int npix          = netrue * nmigra;

    // Get smoothing kernel
    GFft fft_kernel = fft_smooth_kernel2(netrue_padded, nmigra, sigma_etrue, sigma_migra);

    // Loop over all offset angles
    for (int itheta = 0; itheta < ntheta; ++itheta) {

        // Compute base index
        int ibase = itheta * npix;

        // Extract matrix in GNdarray
        GNdarray array(netrue_padded, nmigra);
        for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra) {
            double* ptr = array.data() + imigra * netrue_padded + netrue_pad;
            for (int ietrue = 0; ietrue < netrue; ++ietrue, ++i) {
                *ptr++ = m_edisp(m_inx_matrix, i);
            }
        }

        // FFT of array
        GFft fft_array(array);

        // Smooth array
        GFft fft_smooth = fft_array * fft_kernel;

        // Backward transform array
        GNdarray smooth = fft_smooth.backward();

        // Put back array values in matrix
        for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra) {
            double* ptr = smooth.data() + imigra * netrue_padded + netrue_pad;
            for (int ietrue = 0; ietrue < netrue; ++ietrue, ++i) {
                m_edisp(m_inx_matrix, i) = *ptr++;
            }
        }

    } // endfor: looped over all offset angles

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get FFT of smoothing kernel
 *
 * @param[in] nbins Number of kernel bins.
 * @param[in] sigma Gaussian kernel sigma in bins.
 * @return FFT of smoothing kernel
 *
 * Returns the fast-fourrier transform of a Gaussian smoothing kernel.
 ***************************************************************************/
GFft GCTAEdisp2D::fft_smooth_kernel(const int&    nbins,
                                    const double& sigma) const
{
    // Allocate kernel
    GNdarray kernel(nbins);

    // Initialise sum and compute Gaussian normalisation
    double sum  =  0.0;
    double norm = -0.5 / (sigma * sigma);

    // Set Gaussian kernel
    for (int i = 0; i < nbins; ++i) {
        double value = std::exp(norm*double(i*i));
        kernel(i) += value;
        sum       += value;
        if (i > 0) {
            kernel(nbins-i) += value;
            sum             += value;
        }
    }

    // Normalize kernel
    if (sum > 0.0) {
        for (int i = 0; i < nbins; ++i) {
            kernel(i) /= sum;
        }
    }

    // Return FFT of kernel
    return (GFft(kernel));
}


/***********************************************************************//**
 * @brief Get FFT of smoothing kernel for true energy smoothing
 *
 * @param[in] netrue Number of true energy bins.
 * @param[in] nmigra Number of migration bins bins.
 * @param[in] sigma Gaussian kernel sigma in true energy direction.
 * @return FFT of smoothing kernel
 *
 * Returns the fast-fourrier transform of a Gaussian smoothing kernel in
 * true energy.
 ***************************************************************************/
GFft GCTAEdisp2D::fft_smooth_kernel(const int&    netrue,
                                    const int&    nmigra,
                                    const double& sigma) const
{
    // Allocate kernel
    GNdarray kernel(netrue, nmigra);

    // Initialise sum and compute Gaussian normalisation
    double sum  =  0.0;
    double norm = -0.5 / (sigma * sigma);

    // Set Gaussian kernel
    for (int i = 0; i < netrue; ++i) {
        double value = std::exp(norm*double(i*i));
        kernel(i,0) += value;
        sum         += value;
        if (i > 0) {
            kernel(netrue-i,0) += value;
            sum                += value;
        }
    }

    // Normalize kernel
    if (sum > 0.0) {
        for (int i = 0; i < netrue; ++i) {
            kernel(i,0) /= sum;
        }
    }

    // Return FFT of kernel
    return (GFft(kernel));
}


/***********************************************************************//**
 * @brief Get FFT of smoothing kernel for 2D smoothing
 *
 * @param[in] netrue Number of true energy bins.
 * @param[in] nmigra Number of migration bins bins.
 * @param[in] sigma_etrue Gaussian kernel sigma in true energy direction.
 * @param[in] sigma_migra Gaussian kernel sigma in migration direction.
 * @return FFT of smoothing kernel
 *
 * Returns the fast-fourrier transform of a Gaussian smoothing kernel in
 * true energy and migration.
 ***************************************************************************/
GFft GCTAEdisp2D::fft_smooth_kernel2(const int&    netrue,
                                     const int&    nmigra,
                                     const double& sigma_etrue,
                                     const double& sigma_migra) const
{
    // Allocate kernel
    GNdarray kernel(netrue, nmigra);

    // Initialise sum and compute Gaussian normalisation
    double sum        =  0.0;
    double norm_etrue = -0.5 / (sigma_etrue * sigma_etrue);
    double norm_migra = -0.5 / (sigma_migra * sigma_migra);

    // Set Gaussian kernel
    for (int itrue = 0; itrue < netrue; ++itrue) {
        double value_true = std::exp(norm_etrue*double(itrue*itrue));
        for (int imigra = 0; imigra < nmigra; ++imigra) {
            double value = std::exp(norm_migra*double(imigra*imigra)) * value_true;
            kernel(itrue,imigra) += value;
            sum                  += value;
            if (itrue > 0) {
                kernel(netrue-itrue,imigra) += value;
                sum                         += value;
            }
            if (imigra > 0) {
                kernel(itrue,nmigra-imigra) += value;
                sum                         += value;
            }
            if ((itrue > 0) && (imigra > 0)) {
                kernel(netrue-itrue,nmigra-imigra) += value;
                sum                                += value;
            }
        }
    }

    // Normalize kernel
    if (sum > 0.0) {
        for (int itrue = 0; itrue < netrue; ++itrue) {
            for (int imigra = 0; imigra < nmigra; ++imigra) {
                kernel(itrue,imigra) /= sum;
            }
        }
    }

    // Return FFT of kernel
    return (GFft(kernel));
}


/***********************************************************************//**
 * @brief Clip energy dispersion table
 *
 * @param[in] threshold Clipping threshold in fraction of maximum value.
 *
 * Clips all migration vectors at a threshold that is a fraction of the
 * maximum value in the migration vector.
 ***************************************************************************/
void GCTAEdisp2D::clip_table(const double& threshold)
{
    // Get axes dimensions
    int netrue = m_edisp.axis_bins(m_inx_etrue);
    int nmigra = m_edisp.axis_bins(m_inx_migra);
    int ntheta = m_edisp.axis_bins(m_inx_theta);
    int npix   = netrue * nmigra;

    // Loop over all offset angles
    for (int itheta = 0; itheta < ntheta; ++itheta) {

        // Loop over all true energies
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {

            // Compute base index
            int ibase = itheta * npix + ietrue;

            // Compute maximum value
            double max = 0.0;
            for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
                if (m_edisp(m_inx_matrix, i) > max) {
                    max = m_edisp(m_inx_matrix, i);
                }
            }

            // Clip values below threshold
            double clip_value = max * threshold;
            for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
                if (m_edisp(m_inx_matrix, i) < clip_value) {
                    m_edisp(m_inx_matrix, i) = 0.0;
                }
            }

        } // endfor: looped over true energies

    } // endfor: looped over all offset angles

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clip array
 *
 * @param[in] array Ndarray that should be clipped.
 * @param[in] threshold Clipping threshold in fraction of maximum value.
 * @return Clipped Ndarray
 *
 * Clips all values below a threshold to zero.
 ***************************************************************************/
/*
GNdarray GCTAEdisp2D::clip_array(const GNdarray& array,
                                 const double&   threshold) const
{
    // Set clipped array
    GNdarray clipped = array;

    // Initialise variables
    int nsize = clipped.size();

    // Determine maximum value in array
    double  max  = 0.0;
    double* cptr = clipped.data();
    for (int i = 0; i < nsize; ++i, cptr++) {
        if (*cptr > max) {
            max = *cptr;
        }
    }

    // Clip values below threshold
    double clip_value = max * threshold;
    cptr              = clipped.data();
    for (int i = 0; i < nsize; ++i, cptr++) {
        if (std::abs(*cptr) < clip_value) {
            *cptr = 0.0;
        }
    }

    // Return clipped array
    return clipped;
}
*/


/***********************************************************************//**
 * @brief Denoise energy dispersion table
 *
 * Denoise the energy dispersion table by smoothing the energy dispersion
 * for each offset angle and true energy using the smooth_array() method.
 *
 * True energy values with no energy dispersion information are interpolated
 * from neighbouring energy dispersion values.
 ***************************************************************************/
void GCTAEdisp2D::denoise_table(void)
{
    // Log entrance
    //std::cout << "GCTAEdisp2D::denoise_table in" << std::endl;

    // Get axes dimensions
    int netrue = m_edisp.axis_bins(m_inx_etrue);
    int nmigra = m_edisp.axis_bins(m_inx_migra);
    int ntheta = m_edisp.axis_bins(m_inx_theta);
    int npix   = netrue * nmigra;

    // Get Gaussian parameters for all offset angles and true energies
    GNdarray gpars = get_gaussian_pars(5.0);

    // Loop over all offset angles
    for (int itheta = 0; itheta < ntheta; ++itheta) {

        // Initialise sums
        std::vector<double> sum_array;

        // Loop over all true energies
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {

            // Compute base index
            int ibase = itheta * npix + ietrue;

            // Extract array
            GNdarray array(nmigra);
            for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
                array(imigra) = m_edisp(m_inx_matrix, i);
            }

            // Compute single event value
            double event_value = get_single_event_value(array);

            // Smooth array
            GNdarray smoothed_array(nmigra);
            if (max(array) >= 2.0*event_value) {
                smoothed_array = gaussian_array(gpars(itheta, ietrue, 0),
                                                gpars(itheta, ietrue, 1));
            }

            // Compute total of smoothed array
            double total = sum(smoothed_array);

            // Store sum
            sum_array.push_back(total);

            // Put back smoothed array values in matrix
            if (total > 0.0) {
                for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
                    m_edisp(m_inx_matrix, i) = smoothed_array(imigra) / total;
                }
            }

        } // endfor: looped over all true energies

        // Interpolate gaps
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {

            // Fall through if array is not empty
            if (sum_array[ietrue] > 0.0) {
                continue;
            }

            // Find enclosing indices of non-empty arrays
            int ileft  = ietrue;
            int iright = ietrue;
            for (; ileft >= 0; ileft--) {
                if (sum_array[ileft] > 0.0) {
                    break;
                }
            }
            for (; iright < netrue; ++iright) {
                if (sum_array[iright] > 0.0) {
                    break;
                }
            }

            // If encosing indices are valid then interpolate array
            if ((ileft >= 0) && (iright < netrue)) {

                // Compute weight factors
                double norm   = 1.0 / double(iright-ileft);
                double wleft  = double(iright - ietrue) * norm;
                double wright = double(ietrue - ileft)  * norm;

                // Fill array by linear interpolation
                int idst       = itheta * npix + ietrue;
                int isrc_left  = itheta * npix + ileft;
                int isrc_right = itheta * npix + iright;
                for (int imigra = 0; imigra < nmigra; ++imigra) {
                    m_edisp(m_inx_matrix, idst) =
                            wleft  * m_edisp(m_inx_matrix, isrc_left) +
                            wright * m_edisp(m_inx_matrix, isrc_right);
                    idst       += netrue;
                    isrc_left  += netrue;
                    isrc_right += netrue;
                }

            } // endif: there were valid enclosing indices
        
        } // endfor: looped over all true energies

    } // endfor: looped over all offset angles

    // Log exit
    //std::cout << "GCTAEdisp2D::denoise_table out" << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Smooth array
 *
 * @param[in] array Energy dispersion array.
 * @param[in] mean Mean of Gaussian.
 * @param[in] rms Root mean square of Gaussian.
 * @return Smoothed energy dispersion array.
 ***************************************************************************/
GNdarray GCTAEdisp2D::smooth_array(const GNdarray& array,
                                   const double&   mean,
                                   const double&   rms) const
{
    // Get reference to migration values
    const GNodeArray& migras = m_edisp.axis_nodes(m_inx_migra);
    int               nmigra = migras.size();

    // Initialise empty smoothed array
    GNdarray smoothed_array(nmigra);

    // Set required number of events and minimum and maximum sigma
    double num_events = 200.0;          // required number of events under Gaussian
    double min_events =   2.0;          // at least 2 events are required
    double sigma_min  =   1.0;          // minimum sigma
    double sigma_max  =   0.1 * nmigra; // maximum sigma (< 20)
    if (sigma_max > 20.0) {
       sigma_max = 20.0;
    }

    // Compute single event value
    double event_value = get_single_event_value(array);

    // Continue only there are pixels with the requested minimum number of
    // events (otherwise the returned vector will be zero)
    if (max(array) >= min_events*event_value) {

        // Get Gaussian approximation
        GNdarray gaussian = gaussian_array(mean, rms);

        // Get residual by subtracting Gaussian from array
        GNdarray residual = array - gaussian;

        // Initialise smoothed residual
        GNdarray smoothed_residual(nmigra);

        // Loop over all migration pixels
        for (int imigra = 0; imigra < nmigra; ++imigra) {

            // Get pixel value
            double value = array(imigra);

            // Fall through if pixel value is empty
            if (value == 0.0) {
                continue;
            }

            // Set Gaussian sigma. Make sure that sigma is within the allowed
            // range
            double sigma = num_events * event_value / value;
            if (sigma < sigma_min) {
                sigma = sigma_min;
            }
            else if (sigma > sigma_max) {
                sigma = sigma_max;
            }

            // Convolve residual with a Gaussian function
            double   total =  0.0;
            double   norm  = -0.5 / (sigma * sigma);
            GNdarray work(nmigra);
            for (int k = 0; k < nmigra; ++k) {
                double weight = std::exp(norm*k*k);
                if (weight > 1.0e-30) {
                    if (k == 0) {
                        work(imigra) = weight;
                        total       += weight;
                    }
                    else {
                        int ileft = imigra - k;
                        if (ileft >= 0) {
                            work(ileft) = weight;
                            total      += weight;
                        }
                        int iright = imigra + k;
                        if (iright < nmigra) {
                            work(iright) = weight;
                            total       += weight;
                        }
                    }
                }
            }
            if (total > 0.0) {
                smoothed_residual += work * residual(imigra) / total;
            }

        } // endfor: looped over all migration pixels

        // Set smoothed array
        smoothed_array = gaussian; // + smoothed_residual;

        // Make sure that smooth array is not negative
        for (int imigra = 0; imigra < nmigra; ++imigra) {
            if (smoothed_array(imigra) < 0.0) {
                smoothed_array(imigra) = 0.0;
            }
        }

    } // endif: single event value was positive

    // Return smoothed array
    return smoothed_array;
}


/***********************************************************************//**
 * @brief Return Gaussian approximation of energy dispersion array
 *
 * @param[in] mean Gaussian mean.
 * @param[in] rms Gaussian rms.
 * @return Gaussian approximation of energy dispersion array.
 *
 * Returns a Gaussian approximation of the energy dispersion array by
 * computing the mean migration value and its root mean square and by using
 * these values as the centre and the width of a Gaussian function. The
 * Gaussian function is normalized so that the sum of the output array is
 * unity.
 ***************************************************************************/
GNdarray GCTAEdisp2D::gaussian_array(const double& mean, const double& rms) const
{
    // Initialise empty Gaussian array
    int      nmigra = m_edisp.axis_bins(m_inx_migra);
    GNdarray gaussian(nmigra);

    // Get mean and rms of migration values in array
    /*
    double mean = 0.0;
    double rms  = 0.0;
    get_mean_rms(array, &mean, &rms);
    */

    // If the array contains information then compute a Gaussian
    // approximation
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
            gaussian *= 1.0 / total_gauss;
        }

    } // endif: computed Gaussian approximation

    // Return Gaussian array
    return gaussian;
}


/***********************************************************************//**
 * @brief Compute mean and root mean square of migration array
 *
 * @param[in] sigma Smoothing parameter in true energy.
 * @return Array of mean and rms values for each offset angle and true energy.
 *
 * Computes the mean and the root mean square of the migration array for each
 * offset angle and true energy and applies a Gaussian smoothing in true
 * energy.
 ***************************************************************************/
GNdarray GCTAEdisp2D::get_gaussian_pars(const double& sigma) const
{
    // Get axes dimensions
    int netrue        = m_edisp.axis_bins(m_inx_etrue);
    int nmigra        = m_edisp.axis_bins(m_inx_migra);
    int ntheta        = m_edisp.axis_bins(m_inx_theta);
    int npix          = netrue * nmigra;

    // Initialise result
    GNdarray result(ntheta, netrue, 2);

    // Setup smoothing kernel
    GFft fft_kernel = fft_smooth_kernel(netrue, sigma);

    // Loop over all offset angles
    for (int itheta = 0; itheta < ntheta; ++itheta) {

        // Allocate working arrays
        GNdarray work_mean(netrue);
        GNdarray work_rms(netrue);
        GNdarray work_sum(netrue);

        // Initialise mean and rms values
        double mean_value = 0.0;
        double rms_value  = 0.0;

        // Loop over all true energies
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {

            // Compute base index
            int ibase = itheta * npix + ietrue;

            // Extract array
            GNdarray array(nmigra);
            for (int imigra = 0, i = ibase; imigra < nmigra; ++imigra, i += netrue) {
                array(imigra) = m_edisp(m_inx_matrix, i);
            }

            // Store sum (only results with positive sum will be kept)
            work_sum(ietrue) = sum(array);

            // Compute single event value
            double event_value = get_single_event_value(array);

            // Get Gaussian mean and rms. Use last mean and rms value in case
            // that the array is empty
            if (work_sum(ietrue) > 0) {
                mean_value = 0.0;
                rms_value  = 0.0;
                get_mean_rms(array, &mean_value, &rms_value);
            }

            // Store mean and rms
            work_mean(ietrue) = mean_value;
            work_rms(ietrue)  = rms_value;

        } // endfor: looped over all true energies

        // Pad front
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {
            if (work_mean(ietrue) > 0.0) {
                for (int i = 0; i < ietrue; ++i) {
                    work_mean(i) = work_mean(ietrue);
                    work_rms(i)  = work_rms(ietrue);
                }
                break;
            }
        }

        // Pad back
        for (int ietrue = netrue-1; ietrue >= 0; ietrue--) {
            if (work_mean(ietrue) > 0.0) {
                for (int i = ietrue; i < netrue; ++i) {
                    work_mean(i) = work_mean(ietrue);
                    work_rms(i)  = work_rms(ietrue);
                }
                break;
            }
        }

        // Smooth mean
        GFft fft_mean(work_mean);
        GFft fft_mean_smooth = fft_mean * fft_kernel;
        work_mean = fft_mean_smooth.backward();

        // Smooth rms
        GFft fft_rms(work_rms);
        GFft fft_rms_smooth = fft_rms * fft_kernel;
        work_rms = fft_rms_smooth.backward();

        // Store result
        for (int ietrue = 0; ietrue < netrue; ++ietrue) {
            if (work_sum(ietrue) > 0.0) {
                result(itheta, ietrue, 0) = work_mean(ietrue);
                result(itheta, ietrue, 1) = work_rms(ietrue);
            }
        }

    } // endif: looped over all offset angles

    // Return result
    return result;
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
 * @brief Estimate the value of a single event
 *
 * @param[in] array Energy dispersion array.
 * @return Value of a single event.
 *
 * Returns the minimum non-zero value in an array as an estimate of the value
 * of a single event. If the array is empty the method returns zero.
 ***************************************************************************/
double GCTAEdisp2D::get_single_event_value(const GNdarray& array) const
{
    // Initialise minimum value
    double event_value = 1.0e30;

    // Get minimum non-zero value
    for (int i = 0; i < array.size(); ++i) {
        double value = array(i);
        if ((value > 0.0) && (value < event_value)) {
            event_value = value;
        }
    }

    // If no minimum value was encountered then set single event value to
    // zero
    if (event_value == 1.0e30) {
        event_value = 0.0;
    }

    // Return
    return event_value;
}


/***********************************************************************//**
 * @brief Integration kernel for edisp_kern() class
 *
 * @param[in] x Function value.
 *
 * This method implements the integration kernel needed for the edisp_kern()
 * class.
 ***************************************************************************/
double GCTAEdisp2D::edisp_kern::eval(const double& x)
{
    // Get function value
    double value = m_parent->operator()(x, m_logEsrc, m_theta);

    // Compile option: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTAEdisp2D::edisp_kern::eval";
        std::cout << "(x=" << x << "): ";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}
