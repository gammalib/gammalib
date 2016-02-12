/***************************************************************************
 *            GCTAEdisp2D.cpp - CTA 2D energy dispersion class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Florent Forest                                   *
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
#include "GIntegral.hpp"
#include "GFilename.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GCTAEdisp2D.hpp"

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
GCTAEdisp2D::GCTAEdisp2D(const std::string& filename) : GCTAEdisp()
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
    // Initalize edisp
    double edisp = 0.0;

    // Compute Eobs/Esrc
    double EobsOverEsrc = std::exp((logEobs-logEsrc) * gammalib::ln10);

    // Compute edisp
    edisp = m_edisp(0, logEsrc, EobsOverEsrc, theta);

    // Return
    return edisp;
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
void GCTAEdisp2D::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAEdisp::free_members();

    // Initialise members
    this->GCTAEdisp::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of energy dispersion instance.
 ***************************************************************************/
GCTAEdisp2D* GCTAEdisp2D::clone(void) const
{
    return new GCTAEdisp2D(*this);
}


/***********************************************************************//**
 * @brief Load energy dispersion from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the energy dispersion from a FITS file.
 *
 * If no extension name is provided, the energy dispersion will be loaded
 * from the "ENERGY DISPERSION" extension.
 ***************************************************************************/
void GCTAEdisp2D::load(const std::string& filename)
{
    // Create file name
    GFilename fname(filename);

    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(fname.filename());

    // Get energy dispersion table
    const GFitsTable& table = *file.table(fname.extname("ENERGY DISPERSION"));

    // Read energy dispersion from table
    read(table);

    // Close FITS file
    file.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy dispersion from FITS table
 *
 * @param[in] table FITS table.
 *
 * Reads the energy dispersion form the FITS @p table.
 *
 * The data are stored in m_edisp which is of type GCTAResponseTable. The
 * energy axis will be set to log10, the offset angle axis to radians.
 *
 * The method assures that the energy dispersion is properly normalised.
 ***************************************************************************/
void GCTAEdisp2D::read(const GFitsTable& table)
{
    // Clear response table
    m_edisp.clear();

    // Read energy dispersion table
    m_edisp.read(table);

    // Set true energy axis to logarithmic scale
    m_edisp.axis_log10(0);

    // Set offset angle axis to radians
    m_edisp.axis_radians(2);

    // Get axes dimensions
    int etrue_size = m_edisp.axis_bins(0);
    int migra_size = m_edisp.axis_bins(1);
    int theta_size = m_edisp.axis_bins(2);

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
            double theta = 0.5 * (m_edisp.axis_lo(2,i_theta) +
                                  m_edisp.axis_hi(2,i_theta)) *
                           gammalib::deg2rad;

            // Loop over true photon energy
            for (int i_etrue = 0; i_etrue < etrue_size; ++i_etrue) {

                // Get energy
                double emin    = std::log10(m_edisp.axis_lo(0,i_etrue));
                double emax    = std::log10(m_edisp.axis_hi(0,i_etrue));
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
                        m_edisp(0,i) /= sum;
                    }
                }
            }
        }
    
    } // endfor: made 2 passes

    // Set maximum energy dispersion value
    set_max_edisp();

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
    // Write background table
    m_edisp.write(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy dispersion table into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file? (default: false)
 *
 * Save the energy dispersion table into FITS file.
 *
 * If no extension name is provided, the energy dispersion will be saved into
 * the "ENERGY DISPERSION" extension.
 ***************************************************************************/
void GCTAEdisp2D::save(const std::string& filename, const bool& clobber) const
{
    // Create file name
    GFilename fname(filename);

    // Create binary table
    GFitsBinTable table;
    table.extname(fname.extname("ENERGY DISPERSION"));

    // Write the energy dispersion table
    write(table);

    // Create FITS file, append table, and write into the file
    GFits fits;
    fits.append(table);
    fits.saveto(fname.filename(), clobber);

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
 *
 * Draws observed energy value given a true energy @p logEsrc and offset
 * angle @p theta.
 ***************************************************************************/
GEnergy GCTAEdisp2D::mc(GRan&         ran,
                        const double& logEsrc,
                        const double& theta,
                        const double& phi,
                        const double& zenith,
                        const double& azimuth) const
{
    // Get boundaries for observed energy
    GEbounds ebounds = ebounds_obs(logEsrc, theta, phi, zenith, azimuth);
    double   emin    = ebounds.emin().log10TeV();
    double   emax    = ebounds.emax().log10TeV();

    // Find energy by rejection method
    double ewidth  = emax - emin;
    double logEobs = 0.5*(emin+emax);
    if (m_max_edisp > 0.0) {
        double f       = 0.0;
        double ftest   = 1.0;
        while (ftest > f) {
            logEobs = emin + ewidth * ran.uniform();
            f       = operator()(logEobs, logEsrc, theta, phi, zenith, azimuth);
            ftest   = ran.uniform() * m_max_edisp;
        }
    }

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
    // Compute only if parameters changed
    if (!m_ebounds_obs_computed || theta != m_last_theta_obs) {

        // Set computation flag
        m_ebounds_obs_computed = true;
        m_last_theta_obs       = theta;

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
            double e = m_edisp.axis_lo(0, mid);
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
    // Compute only if parameters changed
    if (!m_ebounds_src_computed || theta != m_last_theta_src) {

        // Set computation flag
        m_ebounds_src_computed = true;
        m_last_theta_src       = theta;

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
            double e = m_edisp.axis_lo(0, mid);
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
 * @brief Print energy dispersion information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy dispersion information.
 ***************************************************************************/
std::string GCTAEdisp2D::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        double emin = m_edisp.axis_lo(0,0);
        double emax = m_edisp.axis_hi(0,m_edisp.axis_bins(0)-1);

        // Compute offset angle boundaries in deg
        double omin = m_edisp.axis_lo(1,0);
        double omax = m_edisp.axis_hi(1,m_edisp.axis_bins(1)-1);

        // Append header
        result.append("=== GCTAEdisp2D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_edisp.axis_bins(0)));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(m_edisp.axis_bins(1)));
        result.append("\n"+gammalib::parformat("Log10(Energy) range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");
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
    m_filename  = edisp.m_filename;
    m_edisp     = edisp.m_edisp;

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

    // Loop over Esrc
    for (int isrc = 0; isrc < m_edisp.axis_bins(0); ++isrc) {

        // Set Esrc
        double Esrc    = std::sqrt(m_edisp.axis_hi(0,isrc) *
                                   m_edisp.axis_lo(0,isrc));
        double logEsrc = std::log10(Esrc);

        // Initialise results
        double logEobsMin = 0.0;
        double logEobsMax = 0.0;
        bool   minFound   = false;
        bool   maxFound   = false;

        // Determine number of MIGRA bins
        int n_migra = m_edisp.axis_bins(1);

        // Find minimum boundary
        for (int i = 0; i < n_migra; ++i) {

            // Compute EobsOverEsrc value
            double EobsOverEsrc = 0.5 * (m_edisp.axis_hi(1,i) +
                                         m_edisp.axis_lo(1,i));

            // Get matrix term
            double edisp = m_edisp(0, logEsrc, EobsOverEsrc, theta);

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
            double EobsOverEsrc = 0.5 * (m_edisp.axis_hi(1,i) +
                                         m_edisp.axis_lo(1,i));

            // Get matrix term
            double edisp = m_edisp(0, logEsrc, EobsOverEsrc, theta);

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

    // Loop over Eobs
    for (int iobs = 0; iobs < m_edisp.axis(0); ++iobs) {

        // Set Eobs
        double Eobs    = std::sqrt(m_edisp.axis_hi(0,iobs) *
                                   m_edisp.axis_lo(0,iobs));
        double logEobs = std::log10(Eobs);

        // Initialise results
        double logEsrcMin = 0.0;
        double logEsrcMax = 0.0;
        bool   minFound   = false;
        bool   maxFound   = false;

        // Determine number of true energy bins
        int n_bins = m_edisp.axis(0);

        // Find minimum boundary
        for (int isrc = 0; isrc < n_bins; ++isrc) {

            // Set Esrc
            double Esrc    = std::sqrt(m_edisp.axis_hi(0,isrc) *
                                       m_edisp.axis_lo(0,isrc));
            double logEsrc = std::log10(Esrc);

            // Get matrix term
            double edisp = operator()(logEobs, logEsrc, theta);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                minFound   = true;
                logEsrcMin = logEsrc;
                break;
            }

        } // endfor: find minimum boundary

        // Find maximum boundary
        for (int isrc = n_bins-1; isrc >= 0; isrc--) {

            // Set Esrc
            double Esrc    = std::sqrt(m_edisp.axis_hi(0,isrc) *
                                       m_edisp.axis_lo(0,isrc));
            double logEsrc = std::log10(Esrc);

            // Get matrix term
            double edisp = operator()(logEobs, logEsrc, theta);

            // Find first non-negligible matrix term
            if (edisp >= eps) {
                maxFound   = true;
                logEsrcMax = logEsrc;
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
 * @brief Set maximum energy dispersion value
 ***************************************************************************/
void GCTAEdisp2D::set_max_edisp(void) const
{
    // Initialise maximum
    m_max_edisp = 0.0;

    // Loop over all response table elements
    for (int i = 0; i < m_edisp.elements(); ++i) {
        double value = m_edisp(0,i);
        if (value > m_max_edisp) {
            m_max_edisp = value;
        }
    }

    // Return
    return;
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
