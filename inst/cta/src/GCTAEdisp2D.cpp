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
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GCTAEdisp2D.hpp"

//#define G_CDF_METHOD
#define G_REJECTION

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
 * Construct instance by loading the energy dispersion information from a
 * FITS file.
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
 * Returns the energy dispersion in units of s^-1 MeV^-1 for a given observed energy,
 * true photon energy and offset angle.
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
 * @param[in] filename FITS file.
 *
 * This method loads the energy dispersion information from a FITS file.
 ***************************************************************************/
void GCTAEdisp2D::load(const std::string& filename)
{

    // Open FITS file
    GFits fits(filename);

    // Read energy dispersion from file
    read(fits);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy dispersion from FITS file
 *
 * @param[in] fits FITS file pointer.
 *
 * Reads the energy dispersion form the FITS file extension
 * "ENERGY DISPERSION". The data are stored in m_edisp which is of type
 * GCTAResponseTable. The energy axis will be set to log10, the offset
 * angle axis to radians.
 *
 * The method assures that the energy dispersion information is properly
 * normalised.
 ***************************************************************************/
void GCTAEdisp2D::read(const GFits& fits)
{

    // Clear response table
    m_edisp.clear();

    // Get energy dispersion table
    const GFitsTable& table = *fits.table("ENERGY DISPERSION");

    // Read energy dispersion table
    m_edisp.read(table);

    // Get axes dimensions
    int etrue_size = m_edisp.axis(0);
    int migra_size = m_edisp.axis(1);
    int theta_size = m_edisp.axis(2);

    // Normalize vectors of migration matrix for all etrue and theta
    for (int i_etrue = 0; i_etrue < etrue_size; ++i_etrue) {
        for (int i_theta = 0; i_theta < theta_size; ++i_theta) {

            // Initialise sum and offset
            double sum    = 0.0;
            int    offset = i_etrue + (i_theta*migra_size) * etrue_size;

            // Compute sum
            for (int i_migra = 0, i = offset; i_migra < migra_size;
                 ++i_migra, i += etrue_size) {

                // Compute delta(Eobs/Etrue) and Eobs/Etrue
                double delta        = m_edisp.axis_hi(1, i_migra) -
                                      m_edisp.axis_lo(1, i_migra);
                double EobsOverEsrc = 0.5 * (m_edisp.axis_hi(1, i_migra) +
                                             m_edisp.axis_lo(1, i_migra));

                // Add dispersion matrix element to sum
                sum += m_edisp(0, i) * delta / EobsOverEsrc;

            }

            // Normalise by sum/log(10.0)
            if (sum != 0.0) {
                for (int i_migra = 0, i = offset; i_migra < migra_size;
                     ++i_migra, i += etrue_size) {
                    m_edisp(0, i) *= (gammalib::ln10 / sum);
                }
            }

        } // endfor: looped over theta
    } // endfor: looped over Etrue

    // Set true energy axis to logarithmic scale
    m_edisp.axis_log10(0);

    // Set offset angle axis to radians
    m_edisp.axis_radians(2);

    // ...
    m_temp.clear();
    int step = migra_size / etrue_size;
    for (int itrue = 0; itrue < etrue_size; ++itrue) {
        double logEsrc = std::log10(0.5 * (m_edisp.axis_hi(0, itrue) +
                                           m_edisp.axis_lo(0, itrue)));
        m_temp.push_back(std::vector<double>());
        for (int iobs = 0; iobs < etrue_size; ++iobs) {
            double EobsOverEsrc = 0.5 * (m_edisp.axis_hi(1, iobs*step) +
                                         m_edisp.axis_lo(1, iobs*step));
            m_temp[itrue].push_back(operator()(logEsrc, EobsOverEsrc));
        }
    }

    // Set fmax
    set_fmax();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA energy dispersion table into FITS binary table object.
 *
 * @param[in] hdu FITS binary table.
 *
 * @todo Add necessary keywords.
 ***************************************************************************/
void GCTAEdisp2D::write(GFitsBinTable& hdu) const
{
    // Write background table
    m_edisp.write(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy dispersion table into FITS file
 *
 * @param[in] filename Energy dispersion table FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * Save the energy dispersion table into a FITS file.
 ***************************************************************************/
void GCTAEdisp2D::save(const std::string& filename, const bool& clobber) const
{
    // Create binary table
    GFitsBinTable table;
    table.extname("ENERGY DISPERSION");

    // Write the Energy dispersion table
    write(table);

    // Create FITS file, append table, and write into the file
    GFits fits;
    fits.append(table);
    fits.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate energy dispersion
 *
 * @param[in] ran Random number generator.
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 *
 * Draws observed energy value.
 ***************************************************************************/
#if defined(G_REJECTION)
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
    double logEobs = emin;
    double f       = 0.0;
    double ftest   = 1.0;
    while (ftest > f) {
        logEobs = emin + ewidth * ran.uniform();
        f      = operator()(logEobs, logEsrc, theta, phi, zenith, azimuth);
        ftest  = ran.uniform() * m_fmax;
    }

    // Set observed energy
    GEnergy energy;
    energy.log10TeV(logEobs);

    // Return energy
    return energy;

}
#else
GEnergy GCTAEdisp2D::mc(GRan&         ran,
                        const double& logEsrc,
                        const double& theta,
                        const double& phi,
                        const double& zenith,
                        const double& azimuth) const
{
    // Compute cumulative probability
    compute_cumul(theta, phi, zenith, azimuth);

    // Find right logEsrc index with bisection search
    int low = 0;
    int high = m_edisp.axis(0);
    double Esrc = std::exp(logEsrc*gammalib::ln10);
    while ((high-low) > 1) {
        int mid = (low+high) / 2;
        if (Esrc < m_edisp.axis_lo(0, mid)) {
            high = mid;
        }
        else {
            low = mid;
        }
    }
    // Index found
    int isrc = low;

    // Draw random number between 0 and 1 from uniform distribution
    double p = ran.uniform();

    // Find right index with bisection search
    low = 0;
    high = m_cumul[isrc].size();
    while ((high-low) > 1) {
        int mid = (low+high) / 2;
        if (p < m_cumul[isrc][mid].second) {
            high = mid;
        }
        else if (m_cumul[isrc][mid].second <= p) {
            low = mid;
        }
    }

    // Index found
    int index = 0;
    if (low >= m_cumul[isrc].size() - 1) {
        index = m_cumul[isrc].size() - 2;
    }
    else {
        index = low;
    }


    // Interpolate EobsOverEsrc value
    double EobsOverEsrc = 0.0;

    if (m_cumul[isrc][index+1].second == m_cumul[isrc][index].second) {
        EobsOverEsrc = m_cumul[isrc][index].first;
    }
    else {
        EobsOverEsrc  =   (m_cumul[isrc][index+1].second - p)*m_cumul[isrc][index].first
                        + (p - m_cumul[isrc][index].second)*m_cumul[isrc][index+1].first;
        EobsOverEsrc /=   (m_cumul[isrc][index+1].second - m_cumul[isrc][index].second);
    }


    // Set energy
    GEnergy energy;
    energy.TeV(EobsOverEsrc * Esrc);

    // Return energy
    return energy;
}
#endif

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
    if (!m_ebounds_obs_computed || theta != m_theta) {

        m_ebounds_obs_computed = true;
        m_theta = theta;

        // Compute ebounds_obs
        compute_ebounds_obs(theta, phi, zenith, azimuth);

    }

    // Search index only if logEsrc has changed
    if (logEsrc != m_logEsrc) {

        m_logEsrc = logEsrc;

        // Find right index with bisection
        int low = 0;
        int high = m_ebounds_obs.size();
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
    if (!m_ebounds_src_computed || theta != m_theta) {

        m_ebounds_src_computed = true;
        m_theta = theta;

        // Compute ebounds_src
        compute_ebounds_src(theta, phi, zenith, azimuth);

    }

    // Search index only if logEobs has changed
    if (logEobs != m_logEobs) {

        m_logEobs = logEobs;

        // Find right index with bisection
        int low = 0;
        int high = m_ebounds_src.size();

        while ((high-low) > 1) {
            int  mid = (low+high) / 2;

            double e = m_edisp.axis_lo(0, mid);
            //double e = m_eobs_axis[mid];
            if (logEobs < std::log10(e)) {
                high = mid;
            }
            else {
                low = mid;
            }
        }
        // Index found
        m_index_src = low;
        //m_index_src = low * m_edisp.axis(1) / m_edisp.axis(0);
        //std::cout << "emin = " << m_ebounds_src[m_index_src].emin().TeV() << std::endl;
        //std::cout << "emax = " << m_ebounds_src[m_index_src].emax().TeV() << std::endl;

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
        double emax = m_edisp.axis_hi(0,m_edisp.axis(0)-1);

        // Compute offset angle boundaries in deg
        double omin = m_edisp.axis_lo(1,0);
        double omax = m_edisp.axis_hi(1,m_edisp.axis(1)-1);

        // Append header
        result.append("=== GCTAEdisp2D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_edisp.axis(0)));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(m_edisp.axis(1)));
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

    // Initialise Monte Carlo cache
    m_theta                =   0.0;
    m_logEsrc              = -30.0;
    m_logEobs              = -30.0;
    m_fmax                 =  -1.0;
    m_index_obs            =     0;
    m_index_src            =     0;
    m_cdf_computed         = false;
    m_ebounds_obs_computed = false;
    m_ebounds_src_computed = false;

    // Initialize vectors
    m_ebounds_obs.clear();
    m_ebounds_src.clear();
    m_cumul.clear();
    m_eobs_axis.clear();
    m_temp.clear();

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

    // Copy Monte Carlo cache
    m_theta                = edisp.m_theta;
    m_logEsrc              = edisp.m_logEsrc;
    m_logEobs              = edisp.m_logEobs;
    m_fmax                 = edisp.m_fmax;
    m_index_obs            = edisp.m_index_obs;
    m_index_src            = edisp.m_index_src;
    m_cdf_computed         = edisp.m_cdf_computed;
    m_ebounds_obs_computed = edisp.m_ebounds_obs_computed;
    m_ebounds_src_computed = edisp.m_ebounds_src_computed;
    m_ebounds_obs          = edisp.m_ebounds_obs;
    m_ebounds_src          = edisp.m_ebounds_src;
    m_cumul                = edisp.m_cumul;
    m_eobs_axis            = edisp.m_eobs_axis;
    m_temp                 = edisp.m_temp;

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
 ***************************************************************************/
void GCTAEdisp2D::compute_ebounds_obs(const double& theta,
                                      const double& phi,
                                      const double& zenith,
                                      const double& azimuth) const
{
    // Clear ebounds_obs vector
    m_ebounds_obs.clear();

    // Set epsilon
    const double eps = 1.0e-06;

    // Loop over Esrc
    for (int isrc = 0; isrc < m_edisp.axis(0); ++isrc) {

        // Set Esrc
        double Esrc   = 0.5 * (m_edisp.axis_hi(0, isrc) +
                               m_edisp.axis_lo(0, isrc));
        double logEsrc = std::log10(Esrc);

        // Initialise results
        double logEobsMin = -30.0;
        double logEobsMax =  10.0;
        bool   minFound   = false;
        bool   maxFound   = false;

        // Find boundaries, loop over EobsOverEsrc
        for (int i = 0; i < m_edisp.axis(1); ++i) {

            // Compute value EobsOverEsrc
            double EobsOverEsrc = 0.5 * (m_edisp.axis_hi(1, i) +
                                         m_edisp.axis_lo(1, i));

            // Find first non-negligible matrix term
            if (!minFound && m_edisp(0, logEsrc, EobsOverEsrc, theta) >= eps) {
                minFound   = true;
                logEobsMin = std::log10(EobsOverEsrc * Esrc);
            }

            // Find last non-negligible matrix term
            else if (minFound && !maxFound &&
                     m_edisp(0, logEsrc, EobsOverEsrc, theta) < eps) {
                maxFound   = true;
                logEobsMax = std::log10(EobsOverEsrc * Esrc);
            }

            // Continue searching
            else if (minFound && maxFound &&
                     m_edisp(0, logEsrc, EobsOverEsrc, theta) >= eps) {
                maxFound = false;
            }

        } // endfor: looped over EobsOverEsrc


        // If energy dispersion has never become negligible until end of loop,
        // reset logEobsMax
        if (!maxFound) {
            logEobsMax = 10.0;
        }

        // Set energy boundaries
        GEnergy emin;
        GEnergy emax;
        emin.log10TeV(logEobsMin);
        emax.log10TeV(logEobsMax);

        // Add energy boundaries to vector
        m_ebounds_obs.push_back(GEbounds(emin, emax));

    }


    // Return
    return;
}

/***********************************************************************//**
 * @brief Compute ebounds_src vector
 *
 ***************************************************************************/
void GCTAEdisp2D::compute_ebounds_src(const double& theta,
                                      const double& phi,
                                      const double& zenith,
                                      const double& azimuth) const
{

    // Clear ebounds_src vector
    m_ebounds_src.clear();

    // Set epsilon
    const double eps = 1.0e-06;

/*
    // Loop over EobsOverEsrc
    for (int iobs = 0; iobs < m_edisp.axis(1); ++iobs) {

        // Set EobsOverEsrc
        double EobsOverEsrc   = 0.5 * (m_edisp.axis_hi(1, iobs) +
                                       m_edisp.axis_lo(1, iobs));

        // Initialise results
        double logEsrcMin = -10.0;
        double logEsrcMax =  30.0;
        bool   minFound   = false;
        bool   maxFound   = false;

        // Find boundaries, loop over Esrc
        for (int i = 0; i < m_edisp.axis(0); ++i) {

            // Compute value of corresponding logEsrc
            double Esrc     = 0.5 * (m_edisp.axis_hi(0, i) +
                                     m_edisp.axis_lo(0, i));
            double logEsrc  = std::log10(Esrc);

            // Find first non-negligible matrix term
            if (!minFound && m_edisp(0, logEsrc, EobsOverEsrc, theta) >= eps) {
                minFound   = true;
                logEsrcMin = logEsrc;
            }

            // Find last non-negligible matrix term
            else if (minFound && !maxFound &&
                     m_edisp(0, logEsrc, EobsOverEsrc, theta) < eps) {
                maxFound   = true;
                logEsrcMax = logEsrc;
            }

            // Continue
            else if (minFound && maxFound &&
                     m_edisp(0, logEsrc, EobsOverEsrc, theta) >= eps) {
                maxFound = false;
            }

        } // endfor: looped over true energy

        // If energy dispersion has never become negligible until end of loop,
        // reset logEobsMax
        if (!maxFound) {
            logEsrcMax = 30.0;
        }

        // Set energy boundaries
        GEnergy emin;
        GEnergy emax;
        emin.log10TeV(logEsrcMin);
        emax.log10TeV(logEsrcMax);

        // Add energy boundaries to vector
        m_ebounds_src.push_back(GEbounds(emin, emax));

    }

    const int step = m_edisp.axis(1) / m_edisp.axis(0);

    // Loop over EobsOverEsrc to create Eobs axis
    for (int i = 0; i < m_edisp.axis(0); ++i) {
        m_eobs_axis.push_back(m_edisp.axis_lo(1, i*step) * m_edisp.axis_lo(0, i));
    }
*/

    // Loop over Eobs
    for (int iobs = 0; iobs < m_edisp.axis(0); ++iobs) {

        double Eobs    = 0.5 * (m_edisp.axis_hi(0, iobs) +
                                m_edisp.axis_lo(0, iobs));
        double logEobs = std::log10(Eobs);

        // Initialise results
        double logEsrcMin = -10.0;
        double logEsrcMax =  30.0;
        bool   minFound   = false;
        bool   maxFound   = false;

        // Find boundaries, loop over Esrc
        for (int isrc = 0; isrc < m_edisp.axis(0); ++isrc) {

            double Esrc = 0.5 * (m_edisp.axis_hi(0, isrc) +
                                 m_edisp.axis_lo(0, isrc));
            double logEsrc = std::log10(Esrc);

            // Find first non-negligible matrix term
            if (!minFound && operator()(logEobs, logEsrc, theta) >= eps) {
                minFound   = true;
                logEsrcMin = logEsrc;
            }

            // Find last non-negligible matrix term
            else if (minFound && !maxFound &&
                     operator()(logEobs, logEsrc, theta) < eps) {
                maxFound   = true;
                logEsrcMax = logEsrc;
            }

            // Continue
            else if (minFound && maxFound &&
                     operator()(logEobs, logEsrc, theta) >= eps) {
                maxFound = false;
            }

        } // endfor: looped over true energy

        // If energy dispersion has never become negligible until end of loop,
        // reset logEobsMax
        if (!maxFound) {
            logEsrcMax = 30.0;
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
 * @brief Compute cumulative probability
 *
 ***************************************************************************/
#if defined(G_CDF_METHOD)
void GCTAEdisp2D::compute_cumul(const double& theta,
                                const double& phi,
                                const double& zenith,
                                const double& azimuth) const
{
    // TODO ? compute also for each theta ??

    if (!m_cdf_computed) {

        m_cumul.clear();

        //m_theta = theta;
        m_cdf_computed = true;

    // Loop over Esrc
    for (int isrc = 0; isrc < m_edisp.axis(0); ++isrc) {

        m_cumul.push_back(std::vector<std::pair<double, double> >());

        //double Esrc = 0.5 * (m_edisp.axis_lo(0, isrc) + m_edisp.axis_hi(0, isrc));
        double Esrc = m_edisp.axis_lo(0, isrc);

        // Initialize cumulative probability
        double sum = 0.0;

        // Divide bin into n sub-bins
        const int n = 10;

        // Loop over EobsOverEtrue
        for (int imigra = 0; imigra < m_edisp.axis(1); ++imigra) {

            // Compute delta(Eobs/Etrue) and Eobs/Etrue values
            double delta        = m_edisp.axis_hi(1, imigra) -
                                  m_edisp.axis_lo(1, imigra);
            double EobsOverEsrcmin = m_edisp.axis_lo(1, imigra);

            for (int i = 0; i < n; ++i) {

                double EobsOverEsrc = EobsOverEsrcmin + i*delta/n;

                // Compute cumulative probability
                double add = m_edisp(0, std::log10(Esrc), EobsOverEsrc, theta) *
                             delta / (n * EobsOverEsrc * gammalib::ln10);

                // Add to sum (and keep sum lower or equal to 1)
                sum = sum+add >= 1.0 ? 1.0 : sum+add;

                // Create pair containing EobsOverEsrc and cumulated probability
                std::pair<double, double> pair(EobsOverEsrc, sum);

                // Add to vector
                m_cumul[isrc].push_back(pair);

            }

        } // endfor: looped over EobsOverEsrc

    }

    }

    // Return
    return;
}
#endif


/***********************************************************************//**
 * @brief Find maximum energy dispersion value
 ***************************************************************************/
void GCTAEdisp2D::set_fmax(void) const
{
    // Initialise maximum
    m_fmax = 0.0;

    // Loop over all response table elements
    for (int i = 0; i < m_edisp.elements(); ++i) {
        double value = m_edisp(0,i);
        if (value > m_fmax) {
            m_fmax = value;
        }
    }

    // Return
    return;
}
