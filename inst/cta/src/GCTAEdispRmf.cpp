/***************************************************************************
 *            GCTAEdispRmf.cpp - CTA RMF energy dispersion class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2018 by Christoph Deil & Ellis Owen                 *
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
 * @file GCTAEdispRmf.cpp
 * @brief CTA RMF energy dispersion class implementation
 * @author Christoph Deil & Ellis Owen
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GFilename.hpp"
#include "GRan.hpp"
#include "GIntegral.hpp"
#include "GFunction.hpp"
#include "GEnergy.hpp"
#include "GRmf.hpp"
#include "GNodeArray.hpp"
#include "GMatrixSparse.hpp"
#include "GCTAEdispRmf.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                               "GCTAEdispRmf::load(GFilename&)"

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
GCTAEdispRmf::GCTAEdispRmf(void) : GCTAEdisp()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Redistribution Matrix File name.
 *
 * Construct instance by loading the Redistribution Matrix File.
 ***************************************************************************/
GCTAEdispRmf::GCTAEdispRmf(const GFilename& filename) : GCTAEdisp()
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
 * @param[in] edisp Energy dispersion
 ***************************************************************************/
GCTAEdispRmf::GCTAEdispRmf(const GCTAEdispRmf& edisp) : GCTAEdisp(edisp)
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
GCTAEdispRmf::~GCTAEdispRmf(void)
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
 * @param[in] edisp Energy dispersion
 * @return Energy dispersion
 ***************************************************************************/
GCTAEdispRmf& GCTAEdispRmf::operator=(const GCTAEdispRmf& edisp)
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
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Energy dispersion (MeV\f$^{-1}\f$)
 *
 * Returns the energy dispersion
 *
 * \f[
 *    E_{\rm disp}(E_{\rm true}, E_{\rm reco}) =
 *    \frac{E_{\rm disp}(E_{\rm true}, \log_{10} E_{\rm reco}, \theta)}
 *         {\log_{10} E_{\rm reco}}
 * \f]
 *
 * in units of MeV\f$^{-1}\f$ where
 * \f$E_{\rm reco}\f$ is the reconstructed energy in units of MeV, and
 * \f$E_{\rm true}\f$ is the true energy in units of MeV.
 ***************************************************************************/
double GCTAEdispRmf::operator()(const GEnergy& ereco,
                                const GEnergy& etrue,
                                const double&  theta,
                                const double&  phi,
                                const double&  zenith,
                                const double&  azimuth) const
{
    // Get log10 of true and reconstructed photon energies
    double logEsrc = etrue.log10TeV();
    double logEobs = ereco.log10TeV();

    // Update indexes and weighting factors for interpolation
    update(logEobs, logEsrc);

    // Perform interpolation
    double edisp =  m_wgt1 * m_matrix(m_itrue1, m_imeas1) +
                    m_wgt2 * m_matrix(m_itrue1, m_imeas2) +
                    m_wgt3 * m_matrix(m_itrue2, m_imeas1) +
                    m_wgt4 * m_matrix(m_itrue2, m_imeas2);

    // Make sure that energy dispersion is not negative
    if (edisp < 0.0) {
        edisp = 0.0;
    }

    // If the energy dispersion is positive, return it per MeV
    else {
        edisp /= (gammalib::ln10 * ereco.MeV());
    }

    // Return energy dispersion
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
void GCTAEdispRmf::clear(void)
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
 * @return Deep copy of instance.
 ***************************************************************************/
GCTAEdispRmf* GCTAEdispRmf::clone(void) const
{
    return new GCTAEdispRmf(*this);
}


/***********************************************************************//**
 * @brief Load energy dispersion from RMF file
 *
 * @param[in] filename of RMF file.
 *
 * This method loads the energy dispersion information from an RMF file.
 ***************************************************************************/
void GCTAEdispRmf::load(const GFilename& filename)
{
    // Load RMF file
    m_rmf.load(filename);

    // Store the filename
    m_filename = filename;

    // Set cache (has to come before set_matrix())
    set_cache();

    // Set matrix
    set_matrix();

    // Set maximum edisp value
    set_max_edisp();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate energy dispersion
 *
 * @param[in] ran Random number generator.
 * @param[in] etrue True photon energy.
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Energy.
 *
 * Draws observed energy value from RMF matrix.
 ***************************************************************************/
GEnergy GCTAEdispRmf::mc(GRan&          ran,
                         const GEnergy& etrue,
                         const double&  theta,
                         const double&  phi,
                         const double&  zenith,
                         const double&  azimuth) const
{
    // Initialise energies
    GEnergy ereco;

    // Get boundaries for rejection method
    GEbounds ebounds = ebounds_obs(etrue.log10TeV(), theta, phi, zenith, azimuth);
    double   emin    = ebounds.emin().log10TeV();
    double   emax    = ebounds.emax().log10TeV();

    // Find energy by rejection method
    double ewidth  = emax - emin;
    double logEobs = 0.5*(emin+emax);
    if (m_max_edisp > 0.0) {
        double f      = 0.0;
        double ftest  = 1.0;
        while (ftest > f) {
            logEobs = emin + ewidth * ran.uniform();
            ereco.log10TeV(logEobs);
            f       = operator()(ereco, etrue, theta, phi, zenith, azimuth);
            ftest   = ran.uniform() * m_max_edisp;
        }
    }

    // Return reconstructed photon energy
    return ereco;
}


/***********************************************************************//**
 * @brief Return observed energy interval that contains the energy dispersion.
 *
 * @param[in] logEsrc Log10 of the true photon energy (\f$\log_{10}\f$ TeV).
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Reconstructed energy boundaries.
 *
 * Returns the band of observed energies outside of which the energy
 * dispersion becomes negligible for a given true energy @p logEsrc.
 ***************************************************************************/
GEbounds GCTAEdispRmf::ebounds_obs(const double& logEsrc,
                                   const double& theta,
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth) const
{
    // Compute only if parameters changed
    if (!m_ebounds_obs_computed || theta != m_last_theta_obs) {

        // Store last theta value
        m_ebounds_obs_computed = true;
        m_last_theta_obs       = theta;

        // Compute ebounds_obs
        compute_ebounds_obs(theta, phi, zenith, azimuth);

    }

    // Search index only if logEsrc has changed
    if (logEsrc != m_last_logEsrc) {

        // Store last photon energy
        m_last_logEsrc = logEsrc;

        // Find right index with bisection
        int low  = 0;
        int high = m_ebounds_obs.size() - 1;
        while ((high-low) > 1) {
            int  mid = (low+high) / 2;
            double e = m_rmf.etrue().emin(mid).log10TeV();
            if (logEsrc < e) {
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
    return (m_ebounds_obs[m_index_obs]);
}


/***********************************************************************//**
 * @brief Return true energy interval that contains the energy dispersion.
 *
 * @param[in] logEobs Log10 of the observed event energy (\f$\log_{10}\f$ TeV).
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return True energy boundaries.
 *
 * Returns the band of true photon energies outside of which the energy
 * dispersion becomes negligible for a given observed energy @p logEobs.
 ***************************************************************************/
GEbounds GCTAEdispRmf::ebounds_src(const double& logEobs,
                                   const double& theta,
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth) const
{
    // Compute only if parameters changed
    if (!m_ebounds_src_computed || theta != m_last_theta_src) {

        // Store last value
        m_ebounds_src_computed = true;
        m_last_theta_src       = theta;

        // Compute ebounds_obs
        compute_ebounds_src(theta, phi, zenith, azimuth);

    }

    // Search index only if logEobs has changed
    if (logEobs != m_last_logEobs) {

        // Store last value
        m_last_logEobs = logEobs;

        // Find right index with bisection
        int low  = 0;
        int high = m_ebounds_src.size() - 1;
        while ((high-low) > 1) {
            int  mid = (low+high) / 2;
            double e = m_rmf.emeasured().emin(mid).log10TeV();
            if (logEobs < e) {
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
    return (m_ebounds_src[m_index_src]);
}


/***********************************************************************//**
 * @brief Return energy dispersion probability for reconstructed energy
 *        interval
 *
 * @param[in] ereco_min Minimum of reconstructed energy interval.
 * @param[in] ereco_max Maximum of reconstructed energy interval.
 * @param[in] etrue True energy.
 * @param[in] theta Offset angle. Not used.
 * @return Integrated energy dispersion probability.
 *
 * Computes
 *
 * \f[
 *    \int_{E_{\rm reco}^{\rm min}}^{E_{\rm reco}^{\rm max}}
 *    E_{\rm disp}(E_{\rm true}, E_{\rm reco}) \, dE_{\rm reco}
 * \f]
 *
 * where
 * \f$E_{\rm reco}\f$ is the reconstructed energy and
 * \f$E_{\rm true}\f$ is the true energy.
 ***************************************************************************/
double GCTAEdispRmf::prob_erecobin(const GEnergy& ereco_min,
                                   const GEnergy& ereco_max,
                                   const GEnergy& etrue,
                                   const double&  theta) const
{
    // Initalize probability
    double prob = 0.0;

    // Get log10 of energies
    double logEsrc     = etrue.log10TeV();
    double logEobs_min = ereco_min.log10TeV();
    double logEobs_max = ereco_max.log10TeV();

    // Set logEsrc for node array interpolation
    m_etrue.set_value(logEsrc);

    // Store indices and weights for interpolation
    int    inx_left  = m_etrue.inx_left();
    int    inx_right = m_etrue.inx_right();
    double wgt_left  = m_etrue.wgt_left();
    double wgt_right = m_etrue.wgt_right();

    // Initialise first and second node
    GEnergy x1 = ereco_min;
    double  f1 = this->operator()(x1, etrue);
    GEnergy x2;
    double  f2 = 0.0;

    // Loop over all measured energy nodes
    for (int i = 0; i < m_emeasured.size(); ++i) {

        // If measured energy is below logEobs_min then skip node
        if (m_emeasured[i] <= logEobs_min) {
            continue;
        }

        // If measured energy is above maximum migration value then use
        // logEobs_max as measured energy
        if (m_emeasured[i] > logEobs_max) {
            x2 = ereco_max;
            f2 = this->operator()(x2, etrue);
        }

        // ... otherwise use measured energy
        else {
            x2.log10TeV(m_emeasured[i]);
            f2 = wgt_left  * m_matrix(inx_left,  i) +
                 wgt_right * m_matrix(inx_right, i);
            if (f2 < 0.0) {
                f2 = 0.0;
            }
        }

        // Compute integral
        prob += 0.5 * (f1 + f2) * (x2.log10TeV() - x1.log10TeV());

        // Set second node as first node
        x1 = x2;
        f1 = f2;

        // If node energy is above migra_max then break now
        if (m_emeasured[i] > logEobs_max) {
            break;
        }

    } // endfor: looped over all nodes

    // If last node energy is below migra_max then compute last part of
    // integral up to emax
    if (x1 < ereco_max) {
        x2    = ereco_max;
        f2    = this->operator()(x2, etrue);
        prob += 0.5 * (f1 + f2) * (x2.log10TeV() - x1.log10TeV());
    }

    // Return probability
    return prob;
}


/***********************************************************************//**
 * @brief Print RMF information
 *
 * @param[in] chatter Chattiness.
 * @return String containing energy dispersion information.
 ***************************************************************************/
std::string GCTAEdispRmf::print(const GChatter& chatter) const
{ 
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAEdispRmf ===");

        // Append energy boundary information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of true energy bins"));
        result.append(gammalib::str(m_rmf.etrue().size()));
        result.append("\n"+gammalib::parformat("Number of measured bins"));
        result.append(gammalib::str(m_rmf.emeasured().size()));
        result.append("\n"+gammalib::parformat("True energy range"));
        result.append(m_rmf.etrue().emin().print());
        result.append(" - ");
        result.append(m_rmf.etrue().emax().print());
        result.append("\n"+gammalib::parformat("Measured energy range"));
        result.append(m_rmf.emeasured().emin().print());
        result.append(" - ");
        result.append(m_rmf.emeasured().emax().print());

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
void GCTAEdispRmf::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_rmf.clear();
    m_matrix.clear();

    // Initialise interpolation cache
    m_etrue.clear();
    m_emeasured.clear();
    m_last_etrue     = 9999.0;
    m_last_emeasured = 9999.0;
    m_itrue1         = 0;
    m_itrue2         = 0;
    m_imeas1         = 0;
    m_imeas2         = 0;
    m_wgt1           = 0.0;
    m_wgt2           = 0.0;
    m_wgt3           = 0.0;
    m_wgt4           = 0.0;

    // Initialise Monte Carlo cache
    m_max_edisp            = 0.0;
    m_last_theta_obs       = -1.0;
    m_last_theta_src       = -1.0;
    m_last_logEsrc         = -30.0;
    m_last_logEobs         = -30.0;
    m_index_obs            = 0;
    m_index_src            = 0;
    m_ebounds_obs_computed = false;
    m_ebounds_src_computed = false;
    m_ebounds_obs.clear();
    m_ebounds_src.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] edisp Energy dispersion
 ***************************************************************************/
void GCTAEdispRmf::copy_members(const GCTAEdispRmf& edisp)
{
    // Copy members
    m_filename = edisp.m_filename;
    m_rmf      = edisp.m_rmf;
    m_matrix   = edisp.m_matrix;

    // Copy interpolation cache
    m_etrue          = edisp.m_etrue;
    m_emeasured      = edisp.m_emeasured;
    m_last_etrue     = edisp.m_last_etrue;
    m_last_emeasured = edisp.m_last_emeasured;
    m_itrue1         = edisp.m_itrue1;
    m_itrue2         = edisp.m_itrue2;
    m_imeas1         = edisp.m_imeas1;
    m_imeas2         = edisp.m_imeas2;
    m_wgt1           = edisp.m_wgt1;
    m_wgt2           = edisp.m_wgt2;
    m_wgt3           = edisp.m_wgt3;
    m_wgt4           = edisp.m_wgt4;

    // Copy Monte Carlo cache
    m_max_edisp            = edisp.m_max_edisp;
    m_last_theta_obs       = edisp.m_last_theta_obs;
    m_last_theta_src       = edisp.m_last_theta_src;
    m_last_logEsrc         = edisp.m_last_logEsrc;
    m_last_logEobs         = edisp.m_last_logEobs;
    m_index_obs            = edisp.m_index_obs;
    m_index_src            = edisp.m_index_src;
    m_ebounds_obs_computed = edisp.m_ebounds_obs_computed;
    m_ebounds_src_computed = edisp.m_ebounds_src_computed;
    m_ebounds_obs          = edisp.m_ebounds_obs;
    m_ebounds_src          = edisp.m_ebounds_src;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEdispRmf::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set redistribution matrix
 ***************************************************************************/
void GCTAEdispRmf::set_matrix(void)
{
    // Determine matrix rows and columns
    int rows    = m_rmf.ntrue();
    int columns = m_rmf.nmeasured();

    // Initialize matrix
    m_matrix = GMatrixSparse(rows, columns);

    // Fill matrix elements
    for (int itrue = 0; itrue < rows; ++itrue) {
        for (int imeasured = 0; imeasured < columns; ++imeasured) {
            m_matrix(itrue, imeasured) = m_rmf(itrue, imeasured);
        }
    }

    // Initialise row sums
    std::vector<double> row_sums;
    row_sums.reserve(rows);

    // Normalize matrix elements
    for (int itrue = 0; itrue < rows; ++itrue) {

        // Get true photon energy
        double logEsrc = m_rmf.etrue().elogmean(itrue).log10TeV();

        // Get integration boundaries
        GEbounds ebounds = ebounds_obs(logEsrc, 0.0);

        // Integrate contributions over energy boundaries
        double sum = 0.0;
        for (int i = 0; i < ebounds.size(); ++i) {
            sum += prob_erecobin(ebounds.emin(i), ebounds.emax(i),
                                 m_rmf.etrue().elogmean(itrue), 0.0);
        }

        // Store matrix sum
        row_sums.push_back(sum);

    } // endfor: looped over all true energies

    // Normalize matrix elements
    for (int itrue = 0; itrue < rows; ++itrue) {
        double sum = row_sums[itrue];
        if (sum > 0.0) {
            for (int imeasured = 0; imeasured < columns; ++imeasured) {
                m_matrix(itrue, imeasured) /= sum;
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set interpolation cache
 *
 * Sets the interpolation cache.
 ***************************************************************************/
void GCTAEdispRmf::set_cache(void) const
{
    // Clear node arrays
    m_etrue.clear();
    m_emeasured.clear();

    // Set log10(Etrue) nodes
    for (int i = 0; i < m_rmf.ntrue(); ++i) {
        m_etrue.append(m_rmf.etrue().elogmean(i).log10TeV());
    }

    // Set log10(Emeasured) nodes
    for (int i = 0; i < m_rmf.nmeasured(); ++i) {
        m_emeasured.append(m_rmf.emeasured().elogmean(i).log10TeV());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set maximum energy dispersion value
 ***************************************************************************/
void GCTAEdispRmf::set_max_edisp(void) const
{
    // Initialise maximum
    m_max_edisp = 0.0;

    // Loop over all response table elements
    for (int i = 0; i < m_matrix.rows(); ++i) {
        for (int k = 0; k < m_matrix.columns(); ++k) {
            double value = m_matrix(i,k);
            if (value > m_max_edisp) {
                m_max_edisp = value;
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update cache
 *
 * @param[in] logEsrc Log10 of true energy (\f$\log_{10}\f$ TeV).
 * @param[in] logEobs Log10 of reconstructed energy (\f$\log_{10}\f$ TeV).
 *
 * Updates the interpolation cache. The interpolation cache is composed
 * of four indices and weights that define 4 data values of the RMF matrix
 * that are used for bilinear interpolation.
 ***************************************************************************/
void GCTAEdispRmf::update(const double& logEsrc, const double& logEobs) const
{
    // Update cache only of arguments have changed
    if (logEsrc != m_last_etrue || logEobs != m_last_emeasured) {

        // Store actual values
        m_last_etrue     = logEsrc;
        m_last_emeasured = logEobs;

        // Set values for node arrays
        m_etrue.set_value(logEsrc);
        m_emeasured.set_value(logEobs);

        // Set indices for bi-linear interpolation
        m_itrue1 = m_etrue.inx_left();
        m_itrue2 = m_etrue.inx_right();
        m_imeas1 = m_emeasured.inx_left();
        m_imeas2 = m_emeasured.inx_right();

        // Set weighting factors for bi-linear interpolation
        m_wgt1 = m_etrue.wgt_left()  * m_emeasured.wgt_left();
        m_wgt2 = m_etrue.wgt_left()  * m_emeasured.wgt_right();
        m_wgt3 = m_etrue.wgt_right() * m_emeasured.wgt_left();
        m_wgt4 = m_etrue.wgt_right() * m_emeasured.wgt_right();

    } // endif: update requested

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute ebounds_obs vector
 *
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 ***************************************************************************/
void GCTAEdispRmf::compute_ebounds_obs(const double& theta,
                                       const double& phi,
                                       const double& zenith,
                                       const double& azimuth) const
{
    // Clear boundaries
    m_ebounds_obs.clear();

    // Loop over Etrue
    for (int i = 0; i < m_rmf.ntrue(); ++i) {

        // Get true photon energy
        double logEsrc = m_rmf.etrue().elogmean(i).log10TeV();

        // Set true energy
        GEnergy etrue;
        etrue.log10TeV(logEsrc);

        // Add ebounds to vector
        m_ebounds_obs.push_back(m_rmf.emeasured(etrue));

    } // endfor: looped over true photon energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute ebounds_src vector
 *
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 ***************************************************************************/
void GCTAEdispRmf::compute_ebounds_src(const double& theta,
                                       const double& phi,
                                       const double& zenith,
                                       const double& azimuth) const
{
    // Clear boundaries
    m_ebounds_src.clear();

    // Loop over Eobs
    for (int i = 0; i < m_rmf.nmeasured(); ++i) {

        // Get observed photon energy
        double logEobs = m_rmf.emeasured().elogmean(i).log10TeV();

        // Set measured energy
        GEnergy emeasured;
        emeasured.log10TeV(logEobs);

        // Add ebounds to vector
        m_ebounds_src.push_back(m_rmf.emeasured(emeasured));

    } // endfor: looped over measured photon energies

    // Return
    return;
}
