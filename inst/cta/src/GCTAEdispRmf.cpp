/***************************************************************************
 *            GCTAEdispRmf.cpp - CTA RMF energy dispersion class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2015 by Christoph Deil & Ellis Owen                 *
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
#include "GEnergy.hpp"
#include "GCTAEdispRmf.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                             "GCTAEdispRmf::load(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_OLD_MC_CODE
//#define G_MC_REJECTION

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
GCTAEdispRmf::GCTAEdispRmf(const std::string& filename) : GCTAEdisp()
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
 * @brief Return energy dispersion.
 *
 * @param[in] logEobs log10 of the observed photon energy (TeV).
 * @param[in] logEsrc log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 *
 * Returns the energy resolution, i.e. the probability density in observed
 * photon energy at a given (log10(E_src), log10(E_obs)).
 * To be precise: energy dispersion = dP / d(log10(E_obs)).
 *
 ***************************************************************************/
double GCTAEdispRmf::operator()(const double& logEobs,
                                const double& logEsrc,
                                const double& theta,
                                const double& phi,
                                const double& zenith,
                                const double& azimuth) const
{
    // Update indexes and weighting factors for interpolation
    update(logEobs, logEsrc);

    // Perform interpolation
    double edisp =  m_wgt1 * m_matrix(m_itrue1, m_imeas1) +
                    m_wgt2 * m_matrix(m_itrue1, m_imeas2) +
                    m_wgt3 * m_matrix(m_itrue2, m_imeas1) +
                    m_wgt4 * m_matrix(m_itrue2, m_imeas2);

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
void GCTAEdispRmf::load(const std::string& filename)
{
    // Load RMF file
    m_rmf.load(filename);

    // Store the filename
    m_filename = filename;

    // Set matrix
    set_matrix();

    // Set cache
    set_cache();
    set_mc_cache();

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
 * Draws observed energy value from RMF matrix.
 *
 * @todo The new method is very computation intensitive as the cumulative
 *       distributions are computed for every event. The old method is
 *       more efficient. Either implement a rejection method for smoothness
 *       or introduce subsampling for CDF cache.
 ***************************************************************************/
#if defined(G_OLD_MC_CODE)
GEnergy GCTAEdispRmf::mc(GRan&         ran,
                         const double& logEsrc,
                         const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth) const
{
    // Set true energy
    GEnergy energy;
    energy.log10TeV(logEsrc);

    // Determine true energy index
    int itrue = m_rmf.etrue().index(energy);

    // Continue only if the true energy index lies within a true energy
    // boundary of the redistribution matrix file
    if (itrue != -1) {

        // Get offset in measured energy. Continue only if offset is valid
        int offset = m_mc_measured_start[itrue];
        if (offset != -1) {

            #if defined(G_MC_REJECTION)
            // Determine measured energy index from Monte-Carlo cache
            int imeasured = ran.cdf(m_mc_measured_cdf[itrue]) + offset;

            // Get boundaries for observed energy and the function values
            // at the edges and the centre of the selected bin
            double emin  = m_rmf.emeasured().emin(imeasured).log10TeV();
            double emax  = m_rmf.emeasured().emax(imeasured).log10TeV();
            double emean = m_rmf.emeasured().emean(imeasured).log10TeV();
            double fmin  = operator()(emin,  logEsrc, theta, phi, zenith, azimuth);
            double fmax  = operator()(emax,  logEsrc, theta, phi, zenith, azimuth);
            double fmean = operator()(emean, logEsrc, theta, phi, zenith, azimuth);

            // Get function maximum
            if (fmean > fmax) {
                fmax = fmean;
            }
            if (fmin > fmax) {
                fmax = fmin;
            }

            // Add energy shift to account for the fact that the true energy
            // for which the redistribution matrix is given is not the same
            // as the specified true energy logEsrc.
            double eshift = logEsrc - m_rmf.etrue().emean(itrue).log10TeV();
            emin         += eshift;
            emax         += eshift;

            // Find energy by rejection method
            double ewidth = emax - emin;
            double e      = emin;
            double f      = 0.0;
            double ftest  = 1.0;
            while (ftest > f) {
                e      = emin + ewidth * ran.uniform();
                f      = operator()(e, logEsrc, theta, phi, zenith, azimuth);
                ftest  = ran.uniform() * fmax;
            }
            energy.log10TeV(e);

            #else
            // Determine measured energy index from Monte-Carlo cache
            int imeasured = ran.cdf(m_mc_measured_cdf[itrue]) + offset;

            // Get log10 energy minimum and bin width
            double emin   = m_rmf.emeasured().emin(imeasured).log10TeV();
            double ewidth = m_rmf.emeasured().emax(imeasured).log10TeV() - emin;

            // Draw a random energy from this interval
            double e = emin + ewidth * ran.uniform();

            // Set interval
            energy.log10TeV(e);
            #endif

        } // endif: offset was valid

    } // endif: there was redistribution information

    // Return energy
    return energy;
}
#else
GEnergy GCTAEdispRmf::mc(GRan&         ran,
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
    int high = m_rmf.ntrue();
    while ((high-low) > 1) {
        int mid = (low+high) / 2;
        if (logEsrc < m_rmf.etrue().elogmean(mid).log10TeV()) {
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
    int index = low;

    // Interpolate Eobs value
    double Eobs =   (m_cumul[isrc][index+1].second - p)*m_cumul[isrc][index].first
                  + (p - m_cumul[isrc][index].second)*m_cumul[isrc][index+1].first;
    Eobs       /=   (m_cumul[isrc][index+1].second - m_cumul[isrc][index].second);

    // Set energy
    GEnergy energy;
    energy.TeV(Eobs);


    // Return energy
    return energy;
}
#endif


/***********************************************************************//**
 * @brief Return observed energy interval that contains the energy dispersion.
 *
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
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
    // Set true energy
    GEnergy etrue;
    etrue.log10TeV(logEsrc);

    // Return measured energy boundaries
    return (m_rmf.emeasured(etrue));
}


/***********************************************************************//**
 * @brief Return true energy interval that contains the energy dispersion.
 *
 * @param[in] logEobs Log10 of the observed event energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
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
    // Set measured energy
    GEnergy emeasured;
    emeasured.log10TeV(logEobs);

    // Return true energy boundaries
    return (m_rmf.etrue(emeasured));
}


/***********************************************************************//**
 * @brief Print RMF information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
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
    m_mc_measured_start.clear();
    m_mc_measured_cdf.clear();
    m_cdf_computed   = false;
    m_theta          = 0.0;
    m_cumul.clear();

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
    m_mc_measured_start = edisp.m_mc_measured_start;
    m_mc_measured_cdf   = edisp.m_mc_measured_cdf;
    m_theta             = edisp.m_theta;
    m_cumul             = edisp.m_cumul;
    m_cdf_computed      = edisp.m_cdf_computed;

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

    // Set matrix elements
    for (int itrue = 0; itrue < rows; ++itrue) {

        // Initialise sum
        double sum = 0.0;

        // Sum all measured energy bins
        for (int imeasured = 0; imeasured < columns; ++imeasured) {
            double ewidth = m_rmf.emeasured().ewidth(imeasured).TeV();
            double emean  = m_rmf.emeasured().emean(imeasured).TeV();
            sum          += m_rmf(itrue, imeasured) *
                            ewidth / (gammalib::ln10 * emean);
        }

        // If sum is positive then scale all measured energy bins
        if (sum > 0.0) {
            for (int imeasured = 0; imeasured < columns; ++imeasured) {
                m_matrix(itrue, imeasured) = m_rmf(itrue, imeasured) / sum;
            }
        }

    } // endfor: looped over all true energies

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
        m_etrue.append(m_rmf.etrue().emean(i).log10TeV());
    }

    // Set log10(Emeasured) nodes
    for (int i = 0; i < m_rmf.nmeasured(); ++i) {
        m_emeasured.append(m_rmf.emeasured().emean(i).log10TeV());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update the Monte-Carlo cache values
 *
 * Sets the following private class members:
 *
 *      m_mc_measured_start: start index in emeasured boundaries for each
 *                           true energy
 *      m_mc_measured_cdf:   CDF for each true energy
 *
 ***************************************************************************/
void GCTAEdispRmf::set_mc_cache(void) const
{
    // Clear MC cache
    m_mc_measured_start.clear();
    m_mc_measured_cdf.clear();

    // Reserve some space
    m_mc_measured_start.reserve(m_rmf.ntrue());
    m_mc_measured_cdf.reserve(m_rmf.ntrue());

    // Loop over true energies
    for (int itrue = 0; itrue < m_rmf.ntrue(); ++itrue) {

        // Determine first non-zero measured energy
        int imeasured_start = -1;
        for (int i = 0; i < m_rmf.nmeasured(); ++i) {
            if (m_rmf(itrue, i) > 0.0) {
                imeasured_start = i;
                break;
            }
        }
        m_mc_measured_start.push_back(imeasured_start);

        // Determine last non-zero measured energy
        int imeasured_stop = -1;
        for (int i = m_rmf.nmeasured()-1; i >= 0; --i) {
            if (m_rmf(itrue, i) > 0.0) {
                imeasured_stop = i;
                break;
            }
        }

        // Determine number of elements
        int num = (imeasured_stop - imeasured_start + 1);
        if (num < 0) {
            num = 0;
        }

        // Allocate vector
        GVector vector(num+1);

        // If there are measured energies for this true energy then build
        // now the CDF ...
        if (imeasured_start != -1 && imeasured_stop != -1) {
            double sum = 0.0;
            vector[0] = 0.0;
            for (int i = imeasured_start, k = 1; i <= imeasured_stop; ++i, ++k) {
                sum      += m_rmf(itrue, i);
                vector[k] = sum;
            }
            if (sum > 0.0) {
                for (int k = 0; k < num; ++k) {
                    vector[k] /= sum;
                }
            }
        }

        // Make sure that last pixel in the cache is >1
        vector[num] = 1.0001;

        // Push vector on cache
        m_mc_measured_cdf.push_back(vector);

    } // endfor: looped over all true energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update cache
 *
 * @param[in] etrue True energy.
 * @param[in] emeasured Measured energy.
 *
 * Updates the interpolation cache. The interpolation cache is composed
 * of four indices and weights that define 4 data values of the RMF matrix
 * that are used for bilinear interpolation.
 ***************************************************************************/
void GCTAEdispRmf::update(const double& etrue, const double& emeasured) const
{
    // Update cache only of arguments have changed
    if (etrue != m_last_etrue || emeasured != m_last_emeasured) {

        // Store actual values
        m_last_etrue     = etrue;
        m_last_emeasured = emeasured;

        // Set values for node arrays
        m_etrue.set_value(etrue);
        m_emeasured.set_value(emeasured);

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
 * @brief Compute cumulative probability
 *
 ***************************************************************************/
void GCTAEdispRmf::compute_cumul(const double& theta,
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
    for (int isrc = 0; isrc < m_rmf.ntrue(); ++isrc) {

        m_cumul.push_back(std::vector<std::pair<double, double> >());

        double logEsrc = m_rmf.etrue().elogmean(isrc).log10TeV();

        // Initialize cumulative probability
        double sum = 0.0;

        // Divide bin into n sub-bins
        const int n = 10;

        // Loop through Eobs
        for (int imeasured = 0; imeasured < m_rmf.nmeasured(); ++imeasured) {

            // Compute deltaEobs and logEobs values
            double deltaEobs   = m_rmf.emeasured().ewidth(imeasured).TeV();
            double Eobsmin     = m_rmf.emeasured().emin(imeasured).TeV();

            for (int i = 0; i < n; ++i) {

                double Eobs = Eobsmin+i*deltaEobs/n;

                // Compute cumulative probability
                double add = GCTAEdispRmf::operator()(logEsrc, std::log10(Eobs), theta) * deltaEobs / (n * Eobs * gammalib::ln10);

                // Add to sum (and keep sum lower or equal to 1)
                sum = sum+add >= 1.0 ? 1.0 : sum+add;

                // Create pair containing Eobs and cumulative probability
                std::pair<double, double> pair(Eobs, sum);

                // Add to vector
                m_cumul[isrc].push_back(pair);
            }

        }
    }

    }
    // Return
    return;
}
