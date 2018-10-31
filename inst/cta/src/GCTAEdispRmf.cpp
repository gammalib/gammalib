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
#define G_MC  "GCTAEdispRmf::mc(GRan&, GEnergy&, double&, double&, double&, "\
                                                                   "double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_LINEAR_ERECO_INTERPOLATION    //!< Linear interpolation for Ereco

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
 *    E_{\rm disp}(E_{\rm reco} | E_{\rm true})
 * \f]
 *
 * in units of MeV\f$^{-1}\f$ where
 * \f$E_{\rm reco}\f$ is the reconstructed energy, and
 * \f$E_{\rm true}\f$ is the true energy.
 ***************************************************************************/
double GCTAEdispRmf::operator()(const GEnergy& ereco,
                                const GEnergy& etrue,
                                const double&  theta,
                                const double&  phi,
                                const double&  zenith,
                                const double&  azimuth) const
{
    // Update indexes and weighting factors for interpolation
    update(ereco, etrue);

    // Perform interpolation
    double edisp =  m_wgt1 * m_matrix(m_itrue1, m_ireco1) +
                    m_wgt2 * m_matrix(m_itrue1, m_ireco2) +
                    m_wgt3 * m_matrix(m_itrue2, m_ireco1) +
                    m_wgt4 * m_matrix(m_itrue2, m_ireco2);

    // Make sure that energy dispersion is not negative
    if (edisp < 0.0) {
        edisp = 0.0;
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
 * @return Reconstructed energy.
 *
 * @exception GException::invalid_return_value
 *            Energy dispersion matrix is empty.
 *
 * Draws reconstructed energy value from RMF matrix given a true energy
 * @p etrue. If no energy dispersion information is available the method
 * will return the true photon energy.
 ***************************************************************************/
GEnergy GCTAEdispRmf::mc(GRan&          ran,
                         const GEnergy& etrue,
                         const double&  theta,
                         const double&  phi,
                         const double&  zenith,
                         const double&  azimuth) const
{
    // Initialise reconstructed event energy with true photon energy
    GEnergy ereco = etrue;

    // Get boundaries for reconstructed energy
    GEbounds ebounds = ereco_bounds(etrue, theta, phi, zenith, azimuth);
    double   emin    = ebounds.emin().TeV();
    double   emax    = ebounds.emax().TeV();

    // Get maximum energy dispersion value (including some margin)
    double max_edisp = 1.5 * m_max_edisp;

    // Throw an exception if maximum energy dispersion is zero
    if (max_edisp <= 0.0) {
        std::string msg = "Energy dispersion matrix is empty. Please provide "
                          "a valid energy dispersion matrix.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Find energy by rejection method
    double ewidth  = emax - emin;
    if (max_edisp > 0.0) {
        double f      = 0.0;
        double ftest  = 1.0;
        while (ftest > f) {
            ereco.TeV(emin + ewidth * ran.uniform());
            f     = operator()(ereco, etrue, theta, phi, zenith, azimuth);
            ftest = ran.uniform() * max_edisp;
        }
    }

    // Return reconstructed photon energy
    return ereco;
}


/***********************************************************************//**
 * @brief Return observed energy interval that contains the energy dispersion.
 *
 * @param[in] etrue True photon energy.
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Reconstructed energy boundaries.
 *
 * Returns the band of observed energies outside of which the energy
 * dispersion becomes negligible for a given true energy @p etrue.
 ***************************************************************************/
GEbounds GCTAEdispRmf::ereco_bounds(const GEnergy& etrue,
                                    const double&  theta,
                                    const double&  phi,
                                    const double&  zenith,
                                    const double&  azimuth) const
{
    // Compute only if parameters changed
    if (!m_ereco_bounds_computed) {

        // Store last theta value
        m_ereco_bounds_computed = true;
        m_last_etrue_bounds.TeV(0.0); // force update

        // Compute m_ereco_bounds
        compute_ereco_bounds();

    }

    // Search index only if logEsrc has changed
    if (etrue != m_last_etrue_bounds) {

        // Store last photon energy
        m_last_etrue_bounds = etrue;

        // Get true photon energy in TeV
        double etrue_TeV = etrue.TeV();

        // Find right index with bisection
        int low  = 0;
        int high = m_ereco_bounds.size() - 1;
        while ((high-low) > 1) {
            int  mid = (low+high) / 2;
            double e = m_rmf.etrue().emin(mid).TeV();
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
    return (m_ereco_bounds[m_index_ereco]);
}


/***********************************************************************//**
 * @brief Return true energy interval that contains the energy dispersion.
 *
 * @param[in] ereco Reconstructed event energy.
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return True energy boundaries.
 *
 * Returns the band of true photon energies outside of which the energy
 * dispersion becomes negligible for a given observed energy @p logEobs.
 ***************************************************************************/
GEbounds GCTAEdispRmf::etrue_bounds(const GEnergy& ereco,
                                    const double&  theta,
                                    const double&  phi,
                                    const double&  zenith,
                                    const double&  azimuth) const
{
    // Compute only if parameters changed
    if (!m_etrue_bounds_computed) {

        // Set computation flag
        m_etrue_bounds_computed = true;
        m_last_ereco_bounds.TeV(0.0); // force update

        // Compute m_etrue_bounds
        compute_etrue_bounds();

    }

    // Search index only if ereco has changed
    if (ereco != m_last_ereco_bounds) {

        // Store last value
        m_last_ereco_bounds = ereco;

        // Get reconstructed energy in TeV
        double ereco_TeV = ereco.TeV();

        // Find right index with bisection
        int low  = 0;
        int high = m_etrue_bounds.size() - 1;
        while ((high-low) > 1) {
            int  mid = (low+high) / 2;
            double e = m_rmf.emeasured().emin(mid).TeV();
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
    return (m_etrue_bounds[m_index_etrue]);
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
 *    E_{\rm disp}(E_{\rm reco} | E_{\rm true}) \, dE_{\rm reco}
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
    #if defined(G_LINEAR_ERECO_INTERPOLATION)
    double logEobs_min = ereco_min.TeV();
    double logEobs_max = ereco_max.TeV();
    #else
    double logEobs_min = ereco_min.log10TeV();
    double logEobs_max = ereco_max.log10TeV();
    #endif

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
    for (int i = 0; i < m_ereco.size(); ++i) {

        // If measured energy is below logEobs_min then skip node
        if (m_ereco[i] <= logEobs_min) {
            continue;
        }

        // If measured energy is above maximum migration value then use
        // logEobs_max as measured energy
        if (m_ereco[i] > logEobs_max) {
            x2 = ereco_max;
            f2 = this->operator()(x2, etrue);
        }

        // ... otherwise use measured energy
        else {
            #if defined(G_LINEAR_ERECO_INTERPOLATION)
            x2.TeV(m_ereco[i]);
            #else
            x2.log10TeV(m_ereco[i]);
            #endif
            f2 = wgt_left  * m_matrix(inx_left,  i) +
                 wgt_right * m_matrix(inx_right, i);
            if (f2 < 0.0) {
                f2 = 0.0;
            }
        }

        // Compute integral
        prob += 0.5 * (f1 + f2) * (x2.MeV() - x1.MeV());

        // Set second node as first node
        x1 = x2;
        f1 = f2;

        // If node energy is above migra_max then break now
        if (m_ereco[i] > logEobs_max) {
            break;
        }

    } // endfor: looped over all nodes

    // If last node energy is below migra_max then compute last part of
    // integral up to emax
    if (x1 < ereco_max) {
        x2    = ereco_max;
        f2    = this->operator()(x2, etrue);
        prob += 0.5 * (f1 + f2) * (x2.MeV() - x1.MeV());
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
    m_max_edisp = 0.0;

    // Initialise interpolation cache
    m_etrue.clear();
    m_ereco.clear();
    m_last_etrue.clear();
    m_last_ereco.clear();
    m_itrue1     = 0;
    m_itrue2     = 0;
    m_ireco1     = 0;
    m_ireco2     = 0;
    m_wgt1       = 0.0;
    m_wgt2       = 0.0;
    m_wgt3       = 0.0;
    m_wgt4       = 0.0;

    // Initialise computation cache
    m_ereco_bounds_computed = false;
    m_etrue_bounds_computed = false;
    m_index_ereco           = 0;
    m_index_etrue           = 0;
    m_last_etrue_bounds.clear();
    m_last_ereco_bounds.clear();
    m_ereco_bounds.clear();
    m_etrue_bounds.clear();

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
    m_filename  = edisp.m_filename;
    m_rmf       = edisp.m_rmf;
    m_matrix    = edisp.m_matrix;
    m_max_edisp = edisp.m_max_edisp;

    // Copy interpolation cache
    m_etrue      = edisp.m_etrue;
    m_ereco      = edisp.m_ereco;
    m_last_etrue = edisp.m_last_etrue;
    m_last_ereco = edisp.m_last_ereco;
    m_itrue1     = edisp.m_itrue1;
    m_itrue2     = edisp.m_itrue2;
    m_ireco1     = edisp.m_ireco1;
    m_ireco2     = edisp.m_ireco2;
    m_wgt1       = edisp.m_wgt1;
    m_wgt2       = edisp.m_wgt2;
    m_wgt3       = edisp.m_wgt3;
    m_wgt4       = edisp.m_wgt4;

    // Copy computation cache
    m_ereco_bounds_computed = edisp.m_ereco_bounds_computed;
    m_etrue_bounds_computed = edisp.m_etrue_bounds_computed;
    m_index_ereco           = edisp.m_index_ereco;
    m_index_etrue           = edisp.m_index_etrue;
    m_last_etrue_bounds     = edisp.m_last_etrue_bounds;
    m_last_ereco_bounds     = edisp.m_last_ereco_bounds;
    m_ereco_bounds          = edisp.m_ereco_bounds;
    m_etrue_bounds          = edisp.m_etrue_bounds;

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

    // Get reconstructued energy bin boundaries
    GEbounds ebounds = m_rmf.emeasured();

    // Fill matrix elements
    for (int ireco = 0; ireco < columns; ++ireco) {

        // Compute reconstructued energy bin width in MeV
        double ewidth = ebounds.emax(ireco).MeV() - ebounds.emin(ireco).MeV();

        // If bin width is positive then fill all matrix bins by deviding
        // the RMF values by the bin width. This assures that the stored
        // values are per MeV of reconstructed energy.
        if (ewidth > 0.0) {
            for (int itrue = 0; itrue < rows; ++itrue) {
                m_matrix(itrue, ireco) = m_rmf(itrue, ireco) / ewidth;
            }
        }

    } // endfor: looped over all reconstructed energies

    // Initialise row sums
    std::vector<double> row_sums(rows, 0.0);

    // Normalize matrix elements
    for (int itrue = 0; itrue < rows; ++itrue) {

        // Get true energy
        GEnergy etrue = m_rmf.etrue().elogmean(itrue);

        // Get integration boundaries
        GEbounds ebounds = ereco_bounds(etrue);

        // Integrate contributions over energy boundaries
        double sum = 0.0;
        for (int i = 0; i < ebounds.size(); ++i) {
            sum += prob_erecobin(ebounds.emin(i), ebounds.emax(i), etrue, 0.0);
        }

        // Store matrix sum
        row_sums.push_back(sum);

    } // endfor: looped over all true energies

    // Normalize matrix elements
    for (int itrue = 0; itrue < rows; ++itrue) {
        double sum = row_sums[itrue];
        if (sum > 0.0) {
            for (int ireco = 0; ireco < columns; ++ireco) {
                m_matrix(itrue, ireco) /= sum;
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
    m_ereco.clear();

    // Set log10(Etrue) nodes
    for (int i = 0; i < m_rmf.ntrue(); ++i) {
        m_etrue.append(m_rmf.etrue().elogmean(i).log10TeV());
    }

    // Set log10(Emeasured) nodes
    for (int i = 0; i < m_rmf.nmeasured(); ++i) {
        #if defined(G_LINEAR_ERECO_INTERPOLATION)
        m_ereco.append(m_rmf.emeasured().elogmean(i).TeV());
        #else
        m_ereco.append(m_rmf.emeasured().elogmean(i).log10TeV());
        #endif
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set maximum energy dispersion value
 ***************************************************************************/
void GCTAEdispRmf::set_max_edisp(void)
{
    // Initialise maximum
    m_max_edisp = 0.0;

    // Loop over all response table elements
    for (int itrue = 0; itrue < m_matrix.rows(); ++itrue) {
        for (int ireco = 0; ireco < m_matrix.columns(); ++ireco) {
            double value = m_matrix(itrue, ireco);
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
 * @param[in] ereco Reconstructed energy.
 * @param[in] etrue True energy.
 *
 * Updates the interpolation cache. The interpolation cache is composed
 * of four indices and weights that define 4 data values of the RMF matrix
 * that are used for bilinear interpolation.
 ***************************************************************************/
void GCTAEdispRmf::update(const GEnergy& ereco, const GEnergy& etrue) const
{
    // Update cache only of arguments have changed
    if (ereco != m_last_etrue || etrue != m_last_ereco) {

        // Store actual values
        m_last_ereco = ereco;
        m_last_etrue = etrue;

        // Get log10 of true and reconstructed photon energies
        double logEsrc = etrue.log10TeV();
        #if defined(G_LINEAR_ERECO_INTERPOLATION)
        double logEobs = ereco.TeV();
        #else
        double logEobs = ereco.log10TeV();
        #endif

        // Set values for node arrays
        m_etrue.set_value(logEsrc);
        m_ereco.set_value(logEobs);

        // Set indices for bi-linear interpolation
        m_itrue1 = m_etrue.inx_left();
        m_itrue2 = m_etrue.inx_right();
        m_ireco1 = m_ereco.inx_left();
        m_ireco2 = m_ereco.inx_right();

        // Set weighting factors for bi-linear interpolation
        m_wgt1 = m_etrue.wgt_left()  * m_ereco.wgt_left();
        m_wgt2 = m_etrue.wgt_left()  * m_ereco.wgt_right();
        m_wgt3 = m_etrue.wgt_right() * m_ereco.wgt_left();
        m_wgt4 = m_etrue.wgt_right() * m_ereco.wgt_right();

    } // endif: update requested

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute m_ereco_bounds vector
 ***************************************************************************/
void GCTAEdispRmf::compute_ereco_bounds(void) const
{
    // Clear boundaries
    m_ereco_bounds.clear();

    // Loop over Etrue
    for (int i = 0; i < m_rmf.ntrue(); ++i) {

        // Get true photon energy
        GEnergy etrue = m_rmf.etrue().elogmean(i);

        // Add ebounds to vector
        m_ereco_bounds.push_back(m_rmf.emeasured(etrue));

    } // endfor: looped over true photon energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute m_etrue_bounds vector
 ***************************************************************************/
void GCTAEdispRmf::compute_etrue_bounds(void) const
{
    // Clear boundaries
    m_etrue_bounds.clear();

    // Loop over Eobs
    for (int i = 0; i < m_rmf.nmeasured(); ++i) {

        // Get reconstructed energy
        GEnergy ereco = m_rmf.emeasured().elogmean(i);

        // Add ebounds to vector
        m_etrue_bounds.push_back(m_rmf.emeasured(ereco));

    } // endfor: looped over measured photon energies

    // Return
    return;
}
