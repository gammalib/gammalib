/***************************************************************************
 *  GCTAEdispPerfTable.cpp - CTA performance table energy dispersion class *
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
 * @file GCTAEdispPerfTable.cpp
 * @brief CTA performance table energy dispersion class implementation
 * @author Christoph Deil & Ellis Owen
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>             // std::fopen, std::fgets, and std::fclose
#include <cmath>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GCTAEdispPerfTable.hpp"
#include "GCTAException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                         "GCTAEdispPerfTable::load(GFilename&)"

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
GCTAEdispPerfTable::GCTAEdispPerfTable(void) : GCTAEdisp()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Performance table file name.
 *
 * Construct instance by loading the energy dispersion information from
 * an ASCII performance table.
 ***************************************************************************/
GCTAEdispPerfTable::GCTAEdispPerfTable(const GFilename& filename) :
                    GCTAEdisp()
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
GCTAEdispPerfTable::GCTAEdispPerfTable(const GCTAEdispPerfTable& edisp) :
                    GCTAEdisp(edisp)
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
GCTAEdispPerfTable::~GCTAEdispPerfTable(void)
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
GCTAEdispPerfTable& GCTAEdispPerfTable::operator=(const GCTAEdispPerfTable& edisp)
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
 *    \frac{1}{\sqrt{2\pi}\sigma(E_{\rm true})}
 *    \exp \left(\frac{-(\log_{10} E_{\rm reco} - \log_{10} E_{\rm true})^2}
 *                    {2 \sigma(E_{\rm true})^2} \right) \times
 *    \frac{1}{\log_{10} E_{\rm reco}}
 * \f]
 *
 * in units of MeV\f$^{-1}\f$ where
 * \f$E_{\rm reco}\f$ is the reconstructed energy in units of MeV,
 * \f$E_{\rm true}\f$ is the true energy in units of MeV, and
 * \f$\sigma(E_{\rm true})\f$ is the standard deviation of the energy
 * dispersion that depends on the true photon energy.
 ***************************************************************************/
double GCTAEdispPerfTable::operator()(const GEnergy& ereco,
                                      const GEnergy& etrue,
                                      const double&  theta,
                                      const double&  phi,
                                      const double&  zenith,
                                      const double&  azimuth) const
{
    // Get log10 of true and reconstructed photon energies
    double logEsrc = etrue.log10TeV();
    double logEobs = ereco.log10TeV();

    // Update the parameter cache
    update(logEsrc);

    // Compute energy dispersion value
    double delta = logEobs - logEsrc;
    double edisp = m_par_scale * std::exp(m_par_width * delta * delta);

    // Compute energy dispersion per MeV
    edisp /= (gammalib::ln10 * ereco.MeV());

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
void GCTAEdispPerfTable::clear(void)
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
GCTAEdispPerfTable* GCTAEdispPerfTable::clone(void) const
{
    return new GCTAEdispPerfTable(*this);
}


/***********************************************************************//**
 * @brief Load energy dispersion from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the energy dispersion information from an ASCII
 * performance table. The energy resolution is stored in the 5th column
 * of the performance table as RMS(ln(Eest/Etrue)). The method converts
 * this internally to a sigma value by multiplying the stored values by
 * 1/ln(10).
 ***************************************************************************/
void GCTAEdispPerfTable::load(const GFilename& filename)
{

    // Clear arrays
    m_logE.clear();
    m_sigma.clear();

    // Allocate line buffer
    const int n = 1000;
    char  line[n];

    // Open performance table readonly
    FILE* fptr = std::fopen(filename.url().c_str(), "r");
    if (fptr == NULL) {
        throw GCTAException::file_open_error(G_LOAD, filename.url());
    }

    // Read lines
    while (std::fgets(line, n, fptr) != NULL) {

        // Split line in elements. Strip empty elements from vector.
        std::vector<std::string> elements = gammalib::split(line, " ");
        for (int i = elements.size()-1; i >= 0; i--) {
            if (gammalib::strip_whitespace(elements[i]).length() == 0) {
                elements.erase(elements.begin()+i);
            }
        }

        // Skip header
        if (elements[0].find("log(E)") != std::string::npos) {
            continue;
        }

        // Break loop if end of data table has been reached
        if (elements[0].find("----------") != std::string::npos) {
            break;
        }

        // Push elements in node array and vector
        m_logE.append(gammalib::todouble(elements[0]));

        // The energy resolution is stored as RMS(ln(Eest/Etrue))
        m_sigma.push_back(gammalib::todouble(elements[4]) * gammalib::inv_ln10);

    } // endwhile: looped over lines

    // Close file
    std::fclose(fptr);

    // Store filename
    m_filename = filename;

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
 * Draws observed energy value from a normal distribution of width
 * m_par_sigma around @p logE.
 ***************************************************************************/
GEnergy GCTAEdispPerfTable::mc(GRan&          ran,
                               const GEnergy& etrue,
                               const double&  theta,
                               const double&  phi,
                               const double&  zenith,
                               const double&  azimuth) const
{
    // Get log10 TeV of true photon energy
    double logEsrc = etrue.log10TeV();

    // Update the parameter cache
    update(logEsrc);

    // Draw log observed energy in TeV
    double logEobs = m_par_sigma * ran.normal() + logEsrc;

    // Set energy
    GEnergy energy;
    energy.log10TeV(logEobs);

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Return observed energy interval that contains the energy dispersion.
 *
 * @param[in] etrue True photon energy.
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Boundaries of reconstruced energies.
 *
 * Returns the band of observed energies outside of which the energy
 * dispersion becomes negligible for a given true energy @p logEsrc. This
 * band is set to \f$\pm 5 \times \sigma\f$, where \f$\sigma\f$ is the
 * Gaussian width of the energy dispersion.
 ***************************************************************************/
GEbounds GCTAEdispPerfTable::ereco_bounds(const GEnergy& etrue,
                                          const double&  theta,
                                          const double&  phi,
                                          const double&  zenith,
                                          const double&  azimuth) const
{
    // Set energy band constant
    const double number_of_sigmas = 5.0;

    // Get log10 of true photon energy in TeV
    double etrue_log10TeV = etrue.log10TeV();

    // Get energy dispersion sigma
    double sigma = m_logE.interpolate(etrue_log10TeV, m_sigma);

    // Compute energy boundaries
    GEnergy emin;
    GEnergy emax;
    emin.log10TeV(etrue_log10TeV - number_of_sigmas * sigma);
    emax.log10TeV(etrue_log10TeV + number_of_sigmas * sigma);

    // Return energy boundaries
    return (GEbounds(emin, emax));
}


/***********************************************************************//**
 * @brief Return true energy interval that contains the energy dispersion.
 *
 * @param[in] ereco Reconstructed event energy.
 * @param[in] theta Offset angle in camera system (radians). Not used.
 * @param[in] phi Azimuth angle in camera system (radians). Not used.
 * @param[in] zenith Zenith angle in Earth system (radians). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (radians). Not used.
 * @return Boundaries of true energies.
 *
 * Returns the band of true photon energies outside of which the energy
 * dispersion becomes negligible for a given observed energy @p logEobs.
 * This band is set to \f$\pm 5 \times \sigma\f$, where \f$\sigma\f$ is the
 * Gaussian width of the energy dispersion.
 ***************************************************************************/
GEbounds GCTAEdispPerfTable::etrue_bounds(const GEnergy& ereco,
                                          const double&  theta,
                                          const double&  phi,
                                          const double&  zenith,
                                          const double&  azimuth) const
{
    // Set energy band constant
    const double number_of_sigmas = 5.0;

    // Get log10 of reconstructed event energy in TeV
    double ereco_log10TeV = ereco.log10TeV();

    // Get energy dispersion sigma
    double sigma = m_logE.interpolate(ereco_log10TeV, m_sigma);

    // Compute energy boundaries
    GEnergy emin;
    GEnergy emax;
    emin.log10TeV(ereco_log10TeV - number_of_sigmas * sigma);
    emax.log10TeV(ereco_log10TeV + number_of_sigmas * sigma);

    // Return energy boundaries
    return (GEbounds(emin, emax));
}


/***********************************************************************//**
 * @brief Return energy dispersion probability for reconstructed energy
 *        interval
 *
 * @param[in] ereco_min Minimum of reconstructed energy interval.
 * @param[in] ereco_max Maximum of reconstructed energy interval.
 * @param[in] etrue True energy.
 * @param[in] theta Offset angle (radians). Not used.
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
double GCTAEdispPerfTable::prob_erecobin(const GEnergy& ereco_min,
                                         const GEnergy& ereco_max,
                                         const GEnergy& etrue,
                                         const double&  theta) const
{
    // Update the parameter cache
    update(etrue.log10TeV());

    // Get normalized energy boundaries in log10 energy
    double xmin = (ereco_min.log10TeV() - etrue.log10TeV()) / m_par_sigma;
    double xmax = (ereco_max.log10TeV() - etrue.log10TeV()) / m_par_sigma;

    // Compute fraction of probability within the energy boundaries
    double prob = gammalib::gauss_integral(xmin, xmax);

    // Return
    return prob;
}


/***********************************************************************//**
 * @brief Print energy dispersion information
 *
 * @param[in] chatter Chattiness.
 * @return String containing energy dispersion information.
 ***************************************************************************/
std::string GCTAEdispPerfTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        int    num  = m_logE.size();
        double emin = std::pow(10.0, m_logE[0]);
        double emax = std::pow(10.0, m_logE[num-1]);

        // Append header
        result.append("=== GCTAEdispPerfTable ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(num));
        result.append("\n"+gammalib::parformat("Log10(Energy) range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");

        /*
        for(int i=0; i < num; ++i) {
            double sigma = m_sigma[i];
            double logE=m_logE[i];
            result.append("\n"+gammalib::str(logE)+"    "+gammalib::str(sigma));
        }
        */

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
void GCTAEdispPerfTable::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_logE.clear();
    m_sigma.clear();
    m_par_logE  = -1.0e30;
    m_par_scale = 1.0;
    m_par_sigma = 0.0;
    m_par_width = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] edisp Energy dispersion
 ***************************************************************************/
void GCTAEdispPerfTable::copy_members(const GCTAEdispPerfTable& edisp)
{
    // Copy members
    m_filename  = edisp.m_filename;
    m_logE      = edisp.m_logE;
    m_sigma     = edisp.m_sigma;
    m_par_logE  = edisp.m_par_logE;
    m_par_scale = edisp.m_par_scale;
    m_par_sigma = edisp.m_par_sigma;
    m_par_width = edisp.m_par_width;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEdispPerfTable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update energy dispersion parameter cache
 *
 * @param[in] logE Log10 of the true photon energy (\f$\log_{10}\f$ TeV).
 *
 * This method updates the energy dispersion parameter cache. As the
 * performance table energy dispersion only depends on true photon energy,
 * the only parameter on which the cache values depend is the true photon
 * energy.
 ***************************************************************************/
void GCTAEdispPerfTable::update(const double& logE) const
{
    // Only compute energy dispersion parameters if arguments have changed
    if (logE != m_par_logE) {

        // Save energy
        m_par_logE = logE;

        // Determine Gaussian sigma and pre-compute Gaussian parameters
        m_par_sigma = m_logE.interpolate(logE, m_sigma);
        m_par_scale = gammalib::inv_sqrt2pi / m_par_sigma;
        m_par_width = -0.5 / (m_par_sigma * m_par_sigma);

    }

    // Return
    return;
}
