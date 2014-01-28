/***************************************************************************
 *  GCTAEdispPerfTable.cpp - CTA performance table energy dispersion class *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAEdispPerfTable.hpp
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
#include "GCTAEdispPerfTable.hpp"
#include "GCTAException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                       "GCTAEdispPerfTable::load(std::string&)"

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
 * Construct instance by loading the point spread function information from
 * an ASCII performance table.
 ***************************************************************************/
GCTAEdispPerfTable::GCTAEdispPerfTable(const std::string& filename) : GCTAEdisp()
{
    // Initialise class members
    init_members();

    // Load point spread function from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] edisp Energy dispersion
 ***************************************************************************/
GCTAEdispPerfTable::GCTAEdispPerfTable(const GCTAEdispPerfTable& edisp) : GCTAEdisp(edisp)
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
 * @brief Return point spread function (in units of sr^-1)
 *
 * @param[in] delta Angular separation between true and measured photon
 *            directions (rad).
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 *
 * Returns the point spread function for a given angular separation in units
 * of sr^-1 for a given energy.
 ***************************************************************************/
double GCTAEdispPerfTable::operator()(const double& logEobs,
                                      const double& logEsrc,
                                      const double& theta,
                                      const double& phi,
                                      const double& zenith,
                                      const double& azimuth) const
{
    // Update the parameter cache
    update(logEsrc);

    // Compute energy dispersion value
    double delta = logEobs - logEsrc;
    double edisp = m_par_scale * std::exp(m_par_width * delta * delta);
    
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
 * @return Deep copy of point spread function instance.
 ***************************************************************************/
GCTAEdispPerfTable* GCTAEdispPerfTable::clone(void) const
{
    return new GCTAEdispPerfTable(*this);
}


/***********************************************************************//**
 * @brief Load point spread function from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the point spread function information from an ASCII
 * performance table.
 *
 * TODO: share file parsing code with GCTAPsfPerfTable instead of duplicating it here.
 ***************************************************************************/
void GCTAEdispPerfTable::load(const std::string& filename)
{
    // Clear arrays
    m_logE.clear();
    m_sigma.clear();

    // Allocate line buffer
    const int n = 1000;
    char  line[n];

    // Expand environment variables
    std::string fname = gammalib::expand_env(filename);

    // Open performance table readonly
    FILE* fptr = std::fopen(fname.c_str(), "r");
    if (fptr == NULL) {
        throw GCTAException::file_open_error(G_LOAD, fname);
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

        // Set conversion factor from 68% containment to 1 sigma
        const double conv = gammalib::gauss68containment2sigma;
        m_sigma.push_back(gammalib::todouble(elements[4])*conv);

    } // endwhile: looped over lines

    // Close file
    std::fclose(fptr);

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Simulate PSF offset (radians)
 *
 * @param[in] ran Random number generator.
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 ***************************************************************************/
GEnergy GCTAEdispPerfTable::mc(GRan&         ran,
                               const double& logE,
                               const double& theta,
                               const double& phi,
                               const double& zenith,
                               const double& azimuth) const
{
    // Update the parameter cache
    update(logE);

    // Draw log observed energy
    double logEobs = ran.normal(0, m_par_sigma);
    double energy = std::pow(10., logEobs);
    
    return GEnergy(energy, "TeV");
}


/***********************************************************************//**
 * @brief Return maximum size of PSF (radians)
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 *
 * Determine the radius beyond which the PSF becomes negligible. This radius
 * is set by this method to \f$5 \times \sigma\f$, where \f$\sigma\f$ is the
 * Gaussian width of the largest PSF component.
 ***************************************************************************/
GEbounds GCTAEdispPerfTable::ebounds(const double& logE,
                                     const double& theta,
                                     const double& phi,
                                     const double& zenith,
                                     const double& azimuth) const
{
    // Update the parameter cache
    update(logE);

    // Compute energy bounds
    const double n_sigmas = 5.0;
    double emin = std::pow(10., logE - n_sigmas * m_par_sigma);
    double emax = std::pow(10., logE + n_sigmas * m_par_sigma);
    GEbounds ebounds(GEnergy(emin, "TeV"), GEnergy(emax, "TeV"));
    
    return ebounds;
}


/***********************************************************************//**
 * @brief Print point spread function information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing point spread function information.
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
        for(int i=0; i < num; ++i) {
          double sigma = m_sigma[i];
          double logE=m_logE[i];
          result.append("\n"+gammalib::str(logE)+"    "+gammalib::str(sigma));
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
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GCTAEdispPerfTable::copy_members(const GCTAEdispPerfTable& psf)
{
    // Copy members
    m_filename  = psf.m_filename;
    m_logE      = psf.m_logE;
    m_sigma     = psf.m_sigma;
    m_par_logE  = psf.m_par_logE;
    m_par_scale = psf.m_par_scale;
    m_par_sigma = psf.m_par_sigma;
    m_par_width = psf.m_par_width;

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
 * @param[in] logE Log10 of the true photon energy (TeV).
 *
 * This method updates the energy dispersion parameter cache.
 * As the performance table energy dispersion only depends on energy,
 * the only parameter on which the cache values depend is the energy.
 ***************************************************************************/
void GCTAEdispPerfTable::update(const double& logE) const
{
    // Only compute energy dispersion parameters if arguments have changed
    if (logE != m_par_logE) {

        // Save energy
        m_par_logE = logE;
    
        // Determine Gaussian sigma
        m_par_sigma = m_logE.interpolate(logE, m_sigma);

        // Derive width=-0.5/(sigma*sigma) and scale=1/(twopi*sigma*sigma)
        double sigma2 = m_par_sigma * m_par_sigma;
        m_par_scale   =  1.0 / (gammalib::twopi * sigma2);
        m_par_width   = -0.5 / sigma2;

    }

    // Return
    return;
}
