/***************************************************************************
 *   GCTABackgroundPerfTable.cpp - CTA performance table background class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTABackgroundPerfTable.cpp
 * @brief CTA performance table background class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>             // std::fopen, std::fgets, and std::fclose
#include "GTools.hpp"
#include "GMath.hpp"
#include "GIntegral.hpp"
#include "GCTABackgroundPerfTable.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                  "GCTABackgroundPerfTable::load(std::string&)"
#define G_MC           "GCTABackgroundPerfTable::mc(GEnergy&, GTime&, GRan&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_MC            //!< Debug Monte Carlo method
#define G_LOG_INTERPOLATION   //!< Energy interpolate log10(background rate)

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTABackgroundPerfTable::GCTABackgroundPerfTable(void) : GCTABackground()
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
 * Construct instance by loading the background information from a
 * performance table.
 ***************************************************************************/
GCTABackgroundPerfTable::GCTABackgroundPerfTable(const std::string& filename) :
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
 ***************************************************************************/
GCTABackgroundPerfTable::GCTABackgroundPerfTable(const GCTABackgroundPerfTable& bgd) :
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
 ***************************************************************************/
GCTABackgroundPerfTable::~GCTABackgroundPerfTable(void)
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
 ***************************************************************************/
GCTABackgroundPerfTable& GCTABackgroundPerfTable::operator=(const GCTABackgroundPerfTable& bgd)
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
 * @brief Return background rate in units of events/s/MeV/sr
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] detx Tangential coord in nominal sys (rad).
 * @param[in] dety Tangential coord in nominal sys (rad).
 * @param[in] etrue Use true energy (true/false) (not used).
 *
 * Returns the background rate in units of events/s/MeV/sr for a given energy
 * and detector coordinates. The method assures that the background rate
 * never becomes negative.
 *
 * The method supports true and reconstructed energies for logE. To access
 * the background rate as function of true energy, specify etrue=true
 * (this is the default). The obtained the background rate as function of
 * reconstructed energy, specify etrue=false.
 ***************************************************************************/
double GCTABackgroundPerfTable::operator()(const double& logE, 
                                           const double& detx, 
                                           const double& dety,
                                           const bool&   etrue) const
{
    // Get background rate
    double rate = m_logE.interpolate(logE, m_background);
    #if defined(G_LOG_INTERPOLATION)
    rate = std::pow(10.0, rate);
    #endif

    // Make sure that background rate is not negative
    if (rate < 0.0) {
        rate = 0.0;
    }

    // Optionally add in Gaussian offset angle dependence
    if (m_sigma != 0.0) {
        double theta  = std::sqrt(detx*detx + dety*dety) * gammalib::rad2deg;
        double arg    = theta * theta / m_sigma;
        double scale  = std::exp(-0.5 * arg * arg);
        rate         *= scale;
    }
    
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
 * This method properly resets the background to an initial state.
 ***************************************************************************/
void GCTABackgroundPerfTable::clear(void)
{
    // Free class members (base and derived classes, derived class first)
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
 * @return Pointer to deep copy of background.
 ***************************************************************************/
GCTABackgroundPerfTable* GCTABackgroundPerfTable::clone(void) const
{
    return new GCTABackgroundPerfTable(*this);
}


/***********************************************************************//**
 * @brief Load background from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GException::file_open_error
 *            File could not be opened for read access.
 *
 * Loads the background information from a performance table.
 ***************************************************************************/
void GCTABackgroundPerfTable::load(const std::string& filename)
{
    // Clear arrays
    m_logE.clear();
    m_background.clear();

    // Allocate line buffer
    const int n = 1000;
    char  line[n];

    // Expand environment variables
    std::string fname = gammalib::expand_env(filename);

    // Open performance table readonly
    FILE* fptr = std::fopen(fname.c_str(), "r");
    if (fptr == NULL) {
        throw GException::file_open_error(G_LOAD, fname);
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

        // Determine on-axis background rate (counts/s/MeV/sr)
        double logE       = gammalib::todouble(elements[0]);
        double r80        = gammalib::todouble(elements[3]) * gammalib::deg2rad;
        double bgrate     = gammalib::todouble(elements[5]);         // in Hz
		double emin       = std::pow(10.0, logE-0.1) * 1.0e6;
		double emax       = std::pow(10.0, logE+0.1) * 1.0e6;
		double ewidth     = emax - emin;                             // in MeV
		double solidangle = gammalib::twopi * (1.0 - std::cos(r80)); // in sr
        if (solidangle > 0.0) {
            bgrate /= (solidangle * ewidth); // counts/s/MeV/sr
        }
        else {
            bgrate = 0.0;
        }

        // Push elements in node array and vector
        m_logE.append(logE);
        #if defined(G_LOG_INTERPOLATION)
        m_background.push_back(std::log10(bgrate));
        #else
        m_background.push_back(bgrate);
        #endif

    } // endwhile: looped over lines

    // Close file
    std::fclose(fptr);

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return CTA instrument direction.
 ***************************************************************************/
GCTAInstDir GCTABackgroundPerfTable::mc(const GEnergy& energy,
                                        const GTime&   time,
                                        GRan&          ran) const
{
    // Simulate theta angle
    #if defined(G_DEBUG_MC)
    int    n_samples = 0;
    #endif
    double sigma_max = 4.0 * std::sqrt(sigma());
    double u_max     = std::sin(sigma_max * gammalib::deg2rad);
    double value     = 0.0;
    double u         = 1.0;
    double theta     = 0.0;
    do {
        theta       = ran.uniform() * sigma_max;
        double arg  = theta * theta / sigma();
        double arg2 = arg * arg;
        value       = std::sin(theta * gammalib::deg2rad) * exp(-0.5 * arg2);
        u           = ran.uniform() * u_max;
        #if defined(G_DEBUG_MC)
        n_samples++;
        #endif
    } while (u > value);
    theta *= gammalib::deg2rad;
    #if defined(G_DEBUG_MC)
    std::cout << "#=" << n_samples << " ";
    #endif

    // Simulate azimuth angle
    double phi = gammalib::twopi * ran.uniform();

	// Compute detx and dety
    double detx(0.0);
    double dety(0.0);
	if (theta > 0.0 ) {
		detx = theta * std::cos(phi);
		dety = theta * std::sin(phi);
	}

    // Set instrument direction (in radians)
    GCTAInstDir dir;
    dir.detx(detx);
    dir.dety(dety);

    // Return instrument direction
    return dir;
}


/***********************************************************************//**
 * @brief Print background information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing background information.
 ***************************************************************************/
std::string GCTABackgroundPerfTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        double emin = std::pow(10.0, m_logE[0]);
        double emax = std::pow(10.0, m_logE[size()-1]);

        // Append header
        result.append("=== GCTABackgroundPerfTable ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Log10(Energy) range"));
        result.append(gammalib::str(emin) +
                      " - " +
                      gammalib::str(emax) +
                      " TeV");

        // Append offset angle dependence
        if (m_sigma == 0) {
            result.append("\n"+gammalib::parformat("Offset angle dependence") +
                          "none");
        }
        else {
            std::string txt = "Fixed sigma=" + gammalib::str(m_sigma);
            result.append("\n"+gammalib::parformat("Offset angle dependence") +
                          txt);
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
void GCTABackgroundPerfTable::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_logE.clear();
    m_background.clear();
    m_sigma = 3.0;

    // Initialise MC cache
    m_mc_spectrum.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bgd Background.
 ***************************************************************************/
void GCTABackgroundPerfTable::copy_members(const GCTABackgroundPerfTable& bgd)
{
    // Copy members
    m_filename   = bgd.m_filename;
    m_logE       = bgd.m_logE;
    m_background = bgd.m_background;
    m_sigma      = bgd.m_sigma;

    // Copy MC cache
    m_mc_spectrum = bgd.m_mc_spectrum;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTABackgroundPerfTable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns integral over radial model (in steradians)
 *
 * Computes
 * \f[\Omega = 2 \pi \int_0^{\pi} \sin \theta f(\theta) d\theta\f]
 * where
 * \f[f(\theta) = \exp \left(-\frac{1}{2}
 *                     \left( \frac{\theta^2}{\sigma} \right)^2 \right)\f]
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$\sigma\f$ is the Gaussian width (in degrees\f$^2\f$).
 *
 * The integration is performed numerically, and the upper integration bound
 * \f$\pi\f$
 * is set to
 * \f$\sqrt(10 \sigma)\f$
 * to reduce inaccuracies in the numerical integration.
 ***************************************************************************/
double GCTABackgroundPerfTable::solidangle(void) const
{
    // Allocate integrand
    GCTABackgroundPerfTable::integrand integrand(sigma() *
                                                 gammalib::deg2rad * 
                                                 gammalib::deg2rad);

    // Allocate intergal
    GIntegral integral(&integrand);

    // Set upper integration boundary
    double offset_max = std::sqrt(10.0*sigma()) * gammalib::deg2rad;
    if (offset_max > gammalib::pi) {
        offset_max = gammalib::pi;
    }

    // Perform numerical integration
    double solidangle = integral.romb(0.0, offset_max) * gammalib::twopi;

    // Return integral
    return solidangle;
}


/***********************************************************************//**
 * @brief Initialise Monte Carlo cache
 *
 * @todo Verify assumption made about the solid angles of the response table
 *       elements.
 * @todo Add optional sampling on a finer spatial grid.
 ***************************************************************************/
void GCTABackgroundPerfTable::init_mc_cache(void) const
{
    // Initialise cache
    m_mc_spectrum.clear();

    // Compute solid angle of model
    double solidangle = this->solidangle();

    // Loop over nodes
    for (int i = 0; i < size(); ++i) {

        // Set energy
        GEnergy energy;
        energy.log10TeV(m_logE[i]);

        // Compute total rate
        #if defined(G_LOG_INTERPOLATION)
        double total_rate = std::pow(10.0, m_background[i]) * solidangle;
        #else
        double total_rate = m_background[i] * solidangle;
        #endif

        // Set node
        if (total_rate > 0.0) {
            m_mc_spectrum.append(energy, total_rate);
        }

    }

    // Return
    return;
}
