/***************************************************************************
 *       GCTAPsfVector.cpp - CTA point spread function vector class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAPsfVector.cpp
 * @brief CTA point spread function vector class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GFilename.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"
#include "GCTAPsfVector.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONTAINMENT_RADIUS     "GCTAPsfVector::containment_radius(double&,"\
                       " double&, double&, double&, double&, double&, bool&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_SMOOTH_PSF

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
GCTAPsfVector::GCTAPsfVector(void) : GCTAPsf()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename PSF FITS file.
 *
 * Constructs point spread function vector from a FITS file.
 ***************************************************************************/
GCTAPsfVector::GCTAPsfVector(const GFilename& filename) : GCTAPsf()
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
 * @param[in] psf Point spread function.
 ***************************************************************************/
GCTAPsfVector::GCTAPsfVector(const GCTAPsfVector& psf) : GCTAPsf(psf)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(psf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAPsfVector::~GCTAPsfVector(void)
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
 * @param[in] psf Point spread function.
 * @return Point spread function.
 ***************************************************************************/
GCTAPsfVector& GCTAPsfVector::operator=(const GCTAPsfVector& psf)
{
    // Execute only if object is not identical
    if (this != &psf) {

        // Copy base class members
        this->GCTAPsf::operator=(psf);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(psf);

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
double GCTAPsfVector::operator()(const double& delta,
                                 const double& logE, 
                                 const double& theta, 
                                 const double& phi,
                                 const double& zenith,
                                 const double& azimuth,
                                 const bool&   etrue) const
{
    #if defined(G_SMOOTH_PSF)
    // Compute offset so that PSF goes to 0 at 5 times the sigma value. This
    // is a kluge to get a PSF that smoothly goes to zero at the edge, which
    // prevents steps or kinks in the log-likelihood function.
    static const double offset = std::exp(-0.5*5.0*5.0);
    #endif

    // Update the parameter cache
    update(logE);

    // Compute PSF value
    #if defined(G_SMOOTH_PSF)
    double psf = m_par_scale * (std::exp(m_par_width * delta * delta) - offset);
    #else
    double psf = m_par_scale * (std::exp(m_par_width * delta * delta));
    #endif

    #if defined(G_SMOOTH_PSF)
    // Make sure that PSF is non-negative
    if (psf < 0.0) {
        psf = 0.0;
    }
    #endif
    
    // Return PSF
    return psf;
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
void GCTAPsfVector::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAPsf::free_members();

    // Initialise members
    this->GCTAPsf::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of point spread function instance.
 ***************************************************************************/
GCTAPsfVector* GCTAPsfVector::clone(void) const
{
    return new GCTAPsfVector(*this);
}


/***********************************************************************//**
 * @brief Load point spread function from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the point spread function from a FITS file.
 *
 * If no extension name is provided, the point spread function will be loaded
 * from the "PSF" extension.
 ***************************************************************************/
void GCTAPsfVector::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get PSFa table
    const GFitsTable& table = *fits.table(filename.extname("PSF"));

    // Read PSF from table
    read(table);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read point spread function vector from FITS table
 *
 * @param[in] table FITS table.
 *
 * Reads a point spread function vector from the FITS @p table.
 *
 * The energies are converted to TeV. Conversion is done based on the units
 * provided for the energy columns. Units that are recognized are 'keV',
 * 'MeV', 'GeV', and 'TeV' (case independent).
 *
 * The Gaussian width parameter may be either given in 68% containment radius
 * (R68) or in sigma (ANGRES40; corresponding to a 38% containment radius).
 * The former format is used for CTA and H.E.S.S. data, the latter for
 * MAGIC data. All these things should be more uniform once we have a well
 * defined format.
 ***************************************************************************/
void GCTAPsfVector::read(const GFitsTable& table)
{
    // Clear arrays
    m_logE.clear();
    m_r68.clear();
    m_sigma.clear();

    // Set conversion factor from 68% containment radius to 1 sigma
    const double conv = 0.6624305 * gammalib::deg2rad;

    // Get pointers to table columns
    const GFitsTableCol* energy_lo = table["ENERG_LO"];
    const GFitsTableCol* energy_hi = table["ENERG_HI"];

    // Handle various data formats (H.E.S.S. and MAGIC)
    const GFitsTableCol* r68;
    double               r68_scale = 1.0;
    if (table.contains("R68")) {
        r68 = table["R68"];
    }
    else {
        r68 = table["ANGRES40"];
        r68_scale = 1.0 / 0.6624305; // MAGIC PSF is already 1 sigma
    }

    // Determine unit conversion factors (default: TeV)
    std::string u_energy_lo = gammalib::tolower(gammalib::strip_whitespace(energy_lo->unit()));
    std::string u_energy_hi = gammalib::tolower(gammalib::strip_whitespace(energy_hi->unit()));
    double c_energy_lo = 1.0;
    double c_energy_hi = 1.0;
    if (u_energy_lo == "kev") {
        c_energy_lo = 1.0e-9;
    }
    else if (u_energy_lo == "mev") {
        c_energy_lo = 1.0e-6;
    }
    else if (u_energy_lo == "gev") {
        c_energy_lo = 1.0e-3;
    }
    if (u_energy_hi == "kev") {
        c_energy_hi = 1.0e-9;
    }
    else if (u_energy_hi == "mev") {
        c_energy_hi = 1.0e-6;
    }
    else if (u_energy_hi == "gev") {
        c_energy_hi = 1.0e-3;
    }

    // Extract number of energy bins
    int num = energy_lo->nrows();

    // Set nodes
    for (int i = 0; i < num; ++i) {
    
        // Compute log10 mean energy in TeV
        double e_min = energy_lo->real(i) * c_energy_lo;
        double e_max = energy_hi->real(i) * c_energy_hi;
        double logE  = 0.5 * (log10(e_min) + log10(e_max));
        
        // Extract r68 value and scale as required
        double r68_value = r68->real(i) * r68_scale;
        
        // Store log10 mean energy and r68 value
        m_logE.append(logE);
        m_r68.push_back(r68_value);
        m_sigma.push_back(r68_value*conv);

    } // endfor: looped over nodes
    
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
double GCTAPsfVector::mc(GRan&         ran,
                         const double& logE, 
                         const double& theta, 
                         const double& phi,
                         const double& zenith,
                         const double& azimuth,
                         const bool&   etrue) const
{
    // Update the parameter cache
    update(logE);

    // Draw offset
    double delta = m_par_sigma * ran.chisq2();
    
    // Return PSF offset
    return delta;
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
double GCTAPsfVector::delta_max(const double& logE, 
                                   const double& theta, 
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth,
                                   const bool&   etrue) const
{
    // Update the parameter cache
    update(logE);

    // Compute maximum PSF radius
    double radius = 5.0 * m_par_sigma;
    
    // Return maximum PSF radius
    return radius;
}


/***********************************************************************//**
 * @brief Return the radius that contains a fraction of the events (radians)
 *
 * @param[in] fraction of events (0.0-1.0)
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 * @return Containment radius (radians).
 *
 * @exception GException::invalid_argument
 *            Invalid fraction specified.
 *
 * Evaluates:
 * 
 * \f[
 *
 * radius = \sqrt{ \frac{\ln{\left( \frac{ fraction * m\_par\_width}
 *          {\pi*m\_par\_scale} +1\right)}}{m\_par\_width} }
 *
 * \f]
 *
 * which is derived by integrating
 *
 * \f[
 *
 * fraction = \int_{0}^{2\pi}\int_{0}^{radius}r * 
 *            e^{ m\_par\_width * r^{2}}dr d\phi
 *
 * \f]
 *
 * and solving for radius.
 *  
 * Calculate the radius from the center that contains 'fraction' percent
 * of the events.  fraction * 100. = Containment % .
 ***************************************************************************/
double GCTAPsfVector::containment_radius(const double& fraction, 
                                         const double& logE, 
                                         const double& theta, 
                                         const double& phi,
                                         const double& zenith,
                                         const double& azimuth,
                                         const bool&   etrue) const
{
    // Check input argument
    if (fraction <= 0.0 || fraction >= 1.0) {
        std::string msg = "Containment fraction "+gammalib::str(fraction)+
                          " must lie in interval ]0,1[. Please specify a "
                          "valid containment fraction";
        throw GException::invalid_argument(G_CONTAINMENT_RADIUS, msg);
    }

    // Update the parameter cache
    update(logE);

    // Compute radius for gaussian
    double arg    = fraction * m_par_width / (gammalib::pi * m_par_scale);
    double radius = std::sqrt(std::log(arg + 1.0) / m_par_width);
    
    // Return containement radius
    return radius;
}


/***********************************************************************//**
 * @brief Print point spread function information
 *
 * @param[in] chatter Chattiness.
 * @return String containing point spread function information.
 ***************************************************************************/
std::string GCTAPsfVector::print(const GChatter& chatter) const
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
        result.append("=== GCTAPsfVector ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(num));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");

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
void GCTAPsfVector::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_logE.clear();
    m_r68.clear();
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
void GCTAPsfVector::copy_members(const GCTAPsfVector& psf)
{
    // Copy members
    m_filename  = psf.m_filename;
    m_logE      = psf.m_logE;
    m_r68       = psf.m_r68;
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
void GCTAPsfVector::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update PSF parameter cache
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 *
 * This method updates the PSF parameter cache. As the performance table PSF
 * only depends on energy, the only parameter on which the cache values
 * depend is the energy.
 ***************************************************************************/
void GCTAPsfVector::update(const double& logE) const
{
    // Only compute PSF parameters if arguments have changed
    if (logE != m_par_logE) {

        // Save energy
        m_par_logE = logE;
    
        // Determine Gaussian sigma in radians
        m_par_sigma = m_logE.interpolate(logE, m_sigma);

        // Derive width=-0.5/(sigma*sigma) and scale=1/(twopi*sigma*sigma)
        double sigma2 = m_par_sigma * m_par_sigma;
        m_par_scale   =  1.0 / (gammalib::twopi * sigma2);
        m_par_width   = -0.5 / sigma2;

    }

    // Return
    return;
}
