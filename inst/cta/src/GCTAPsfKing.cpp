/***************************************************************************
 *      GCTAPsfKing.cpp - King profile CTA point spread function class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Michael Mayer                               *
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
 * @file GCTAPsfKing.hpp
 * @brief King profile CTA point spread function class definition
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GException.hpp"
#include "GRan.hpp"
#include "GFilename.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GCTAPsfKing.hpp"
#include "GCTAException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                               "GCTAPsfKing::read(GFitsTable&)"
#define G_CONTAINMENT_RADIUS       "GCTAPsfKing::containment_radius(double&,"\
                       " double&, double&, double&, double&, double&, bool&)"
#define G_UPDATE                      "GCTAPsfKing::update(double&, double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_FIX_DELTA_MAX
#define G_SMOOTH_PSF

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
const double r_max = 0.7 * gammalib::deg2rad;  // Maximum delta for fixed 
                                               // delta_max computation


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs empty point spread function.
 ***************************************************************************/
GCTAPsfKing::GCTAPsfKing(void) : GCTAPsf()
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
 * Constructs point spread function from a FITS file.
 ***************************************************************************/
GCTAPsfKing::GCTAPsfKing(const GFilename& filename) : GCTAPsf()
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
 *
 * Constructs point spread function by copying from another point spread
 * function.
 ***************************************************************************/
GCTAPsfKing::GCTAPsfKing(const GCTAPsfKing& psf) : GCTAPsf(psf)
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
 *
 * Destructs point spread function.
 ***************************************************************************/
GCTAPsfKing::~GCTAPsfKing(void)
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
 *
 * Assigns point spread function.
 ***************************************************************************/
GCTAPsfKing& GCTAPsfKing::operator=(const GCTAPsfKing& psf)
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
 * @param[in] theta Offset angle in camera system (rad).
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 *
 * Returns the point spread function for a given angular separation in units
 * of sr^-1 for a given energy.
 ***************************************************************************/
double GCTAPsfKing::operator()(const double& delta,
                               const double& logE, 
                               const double& theta, 
                               const double& phi,
                               const double& zenith,
                               const double& azimuth,
                               const bool&   etrue) const
{
    #if defined(G_FIX_DELTA_MAX)
    #if defined(G_SMOOTH_PSF)
    // Set ramp down radius
    static const double ramp_down = 0.95 * r_max;
    static const double norm_down = 1.0 / (r_max - ramp_down);
    #endif
    #endif

    // Initialise PSF value
    double psf = 0.0;

    // Compile option: set PSF to zero outside delta_max
    #if defined(G_FIX_DELTA_MAX)
    if (delta <= r_max) {
    #endif

    // Update the parameter cache
    update(logE, theta);

    // Continue only if normalization is positive
    if (m_par_norm > 0.0) {

        // Compute PSF value
        double arg  = delta / m_par_sigma;
        double arg2 = arg * arg;
        psf = m_par_norm * 
        std::pow((1.0 + 1.0 / (2.0 * m_par_gamma) * arg2), -m_par_gamma);

        // If we are at large offset angles, add a smooth ramp down to
        // avoid steps in the log-likelihood computation
        #if defined(G_FIX_DELTA_MAX)
        #if defined(G_SMOOTH_PSF)
        if (delta > ramp_down) {
            double x = norm_down * (delta - ramp_down);
            psf     *= 1.0 - x * x;
        }
        #endif
        #endif

    } // endif: normalization was positive

    // Compile option: set PSF to zero outside delta_max
    #if defined(G_FIX_DELTA_MAX)
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
 * @brief Clear point spread function
 *
 * Clears point spread function.
 ***************************************************************************/
void GCTAPsfKing::clear(void)
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
 * @brief Clone point spread functions
 *
 * @return Deep copy of point spread function.
 *
 * Returns a pointer to a deep copy of the point spread function.
 ***************************************************************************/
GCTAPsfKing* GCTAPsfKing::clone(void) const
{
    return new GCTAPsfKing(*this);
}


/***********************************************************************//**
 * @brief Read point spread function from FITS table
 *
 * @param[in] table FITS table.
 *
 * @exception GException::invalid_value
 *            Response table is not two-dimensional.
 *
 * Reads the point spread function form the FITS @p table. The following
 * column names are mandatory:
 *
 *     ENERG_LO - Energy lower bin boundaries
 *     ENERG_HI - Energy upper bin boundaries
 *     THETA_LO - Offset angle lower bin boundaries
 *     THETA_HI - Offset angle upper bin boundaries
 *     GAMMA    - Gamma
 *     SIGMA    - Sigma
 *
 * The data are stored in the m_psf member. The energy axis will be set to
 * log10, the offset angle axis to radians.
 ***************************************************************************/
void GCTAPsfKing::read(const GFitsTable& table)
{
    // Clear response table
    m_psf.clear();

    // Read PSF table
    m_psf.read(table);

    // Get mandatory indices (throw exception if not found)
    m_inx_energy = m_psf.axis("ENERG");
    m_inx_theta  = m_psf.axis("THETA");
    m_inx_gamma  = m_psf.table("GAMMA");
    m_inx_sigma  = m_psf.table("SIGMA");

    // Throw an exception if the table is not two-dimensional
    if (m_psf.axes() != 2) {
        std::string msg = "Expected two-dimensional point spread function "
                          "response table but found "+
                          gammalib::str(m_psf.axes())+
                          " dimensions. Please specify a two-dimensional "
                          "point spread function.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Set energy axis to logarithmic scale
    m_psf.axis_log10(m_inx_energy);

    // Set offset angle axis to radians
    m_psf.axis_radians(m_inx_theta);

    // Convert sigma parameters to radians
    m_psf.scale(m_inx_sigma, gammalib::deg2rad);
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write point spread function into FITS table
 *
 * @param[in] table FITS binary table.
 *
 * Writes point spread function into a FITS binary @p table.
 *
 * @todo Add keywords.
 ***************************************************************************/
void GCTAPsfKing::write(GFitsBinTable& table) const
{
    // Create a copy of the response table
    GCTAResponseTable psf(m_psf);

    // Convert sigma parameters back to degrees
    psf.scale(m_inx_sigma, gammalib::rad2deg);

    // Write response table
    psf.write(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load point spread function from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the point spread function from a FITS file.
 *
 * If no extension name is provided, the point spread function will be loaded
 * from the "POINT SPREAD FUNCTION" extension.
 ***************************************************************************/
void GCTAPsfKing::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get PSFa table
    const GFitsTable& table = *fits.table(filename.extname("POINT SPREAD FUNCTION"));

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
 * @brief Save point spread function table into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file? (default: false)
 *
 * Saves point spread function into a FITS file. If a file with the given
 * @p filename does not yet exist it will be created, otherwise the method
 * opens the existing file. The method will create a (or replace an existing)
 * point spread function extension. The extension name can be specified as
 * part of the @p filename, or if no extension name is given, is assumed to
 * be "POINT SPREAD FUNCTION".
 *
 * An existing file will only be modified if the @p clobber flag is set to
 * true.
 ***************************************************************************/
void GCTAPsfKing::save(const GFilename& filename, const bool& clobber) const
{
    // Get extension name
    std::string extname = filename.extname("POINT SPREAD FUNCTION");

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
 * @brief Simulate PSF offset (radians)
 *
 * @param[in] ran Random number generator.
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad).
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 *
 * Draws a random offset for the King Profile
 ***************************************************************************/
double GCTAPsfKing::mc(GRan&         ran,
                       const double& logE, 
                       const double& theta, 
                       const double& phi,
                       const double& zenith,
                       const double& azimuth,
                       const bool&   etrue) const
{
	// Initialise random offset
	double delta = 0.0;

    // Update the parameter cache
    update(logE, theta);

    // Compute exponent
    double exponent = 1.0 / (1.0-m_par_gamma);

    // Compile option: sample until delta <= r_max
    #if defined(G_FIX_DELTA_MAX)
    do {
    #endif

    // Get uniform random number
    double u = ran.uniform();

    // Draw random offset using inversion sampling
    double u_max = (std::pow((1.0 - u), exponent) - 1.0) * m_par_gamma;
    delta = m_par_sigma * std::sqrt(2.0 * u_max);

    // Compile option: sample until delta <= r_max
    #if defined(G_FIX_DELTA_MAX)
    } while (delta > r_max);
    #endif

    // Return PSF offset
    return delta;
}


/***********************************************************************//**
 * @brief Return maximum size of PSF (radians)
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad).
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 *
 * Determine the radius beyond which the PSF becomes negligible. This radius
 * is set by this method to where the containment fraction become 99.995% 
 * which equals \f$5 \times \sigma\f$ of a Gaussian width.
 ***************************************************************************/
double GCTAPsfKing::delta_max(const double& logE,
                              const double& theta, 
                              const double& phi,
                              const double& zenith,
                              const double& azimuth,
                              const bool&   etrue) const
{
    // Compile option: use fixed maximum delta
    #if defined(G_FIX_DELTA_MAX)
    double radius = r_max;
    #else

    // Update the parameter cache
    update(logE, theta);

    // Compute maximum PSF radius (99.995% containment)
    double F      = 0.99995;
    double u_max  = (std::pow((1.0 - F), (1.0/(1.0-m_par_gamma))) - 1.0) * 
                    m_par_gamma;
    double radius = m_par_sigma * std::sqrt(2.0 * u_max);
    #endif

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
 * Calculate the radius from the center that contains 'fraction' percent
 * of the events.  fraction * 100. = Containment %.
 ***************************************************************************/
double GCTAPsfKing::containment_radius(const double& fraction, 
                                       const double& logE, 
                                       const double& theta, 
                                       const double& phi,
                                       const double& zenith,
                                       const double& azimuth,
                                       const bool&   etrue) const
{
    // Check input argument
    if (fraction <= 0.0 || fraction >= 1.0) {
        std::string message = "Containment fraction "+
                              gammalib::str(fraction)+" must be between " +
                              "0.0 and 1.0, not inclusive.";
        throw GException::invalid_argument(G_CONTAINMENT_RADIUS, message);
    }

    // Update the parameter cache
    update(logE, theta);

    // Use analytic calculation
    double arg    = std::pow(1.0 - fraction, 1.0/(1.0-m_par_gamma));
    double radius = m_par_sigma * std::pow(2.0 * m_par_gamma * (arg - 1.0), 0.5);
    
    // Return radius containing fraction of events
    return radius;
}


/***********************************************************************//**
 * @brief Print point spread function information
 *
 * @return Content of point spread function instance.
 * @return String containing point spread function information.
 ***************************************************************************/
std::string GCTAPsfKing::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        double emin = m_psf.axis_lo(m_inx_energy,0);
        double emax = m_psf.axis_hi(m_inx_energy,
                                    m_psf.axis_bins(m_inx_energy)-1);

        // Compute offset angle boundaries in deg
        double omin = m_psf.axis_lo(m_inx_theta,0);
        double omax = m_psf.axis_hi(m_inx_theta,
                                    m_psf.axis_bins(m_inx_theta)-1);

        // Append header
        result.append("=== GCTAPsfKing ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_psf.axis_bins(m_inx_energy)));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(m_psf.axis_bins(m_inx_theta)));
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
void GCTAPsfKing::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_psf.clear();
    m_inx_energy = 0;
    m_inx_theta  = 1;
    m_inx_gamma  = 0;
    m_inx_sigma  = 1;
    m_par_logE   = -1.0e30;
    m_par_theta  = -1.0;
    m_par_norm   = 0.0;
    m_par_gamma  = 0.0;
    m_par_sigma  = 0.0;
    m_par_sigma2 = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GCTAPsfKing::copy_members(const GCTAPsfKing& psf)
{
    // Copy members
    m_filename   = psf.m_filename;
    m_psf        = psf.m_psf;
    m_inx_energy = psf.m_inx_energy;
    m_inx_theta  = psf.m_inx_theta;
    m_inx_gamma  = psf.m_inx_gamma;
    m_inx_sigma  = psf.m_inx_sigma;
    m_par_logE   = psf.m_par_logE;
    m_par_theta  = psf.m_par_theta;
    m_par_norm   = psf.m_par_norm;
    m_par_gamma  = psf.m_par_gamma;
    m_par_sigma  = psf.m_par_sigma;
    m_par_sigma2 = psf.m_par_sigma2;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAPsfKing::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update PSF parameter cache
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle.
 *
 * @exception GException::invalid_value
 *            No valid point spread function information has been found.
 *
 * This method updates the PSF parameter cache. As the performance table PSF
 * only depends on energy, the only parameter on which the cache values
 * depend is the energy.
 ***************************************************************************/
void GCTAPsfKing::update(const double& logE, const double& theta) const
{
    // Only compute PSF parameters if arguments have changed
    if (logE != m_par_logE || theta != m_par_theta) {

        // Save parameters
        m_par_logE  = logE;
        m_par_theta = theta;
    
        // Determine sigma and gamma by interpolating between nodes
        std::vector<double> pars = (m_inx_energy == 0) ? m_psf(logE, theta) :
                                                         m_psf(theta, logE);

        // Set parameters
        m_par_gamma  = pars[m_inx_gamma];
        m_par_sigma  = pars[m_inx_sigma];
        m_par_sigma2 = m_par_sigma * m_par_sigma;

        // Check for parameter sanity
        if (m_par_gamma <= 0.0 || m_par_sigma <= 0.0) {
            m_par_norm = 0.0;
            std::string msg = "King function parameters m_par_gamma and"
                              " m_par_sigma are zero (for "
                              " parameter space logE=" +
                              gammalib::str(logE) + " and theta=" + 
                              gammalib::str(theta) + 
                              "), setting normalization m_norm to zero."; 
            gammalib::warning(G_UPDATE, msg);
        }
        else {   
            // Determine normalisation for given parameters
            m_par_norm = 1.0 / gammalib::twopi * (1.0 - 1.0 / m_par_gamma) /
                         m_par_sigma2;
        }

        // Optionally correct for fixed delta_max
        #if defined(G_FIX_DELTA_MAX)
        double u_max = (r_max*r_max) / (2.0 * m_par_sigma2);
        double norm  = 1.0 - std::pow((1.0 + u_max/m_par_gamma), 1.0-m_par_gamma);
        m_par_norm /= norm;
        #endif

    }

    // Return
    return;
}
