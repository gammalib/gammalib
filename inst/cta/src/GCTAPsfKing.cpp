/***************************************************************************
 *      GCTAPsfKing.cpp - King profile CTA point spread function class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Michael Mayer                               *
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
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GCTAPsfKing.hpp"
#include "GCTAException.hpp"

/* __ Method name definitions ____________________________________________ */
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
 * Construct instance by loading the point spread function information from
 * a FITS file that contains the PSF information in form of a column.
 ***************************************************************************/
GCTAPsfKing::GCTAPsfKing(const std::string& filename) : GCTAPsf()
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
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
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
 * @brief Clone instance
 *
 * @return Deep copy of point spread function instance.
 ***************************************************************************/
GCTAPsfKing* GCTAPsfKing::clone(void) const
{
    return new GCTAPsfKing(*this);
}


/***********************************************************************//**
 * @brief Load point spread function from binary table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the point spread function information from a PSF
 * response table.
 ***************************************************************************/
void GCTAPsfKing::load(const std::string& filename)
{
    // Open PSF FITS file
    GFits file(filename);

    // Read fits file
    read(file);

    // Close PSF FITS file
    file.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Read PSF from FITS file
 *
 * @param[in] fits FITS file pointer.
 *
 * Reads the PSF from the FITS file extension "POINT SPREAD FUNCTION". The data
 * are stored in m_psf which is of type GCTAResponseTable. 
 * The energy axis will be set to log10. The offset angle axis and 
 * sigma parameter columns will be set to radians.
 ***************************************************************************/
void GCTAPsfKing::read(const GFits& fits)
{
    // Clear response table
    m_psf.clear();

    // Get PSF table
    const GFitsTable& table = *fits.table("POINT SPREAD FUNCTION");

    // Read PSF table
    m_psf.read(table);

    // Set energy axis to logarithmic scale
    m_psf.axis_log10(0);

    // Set offset angle axis to radians
    m_psf.axis_radians(1);

    // Convert sigma parameters to radians
    m_psf.scale(1, gammalib::deg2rad);
   
    // Return
    return;
}

/***********************************************************************//**
 * @brief Write CTA PSF table into FITS binary table object.
 *
 * @param[in] hdu FITS binary table.
 *
 * @todo Add necessary keywords.
 ***************************************************************************/
void GCTAPsfKing::write(GFitsBinTable& hdu) const
{
    // Write background table
    m_psf.write(hdu);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Save PSF table into FITS file
 *
 * @param[in] filename PSF table FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * Save the PSF table into a FITS file.
 ***************************************************************************/
void GCTAPsfKing::save(const std::string& filename, const bool& clobber) const
{
    // Create binary table
    GFitsBinTable table;
    table.extname("POINT SPREAD FUNCTION");

    // Write the PSF table
    write(table);

    // Create FITS file, append table, and write into the file
    GFits fits;
    fits.append(table);
    fits.saveto(filename, clobber);

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
 * @brief Print point spread function information
 *
 * @return Content of point spread function instance.
 ***************************************************************************/
std::string GCTAPsfKing::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    if (chatter != SILENT) {

		// Append information
		result.append("=== GCTAPsfKing ===");
		result.append("\n"+gammalib::parformat("Filename")+m_filename);
		result.append("\n"+m_psf.print(chatter));
    }

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
        std::vector<double> pars = m_psf(logE,theta);

        // Throw an exception if there are not 2 parameters
        if (pars.size() != 2) {
            std::string msg = gammalib::str(pars.size()) + " parameters have"
                              " been found in the response table of the"
                              " King profile response function while 2"
                              " parameters are expected.\n"
                              "Possibly, the point spread function information"
                              " has not yet been loaded. Please load the point"
                              " spread function before using it.";
            throw GException::invalid_value(G_UPDATE, msg);
        }

        // Set parameters
        m_par_gamma  = pars[0];
        m_par_sigma  = pars[1];
        m_par_sigma2 = m_par_sigma * m_par_sigma;

        // Determine normalisation for given parameters
        m_par_norm = 1.0 / gammalib::twopi * (1.0 - 1.0 / m_par_gamma) /
                     m_par_sigma2;

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
