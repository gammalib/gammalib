/***************************************************************************
 *           GCTAPsf2D.cpp - CTA 2D point spread function class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAPsf2D.hpp
 * @brief CTA 2D point spread function class implementation
 * @author Juergen Knoedlseder
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
#include "GCTAPsf2D.hpp"
#include "GCTAException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                                "GCTAPsf2D::load(std::string&)"

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
GCTAPsf2D::GCTAPsf2D(void) : GCTAPsf()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename PSF FITS file name.
 *
 * Construct instance by loading the point spread function information from
 * a PSF response table.
 ***************************************************************************/
GCTAPsf2D::GCTAPsf2D(const std::string& filename) : GCTAPsf()
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
GCTAPsf2D::GCTAPsf2D(const GCTAPsf2D& psf) : GCTAPsf(psf)
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
GCTAPsf2D::~GCTAPsf2D(void)
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
GCTAPsf2D& GCTAPsf2D::operator=(const GCTAPsf2D& psf)
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
 * @param[in] theta Offset angle in camera system (rad). Defaults to 0.0.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 *
 * Returns the point spread function for a given angular separation in units
 * of sr^-1 for a given energy and offset angle.
 ***************************************************************************/
double GCTAPsf2D::operator()(const double& delta,
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

    // Initialise PSF value
    double psf = 0.0;

    // Update the parameter cache
    update(logE, theta);

    // Continue only if normalization is positive
    if (m_norm > 0.0) {

        // Compute distance squared
        double delta2 = delta * delta;

        // Compute Psf value
        #if defined(G_SMOOTH_PSF)
        psf = std::exp(m_width1 * delta2) - offset;
        if (m_norm2 > 0.0) {
            psf += (std::exp(m_width2 * delta2)-offset) * m_norm2;
        }
        if (m_norm3 > 0.0) {
            psf += (std::exp(m_width3 * delta2)-offset) * m_norm3;
        }
        #else
        psf = std::exp(m_width1 * delta2);
        if (m_norm2 > 0.0) {
            psf += std::exp(m_width2 * delta2) * m_norm2;
        }
        if (m_norm3 > 0.0) {
            psf += std::exp(m_width3 * delta2) * m_norm3;
        }
        #endif
        psf *= m_norm;

        #if defined(G_SMOOTH_PSF)
        // Make sure that PSF is non-negative
        if (psf < 0.0) {
            psf = 0.0;
        }
        #endif

    } // endif: normalization was positive
    
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
void GCTAPsf2D::clear(void)
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
GCTAPsf2D* GCTAPsf2D::clone(void) const
{
    return new GCTAPsf2D(*this);
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
void GCTAPsf2D::read(const GFits& fits)
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
    m_psf.scale(3, gammalib::deg2rad);
    m_psf.scale(5, gammalib::deg2rad);

    // Return
    return;
}





/***********************************************************************//**
 * @brief Load point spread function from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the point spread function information from a PSF
 * response table.
 ***************************************************************************/
void GCTAPsf2D::load(const std::string& filename)
{
    // Open PSF FITS file
    GFits fits(filename);

    // Read fits file
    read(fits);

    // Close PSF FITS file
    fits.close();

    // Store filename
    m_filename = filename;

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
void GCTAPsf2D::write(GFitsBinTable& hdu) const
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
void GCTAPsf2D::save(const std::string& filename, const bool& clobber) const
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
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 *
 * Draws a random offset for a three component Gaussian PSF.
 ***************************************************************************/
double GCTAPsf2D::mc(GRan&         ran,
                     const double& logE,
                     const double& theta,
                     const double& phi,
                     const double& zenith,
                     const double& azimuth,
                     const bool&   etrue) const
{
    // Update the parameter cache
    update(logE, theta);

    // Select in which Gaussian we are
    double sigma = m_sigma1;
    double sum1  = m_sigma1;
    double sum2  = m_sigma2 * m_norm2;
    double sum3  = m_sigma3 * m_norm3;
    double sum   = sum1 + sum2 + sum3;
    double u     = ran.uniform() * sum;
    if (sum2 > 0.0 && u >= sum2) {
        sigma = m_sigma3;
    }
    else if (sum1 > 0.0 && u >= sum1) {
        sigma = m_sigma2;
    }

    // Now draw from the selected Gaussian
    double delta = sigma * ran.chisq2();
    
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
double GCTAPsf2D::delta_max(const double& logE, 
                            const double& theta, 
                            const double& phi,
                            const double& zenith,
                            const double& azimuth,
                            const bool&   etrue) const
{
    // Update the parameter cache
    update(logE, theta);

    // Compute maximum sigma
    double sigma = m_sigma1;
    if (m_sigma2 > sigma) sigma = m_sigma2;
    if (m_sigma3 > sigma) sigma = m_sigma3;

    // Compute maximum PSF radius
    double radius = 5.0 * sigma;
    
    // Return maximum PSF radius
    return radius;
}


/***********************************************************************//**
 * @brief Print point spread function information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing point spread function information.
 ***************************************************************************/
std::string GCTAPsf2D::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        double emin = m_psf.axis_lo(0,0);
        double emax = m_psf.axis_hi(0,m_psf.axis(0)-1);

        // Compute offset angle boundaries in deg
        double omin = m_psf.axis_lo(1,0);
        double omax = m_psf.axis_hi(1,m_psf.axis(1)-1);

        // Append header
        result.append("=== GCTAPsf2D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_psf.axis(0)));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(m_psf.axis(1)));
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
void GCTAPsf2D::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_psf.clear();
    m_par_logE  = -1.0e30;
    m_par_theta = -1.0;
    m_norm      = 1.0;
    m_norm2     = 0.0;
    m_norm3     = 0.0;
    m_sigma1    = 0.0;
    m_sigma2    = 0.0;
    m_sigma3    = 0.0;
    m_width1    = 0.0;
    m_width2    = 0.0;
    m_width3    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GCTAPsf2D::copy_members(const GCTAPsf2D& psf)
{
    // Copy members
    m_filename  = psf.m_filename;
    m_psf       = psf.m_psf;
    m_par_logE  = psf.m_par_logE;
    m_par_theta = psf.m_par_theta;
    m_norm      = psf.m_norm;
    m_norm2     = psf.m_norm2;
    m_norm3     = psf.m_norm3;
    m_sigma1    = psf.m_sigma1;
    m_sigma2    = psf.m_sigma2;
    m_sigma3    = psf.m_sigma3;
    m_width1    = psf.m_width1;
    m_width2    = psf.m_width2;
    m_width3    = psf.m_width3;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAPsf2D::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update PSF parameter cache
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad).
 *
 * This method updates the PSF parameter cache.
 ***************************************************************************/
void GCTAPsf2D::update(const double& logE, const double& theta) const
{
    // Only compute PSF parameters if arguments have changed
    if (logE != m_par_logE || theta != m_par_theta) {

        // Save parameters
        m_par_logE  = logE;
        m_par_theta = theta;

        // Interpolate response parameters
        std::vector<double> pars = m_psf(logE, theta);

        // Set Gaussian sigmas
        m_sigma1 = pars[1];
        m_sigma2 = pars[3];
        m_sigma3 = pars[5];

        // Set width parameters
        double sigma1 = m_sigma1 * m_sigma1;
        double sigma2 = m_sigma2 * m_sigma2;
        double sigma3 = m_sigma3 * m_sigma3;

        // Compute Gaussian 1
        if (sigma1 > 0.0) {
            m_width1 = -0.5 / sigma1;
        }
        else {
            m_width1 = 0.0;
        }

        // Compute Gaussian 2
        if (sigma2 > 0.0) {
            m_width2 = -0.5 / sigma2;
            m_norm2  = pars[2];
        }
        else {
            m_width2 = 0.0;
            m_norm2  = 0.0;
        }

        // Compute Gaussian 3
        if (sigma3 > 0.0) {
            m_width3 = -0.5 / sigma3;
            m_norm3  = pars[4];
        }
        else {
            m_width3 = 0.0;
            m_norm3  = 0.0;
        }

        // Compute global normalization parameter
        double integral = gammalib::twopi * (sigma1 + sigma2*m_norm2 + sigma3*m_norm3);
        m_norm = (integral > 0.0) ? 1.0 / integral : 0.0;

    }

    // Return
    return;
}
