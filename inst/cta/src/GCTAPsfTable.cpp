/***************************************************************************
 *         GCTAPsfTable.cpp - CTA point spread function table class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAPsfTable.cpp
 * @brief CTA point spread function table class implementation
 * @author Juergen Knoedlseder
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
#include "GCTAPsfTable.hpp"
#include "GCTASupport.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                              "GCTAPsfTable::read(GFitsTable&)"
#define G_CONTAINMENT_RADIUS      "GCTAPsfTable::containment_radius(double&,"\
                       " double&, double&, double&, double&, double&, bool&)"

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
 *
 * Constructs empty point spread function.
 ***************************************************************************/
GCTAPsfTable::GCTAPsfTable(void) : GCTAPsf()
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
GCTAPsfTable::GCTAPsfTable(const GFilename& filename) : GCTAPsf()
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
GCTAPsfTable::GCTAPsfTable(const GCTAPsfTable& psf) : GCTAPsf(psf)
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
GCTAPsfTable::~GCTAPsfTable(void)
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
GCTAPsfTable& GCTAPsfTable::operator=(const GCTAPsfTable& psf)
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
 * of sr^-1 for a given energy and offset angle. The PSF value will be
 * determined by trilinear interpolation (and extrapolation) in the PSF
 * histograms. If the interpolation or extrapolation would lead to negative
 * PSF values, the value of zero is returned. A zero value is also returned
 * if the @p delta angle is larger than the largest angle in the response
 * table.
 ***************************************************************************/
double GCTAPsfTable::operator()(const double& delta,
                                const double& logE,
                                const double& theta,
                                const double& phi,
                                const double& zenith,
                                const double& azimuth,
                                const bool&   etrue) const
{
    // Initialise PSF value
    double psf = 0.0;

    // Continue only if delta is not larger than delta_max
    if (delta <= m_delta_max) {

        // Setup argument vector
        double arg[3];
        arg[m_inx_energy] = logE;
        arg[m_inx_theta]  = theta;
        arg[m_inx_delta]  = delta;

        // Compute point spread function by trilinear interpolation
        psf = m_psf(m_inx_rpsf, arg[0], arg[1], arg[2]);

        // Constrain PSF value to non-negative values
        if (psf < 0.0) {
            psf = 0.0;
        }

    } // endif: delta was within validity range

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
void GCTAPsfTable::clear(void)
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
GCTAPsfTable* GCTAPsfTable::clone(void) const
{
    return new GCTAPsfTable(*this);
}


/***********************************************************************//**
 * @brief Read point spread function from FITS table
 *
 * @param[in] table FITS table.
 *
 * @exception GException::invalid_value
 *            Response table is not three-dimensional.
 *
 * Reads the point spread function form the FITS @p table. The following
 * column names are mandatory:
 *
 *     ENERG_LO - Energy lower bin boundaries
 *     ENERG_HI - Energy upper bin boundaries
 *     THETA_LO - Offset angle lower bin boundaries
 *     THETA_HI - Offset angle upper bin boundaries
 *     RAD_LO   - Delta angle lower bin boundaries
 *     RAD_HI   - Delta angle upper bin boundaries
 *     RPSF     - PSF histogram
 *
 * The data are stored in the m_psf member. The energy axis will be set to
 * log10, the offset and delta angle axes to radians.
 ***************************************************************************/
void GCTAPsfTable::read(const GFitsTable& table)
{
    // Clear response table
    m_psf.clear();

    // Read PSF table
    m_psf.read(table);

    // Get mandatory indices (throw exception if not found)
    m_inx_energy = m_psf.axis("ENERG");
    m_inx_theta  = m_psf.axis("THETA");
    m_inx_delta  = m_psf.axis("RAD");
    m_inx_rpsf   = m_psf.table("RPSF");

    // Throw an exception if the table is not three-dimensional
    if (m_psf.axes() != 3) {
        std::string msg = "Expected three-dimensional point spread function "
                          "response table but found "+
                          gammalib::str(m_psf.axes())+
                          " dimensions. Please specify a three-dimensional "
                          "point spread function.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Set energy axis to logarithmic scale
    m_psf.axis_log10(m_inx_energy);

    // Set offset angle axis to radians
    m_psf.axis_radians(m_inx_theta);

    // Set delta angle axis to radians
    m_psf.axis_radians(m_inx_delta);

    // Precompute PSF information
    precompute();

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
void GCTAPsfTable::write(GFitsBinTable& table) const
{
    // Create a copy of the response table
    GCTAResponseTable psf(m_psf);

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
 * If no extension name is given the method scans the `HDUCLASS` keywords
 * of all extensions and loads the background from the first extension
 * for which `HDUCLAS4=PSF_TABLE`.
 *
 * Otherwise, the background will be loaded from the `POINT SPREAD FUNCTION`
 * extension.
 ***************************************************************************/
void GCTAPsfTable::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get the default extension name. If no GADF compliant name was found
    // then set the default extension name to "POINT SPREAD FUNCTION".
    std::string extname = gammalib::gadf_hduclas4(fits, "PSF_TABLE");
    if (extname.empty()) {
        extname = gammalib::extname_cta_psftable;
    }

    // Get PSF table
    const GFitsTable& table = *fits.table(filename.extname(extname));

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
 * @param[in] clobber Overwrite existing file?
 *
 * Saves point spread function into a FITS file. If a file with the given
 * @p filename does not yet exist it will be created, otherwise the method
 * opens the existing file. The method will create a (or replace an existing)
 * point spread function extension. The extension name can be specified as
 * part of the @p filename, or if no extension name is given, is assumed to
 * be `POINT SPREAD FUNCTION`.
 *
 * An existing file will only be modified if the @p clobber flag is set to
 * true.
 ***************************************************************************/
void GCTAPsfTable::save(const GFilename& filename, const bool& clobber) const
{
    // Get extension name
    std::string extname = filename.extname(gammalib::extname_cta_psftable);

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
 * Draws a random offset from the point spread function using a rejection
 * method.
 ***************************************************************************/
double GCTAPsfTable::mc(GRan&         ran,
                        const double& logE,
                        const double& theta,
                        const double& phi,
                        const double& zenith,
                        const double& azimuth,
                        const bool&   etrue) const
{
	// Initialise random offset
	double delta = 0.0;

    // Simulate PSF only if maximum PSF radius is positive
    if (m_delta_max > 0.0 && m_psf_max > 0.0) {

        // Pre-compute 1-cos(delta_max)
        double cosrad = 1.0 - std::cos(m_delta_max);

        // Draw random positions
        while (true) {
    
            // Simulate offset angle
            delta = std::acos(1.0 - ran.uniform() * cosrad);

            // Get corresponding PSF value
            double value = this->operator()(delta, logE, theta, phi, zenith,
                                            azimuth, etrue);

            // Get uniform random number up to the maximum
            double uniform = ran.uniform() * m_psf_max;

            // Exit loop if the random number is not larger than the PSF value
            if (uniform <= value) {
                break;
            }

        } // endwhile: rejection method

    } // endif: PSF radius and maximum PSF value was positive

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
 * Returns the maximum PSF radius.
 ***************************************************************************/
double GCTAPsfTable::delta_max(const double& logE,
                               const double& theta,
                               const double& phi,
                               const double& zenith,
                               const double& azimuth,
                               const bool&   etrue) const
{
    // Return maximum PSF radius
    return m_delta_max;
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
double GCTAPsfTable::containment_radius(const double& fraction, 
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

    // Determine number of delta bins
    int ndelta = m_psf.axis_bins(m_inx_delta);

    // Initialise containment radius and PSF integral
    double delta = 0.0;
    double sum   = 0.0;

    // Integrate delta array by computing the sum over the pixels times the
    // pixel size times the delta angle
    for (int idelta = 0; idelta < ndelta; ++idelta) {

        // Compute delta value (radians)
        delta = m_psf.axis_nodes(m_inx_delta)[idelta];

        // Compute delta bin width (radians)
        double width = (m_psf.axis_hi(m_inx_delta, idelta) -
                        m_psf.axis_lo(m_inx_delta, idelta)) *
                        gammalib::deg2rad;

        // Integrate PSF
        sum += this->operator()(delta, logE, theta, phi, zenith, azimuth, etrue) *
               std::sin(delta) * width * gammalib::twopi;

        // Break of the containment fraction is reached or exceeded
        if (sum >= fraction) {
            break;
        }

    } // endfor: loop over delta array
    
    // Return containing radius
    return delta;
}


/***********************************************************************//**
 * @brief Print point spread function information
 *
 * @param[in] chatter Chattiness.
 * @return String containing point spread function information.
 ***************************************************************************/
std::string GCTAPsfTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAPsfTable ===");
        result.append("\n"+gammalib::parformat("Filename")+m_filename);

        // Initialise information
        int    nebins     = 0;
        int    nthetabins = 0;
        int    ndeltabins = 0;
        double emin       = 0.0;
        double emax       = 0.0;
        double omin       = 0.0;
        double omax       = 0.0;
        double dmin       = 0.0;
        double dmax       = 0.0;

        // Extract information if there are axes in the response table
        if (m_psf.axes() > 0) {
            nebins     = m_psf.axis_bins(m_inx_energy);
            nthetabins = m_psf.axis_bins(m_inx_theta);
            ndeltabins = m_psf.axis_bins(m_inx_delta);
            emin       = m_psf.axis_lo(m_inx_energy,0);
            emax       = m_psf.axis_hi(m_inx_energy,nebins-1);
            omin       = m_psf.axis_lo(m_inx_theta,0);
            omax       = m_psf.axis_hi(m_inx_theta,nthetabins-1);
            dmin       = m_psf.axis_lo(m_inx_delta,0);
            dmax       = m_psf.axis_hi(m_inx_delta,ndeltabins-1);
        }


        // Append information
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(nebins));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(nthetabins));
        result.append("\n"+gammalib::parformat("Number of delta bins") +
                      gammalib::str(ndeltabins));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");
        result.append("\n"+gammalib::parformat("Offset angle range"));
        result.append(gammalib::str(omin)+" - "+gammalib::str(omax)+" deg");
        result.append("\n"+gammalib::parformat("Delta angle range"));
        result.append(gammalib::str(dmin)+" - "+gammalib::str(dmax)+" deg");

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
void GCTAPsfTable::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_psf.clear();
    m_inx_energy = 0;
    m_inx_theta  = 1;
    m_inx_delta  = 2;
    m_inx_rpsf   = 0;
    m_delta_max  = 0.0;
    m_psf_max    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GCTAPsfTable::copy_members(const GCTAPsfTable& psf)
{
    // Copy members
    m_filename   = psf.m_filename;
    m_psf        = psf.m_psf;
    m_inx_energy = psf.m_inx_energy;
    m_inx_theta  = psf.m_inx_theta;
    m_inx_delta  = psf.m_inx_delta;
    m_inx_rpsf   = psf.m_inx_rpsf;
    m_delta_max  = psf.m_delta_max;
    m_psf_max    = psf.m_psf_max;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAPsfTable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Performs precomputations for point spread function
 *
 * Replaces any invalid PSF histogram by zeroes, normalises the PSF
 * histograms to unity, and computes maximum PSF value and maximum delta
 * value.
 ***************************************************************************/
void GCTAPsfTable::precompute(void)
{
    // Determine PSF dimensions
    int neng   = m_psf.axis_bins(m_inx_energy);
    int ntheta = m_psf.axis_bins(m_inx_theta);
    int ndelta = m_psf.axis_bins(m_inx_delta);

    // Determine delta element increment
    int inx_inc = element(0,0,1) - element(0,0,0);

    // Replace any negative or invalid histogram values by zero values
    for (int i = 0; i < m_psf.elements(); ++i) {
        double element = m_psf(m_inx_rpsf, i);
        if (element < 0.0 ||
            gammalib::is_notanumber(element) ||
            gammalib::is_infinite(element)) {
            m_psf(m_inx_rpsf, i) = 0.0;
        }
    }

    // Normalise PSF to unity
    for (int ieng = 0; ieng < neng; ++ieng) {
        for (int itheta = 0; itheta < ntheta; ++itheta) {

            // Initialise PSF integral
            double sum = 0.0;

            // Set start element
            int inx = element(ieng, itheta, 0);

            // Integrate delta array by computing the sum over the pixels
            // times the pixel size time the delta angle
            for (int idelta = 0; idelta < ndelta; ++idelta, inx += inx_inc) {

                // Compute delta value (radians)
                double delta = m_psf.axis_nodes(m_inx_delta)[idelta];

                // Compute delta bin width (radians)
                double width = (m_psf.axis_hi(m_inx_delta, idelta) -
                                m_psf.axis_lo(m_inx_delta, idelta)) *
                                gammalib::deg2rad;

                // Integrate PSF
                sum += m_psf(m_inx_rpsf, inx) * std::sin(delta) * width *
                       gammalib::twopi;

            }

            // If integral is positive then divide PSF by integral so that it
            // normalises to unity
            if (sum > 0.0) {

                // Set start element
                int inx = element(ieng, itheta, 0);

                // Normalise PSF
                for (int idelta = 0; idelta < ndelta; ++idelta, inx += inx_inc) {
                    m_psf(m_inx_rpsf, inx) /= sum;
                }

            } // endif: PSF integral was positive

        } // endfor: looped over all theta angles
    } // endfor: looped over all PSF values

    // Compute maximum PSF value
    m_psf_max = 0.0;
    for (int i = 0; i < m_psf.elements(); ++i) {
        if (m_psf(m_inx_rpsf, i) > m_psf_max) {
            m_psf_max = m_psf(m_inx_rpsf, i);
        }
    }

    // Set maximum PSF radius (radians)
    m_delta_max = m_psf.axis_hi(m_inx_delta, m_psf.axis_bins(m_inx_delta)-1) *
                  gammalib::deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return element index
 *
 * @param[in] ieng Energy index [0,...,neng-1]
 * @param[in] itheta Offset angle index [0,...,ntheta-1]
 * @param[in] idelta Delta angle index [0,...,ndelta-1]
 ***************************************************************************/
int GCTAPsfTable::element(const int& ieng, const int& itheta, const int& idelta)
{
    // Setup index vector for requested element
    int index[3];
    index[m_inx_energy] = ieng;
    index[m_inx_theta]  = itheta;
    index[m_inx_delta]  = idelta;
    
    // Get element index
    int element= index[0] + (index[1] + index[2]*m_psf.axis_bins(1)) *
                             m_psf.axis_bins(0);

    // Return element index
    return element;
}
