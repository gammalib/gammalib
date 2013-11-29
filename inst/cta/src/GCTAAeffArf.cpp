/***************************************************************************
 *                GCTAAeffArf.cpp - CTA ARF effective area class           *
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
 * @file GCTAAeffArf.hpp
 * @brief CTA ARF effective area class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GIntegral.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"
#include "GCTAAeffArf.hpp"
#include "GCTAResponse.hpp"
#include "GCTAResponse_helpers.hpp"
#include "GCTAException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                              "GCTAAeffArf::load(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_APPLY_THETACUT              //!< Debug thetacut application

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAAeffArf::GCTAAeffArf(void) : GCTAAeff()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename ARF FITS file name.
 *
 * Construct instance by loading the effective area information from an
 * ARF FITS file.
 ***************************************************************************/
GCTAAeffArf::GCTAAeffArf(const std::string&  filename) : GCTAAeff()
{
    // Initialise class members
    init_members();

    // Load ARF
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] aeff Effective area.
 ***************************************************************************/
GCTAAeffArf::GCTAAeffArf(const GCTAAeffArf& aeff) : GCTAAeff(aeff)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(aeff);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAAeffArf::~GCTAAeffArf(void)
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
 * @param[in] aeff Effective area.
 * @return Effective area.
 ***************************************************************************/
GCTAAeffArf& GCTAAeffArf::operator= (const GCTAAeffArf& aeff)
{
    // Execute only if object is not identical
    if (this != &aeff) {

        // Copy base class members
        this->GCTAAeff::operator=(aeff);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(aeff);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return effective area in units of cm2
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Defaults to 0.0.
 * @param[in] phi Azimuth angle in camera system (rad). Not used in this method.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used in this method.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used in this method.
 * @param[in] etrue Use true energy (true/false). Not used.
 *
 * Returns the effective area in units of cm2 for a given energy and
 * offset angle. The effective area is linearily interpolated in
 * log10(energy). The method assures that the effective area value never
 * becomes negative.
 ***************************************************************************/
double GCTAAeffArf::operator()(const double& logE, 
                               const double& theta, 
                               const double& phi,
                               const double& zenith,
                               const double& azimuth,
                               const bool&   etrue) const
{
    // Get effective area value in cm2
    double aeff = m_logE.interpolate(logE, m_aeff);

    // Make sure that effective area is not negative
    if (aeff < 0.0) {
        aeff = 0.0;
    }

    // Optionally add in Gaussian offset angle dependence
    if (m_sigma != 0.0) {
        double offset = theta * gammalib::rad2deg;
        double arg    = offset * offset / m_sigma;
        double scale  = exp(-0.5 * arg * arg);
        aeff         *= scale;
    }
    
    // Return effective area value
    return aeff;
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
void GCTAAeffArf::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAAeff::free_members();

    // Initialise members
    this->GCTAAeff::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of effective area instance.
 ***************************************************************************/
GCTAAeffArf* GCTAAeffArf::clone(void) const
{
    return new GCTAAeffArf(*this);
}


/***********************************************************************//**
 * @brief Load effective area from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the effective area information from an ASCII
 * performance table.
 ***************************************************************************/
void GCTAAeffArf::load(const std::string& filename)
{
    // Open ARF FITS file
    GFits file(filename);

    // Get ARF table
    GFitsTable* table = file.table("SPECRESP");

    // Read ARF
    read_arf(table);

    // Close ARF FITS file
    file.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which effective area was loaded
 ***************************************************************************/
std::string GCTAAeffArf::filename(void) const
{
    // Return filename
    return m_filename;
}


/***********************************************************************//**
 * @brief Read CTA ARF vector
 *
 * @param[in] hdu FITS table pointer.
 *
 * This method reads a CTA ARF vector from the FITS HDU. Note that the
 * energies are converted to TeV and the effective area is converted to cm2.
 * Conversion is done based on the units provided for the energy and
 * effective area columns. Units that are recognized are 'keV', 'MeV', 'GeV',
 * 'TeV', 'm^2', 'm2', 'cm^2' and 'cm^2' (case independent).
 *
 * @todo Assign appropriate theta angle for PSF. So far we use onaxis.
 *       For appropriate theta angle assignment, we would need this
 *       information in the response header.
 ***************************************************************************/
void GCTAAeffArf::read_arf(const GFitsTable* hdu)
{
    // Clear arrays
    m_logE.clear();
    m_aeff.clear();

    // Get pointers to table columns
    const GFitsTableCol* energy_lo = (*hdu)["ENERG_LO"];
    const GFitsTableCol* energy_hi = (*hdu)["ENERG_HI"];
    const GFitsTableCol* specresp  = (*hdu)["SPECRESP"];

    // Determine unit conversion factors (default: TeV and cm^2)
    std::string u_energy_lo = gammalib::tolower(gammalib::strip_whitespace(energy_lo->unit()));
    std::string u_energy_hi = gammalib::tolower(gammalib::strip_whitespace(energy_hi->unit()));
    std::string u_specresp  = gammalib::tolower(gammalib::strip_whitespace(specresp->unit()));
    double c_energy_lo = 1.0;
    double c_energy_hi = 1.0;
    double c_specresp  = 1.0;
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
    if (u_specresp == "m^2" || u_specresp == "m2") {
        c_specresp = 10000.0;
    }

    // Extract number of energy bins
    int num = energy_lo->length();

    // Set nodes
    for (int i = 0; i < num; ++i) {
    
        // Compute log10 mean energy in TeV
        double e_min = energy_lo->real(i) * c_energy_lo;
        double e_max = energy_hi->real(i) * c_energy_hi;
        double logE  = 0.5 * (log10(e_min) + log10(e_max));

        // Initialise scale factor
        double scale = m_scale;

        // Compute effective area in cm2
        double aeff = specresp->real(i) * c_specresp * scale;
        
        // Store log10 mean energy and effective area value
        m_logE.append(logE);
        m_aeff.push_back(aeff);

    } // endfor: looped over nodes
    
    // Disable offset angle dependence
    m_sigma = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove thetacut
 *
 * @param[in] rsp CTA response.
 *
 * Removes thetacut from Aeff values read from a FITS file. Note that this
 * method should only be called once directly after loading all response
 * components.
 ***************************************************************************/
void GCTAAeffArf::remove_thetacut(const GCTAResponse& rsp)
{
    // Continue only if thetacut value has been set
    if (m_thetacut > 0.0) {

        // Get maximum integration radius
        double rmax = m_thetacut * gammalib::deg2rad;

        // Loop over Aeff array
        for (int i = 0; i < m_aeff.size(); ++i) {
    
            // Setup integration kernel for on-axis PSF
            cta_npsf_kern_rad_azsym integrand(rsp,
                                              rmax,
                                              0.0,
                                              m_logE[i],
                                              0.0,
                                              0.0,
                                              0.0,
                                              0.0);

            // Setup integration
            GIntegral integral(&integrand);
            integral.eps(rsp.eps());

            // Perform integration
            double fraction = integral.romb(0.0, rmax);

            // Set scale factor
            double scale = 1.0;
            if (fraction > 0.0) {
                scale /= fraction;
                #if defined(G_DEBUG_APPLY_THETACUT)
                std::cout << "GCTAAeffArf::apply_thetacut:";
                std::cout << " logE=" << m_logE[i];
                std::cout << " scale=" << scale;
                std::cout << " fraction=" << fraction;
                std::cout << std::endl;
                #endif
            }
            else {
                std::cout << "WARNING: GCTAAeffArf::apply_thetacut:";
                std::cout << " Non-positive integral occured in";
                std::cout << " PSF integration.";
                std::cout << std::endl;
            }

            // Apply scaling factor
            if (scale != 1.0) {
                m_aeff[i] *= scale;
            }

        } // endfor: looped over Aeff array

    } // endif: thetacut value was set

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print effective area information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing effective area information.
 ***************************************************************************/
std::string GCTAAeffArf::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        double emin = std::pow(10.0, m_logE[0]);
        double emax = std::pow(10.0, m_logE[size()-1]);

        // Append header
        result.append("=== GCTAAeffArf ===");

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
void GCTAAeffArf::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_logE.clear();
    m_aeff.clear();
    m_sigma    = 0.0;
    m_thetacut = 0.0;
    m_scale    = 1.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] aeff Effective area.
 ***************************************************************************/
void GCTAAeffArf::copy_members(const GCTAAeffArf& aeff)
{
    // Copy members
    m_filename = aeff.m_filename;
    m_logE     = aeff.m_logE;
    m_aeff     = aeff.m_aeff;
    m_sigma    = aeff.m_sigma;
    m_thetacut = aeff.m_thetacut;
    m_scale    = aeff.m_scale;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAAeffArf::free_members(void)
{
    // Return
    return;
}
