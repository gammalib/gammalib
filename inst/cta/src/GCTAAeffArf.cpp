/***************************************************************************
 *                GCTAAeffArf.cpp - CTA ARF effective area class           *
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
#include "GFilename.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"
#include "GArf.hpp"
#include "GCTAAeffArf.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAResponse_helpers.hpp"

/* __ Method name definitions ____________________________________________ */

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
 * Constructs effective area from an ARF FITS file.
 ***************************************************************************/
GCTAAeffArf::GCTAAeffArf(const GFilename&  filename) : GCTAAeff()
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
GCTAAeffArf& GCTAAeffArf::operator=(const GCTAAeffArf& aeff)
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
 * @param[in] theta Offset angle in camera system (rad).
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 * @param[in] etrue Use true energy (true/false). Not used.
 * @return Effective area in cm2.
 *
 * Returns the effective area in units of cm2 for a given energy and
 * offset angle. The effective area is linearily interpolated in
 * log10(energy). The method assures that the effective area value never
 * becomes negative.
 *
 * Outside the energy range that is covered by the ARF vector the effective
 * area will be set to zero.
 ***************************************************************************/
double GCTAAeffArf::operator()(const double& logE, 
                               const double& theta, 
                               const double& phi,
                               const double& zenith,
                               const double& azimuth,
                               const bool&   etrue) const
{
    // Initialise effective area
    double aeff = 0.0;

    // Continue only if logE is in validity range
    if ((logE  >= m_logE_min)  && (logE  <= m_logE_max)) {

        // Get effective area value in cm2
        aeff = m_logE.interpolate(logE, m_aeff);

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

    } // endif: logE in validity range

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
 * @brief Load effective area from ARF FITS file
 *
 * @param[in] filename ARF FITS file name.
 *
 * Loads the effective area from an ARF FITS file.
 *
 * If no extension name is provided, the effective area will be loaded from
 * the `SPECRESP` extension.
 ***************************************************************************/
void GCTAAeffArf::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get ARF table
    const GFitsTable& table = *fits.table(filename.extname(gammalib::extname_arf));

    // Read ARF from table
    read(table);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA ARF vector
 *
 * @param[in] table FITS table.
 *
 * Reads a CTA ARF vector from the FITS @p table.
 *
 * The energies are converted to TeV and the effective area is converted to
 * cm2. Conversion is done based on the units provided for the energy and
 * effective area columns. Units that are recognized are 'keV', 'MeV', 'GeV',
 * 'TeV', 'm^2', 'm2', 'cm^2' and 'cm^2' (case independent).
 *
 * @todo Assign appropriate theta angle for PSF. So far we use onaxis.
 *       For appropriate theta angle assignment, we would need this
 *       information in the response header.
 ***************************************************************************/
void GCTAAeffArf::read(const GFitsTable& table)
{
    // Clear arrays
    m_logE.clear();
    m_aeff.clear();
    m_ebounds.clear();

    // Get pointers to table columns
    const GFitsTableCol* energy_lo = table["ENERG_LO"];
    const GFitsTableCol* energy_hi = table["ENERG_HI"];
    const GFitsTableCol* specresp  = table["SPECRESP"];

    // Determine unit conversion factors (default: TeV and cm^2)
    std::string u_specresp  = gammalib::tolower(
                              gammalib::strip_whitespace(specresp->unit()));
    double c_specresp  = 1.0;
    if (u_specresp == "m^2" || u_specresp == "m2") {
        c_specresp = 10000.0;
    }

    // Extract number of energy bins
    int num = energy_lo->nrows();

    // Set nodes
    for (int i = 0; i < num; ++i) {

        // Compute energy boundaries
        GEnergy emin(energy_lo->real(i), energy_lo->unit());
        GEnergy emax(energy_hi->real(i), energy_hi->unit());

        // Compute log10 of mean energy in TeV
        double logE = 0.5 * (emin.log10TeV() + emax.log10TeV());

        // Initialise scale factor
        double scale = m_scale;

        // Compute effective area in cm2
        double aeff = specresp->real(i) * c_specresp * scale;
        
        // Store log10 mean energy and effective area value
        m_logE.append(logE);
        m_aeff.push_back(aeff);

    } // endfor: looped over nodes

    // Set energy boundaries
    GEnergy emin(energy_lo->real(0),     energy_lo->unit());
    GEnergy emax(energy_hi->real(num-1), energy_hi->unit());
    m_logE_min = emin.log10TeV();
    m_logE_max = emax.log10TeV();
    m_ebounds.append(emin, emax);
    
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
void GCTAAeffArf::remove_thetacut(const GCTAResponseIrf& rsp)
{
    // Set number of iterations for Romberg integration.
    static const int iter = 6;

    // Continue only if thetacut value has been set
    if (m_thetacut > 0.0) {

        // Get maximum integration radius
        double rmax = m_thetacut * gammalib::deg2rad;

        // Loop over Aeff array
        for (int i = 0; i < m_aeff.size(); ++i) {
    
            // Setup integration kernel for on-axis PSF
            cta_npsf_kern_rad_azsym integrand(&rsp,
                                              rmax,
                                              0.0,
                                              m_logE[i],
                                              0.0,
                                              0.0,
                                              0.0,
                                              0.0);

            // Setup integration
            GIntegral integral(&integrand);

            // Set fixed number of iterations
            integral.fixed_iter(iter);

            // Perform integration
            double fraction = integral.romberg(0.0, rmax);

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
 * @brief Return maximum effective area at a given energy
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] zenith Zenith angle in Earth system (rad). Not used in this method.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used in this method.
 * @param[in] etrue Use true energy (true/false). Not used.
 * @return Maximum effective area (cm2).
 ***************************************************************************/
double GCTAAeffArf::max(const double& logE,
                        const double& zenith,
                        const double& azimuth,
                        const bool&   etrue) const
{
    // Get effective area value in cm2
    double aeff_max = m_logE.interpolate(logE, m_aeff);

    // Make sure that effective area is not negative
    if (aeff_max < 0.0) {
        aeff_max = 0.0;
    }

    // Return result
    return aeff_max;
}


/***********************************************************************//**
 * @brief Print effective area information
 *
 * @param[in] chatter Chattiness.
 * @return String containing effective area information.
 ***************************************************************************/
std::string GCTAAeffArf::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAAeffArf ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(m_ebounds.emin().print() + " - " +
                      m_ebounds.emax().print());

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
    m_ebounds.clear();
    m_sigma    = 0.0;
    m_thetacut = 0.0;
    m_scale    = 1.0;
    m_logE_min = 0.0;
    m_logE_max = 0.0;

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
    m_ebounds  = aeff.m_ebounds;
    m_sigma    = aeff.m_sigma;
    m_thetacut = aeff.m_thetacut;
    m_scale    = aeff.m_scale;
    m_logE_min = aeff.m_logE_min;
    m_logE_max = aeff.m_logE_max;

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
