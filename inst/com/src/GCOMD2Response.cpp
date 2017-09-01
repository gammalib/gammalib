/***************************************************************************
 *          GCOMD2Response.cpp - COMPTEL D2 module response class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMD2Response.cpp
 * @brief COMPTEL D2 module response class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GMath.hpp"
#include "GIntegral.hpp"
#include "GCOMD2Response.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_UPDATE_RESPONSE_VECTOR_NO_WARNINGS

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_UPDATE_RESPONSE_VECTOR

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates an empty COMPTEL D2 module response.
 ***************************************************************************/
GCOMD2Response::GCOMD2Response(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp COMPTEL D2 module response.
 **************************************************************************/
GCOMD2Response::GCOMD2Response(const GCOMD2Response& rsp)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response constructor
 *
 * @param[in] caldb Calibration database.
 * @param[in] sdbname SDA response name.
 *
 * Create COMPTEL D2 module response by loading an SDB file from a
 * calibration database.
 ***************************************************************************/
GCOMD2Response::GCOMD2Response(const GCaldb& caldb, const std::string& sdbname)
{
    // Initialise members
    init_members();

    // Set calibration database
    this->caldb(caldb);

    // Load D2 module response
    this->load(sdbname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys instance of COMPTEL response object.
 ***************************************************************************/
GCOMD2Response::~GCOMD2Response(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] rsp COMPTEL D2 module response.
 * @return COMPTEL D2 module response.
 ***************************************************************************/
GCOMD2Response& GCOMD2Response::operator=(const GCOMD2Response& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(rsp);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief D2 module response evaluation operator
 *
 * @param[in] etrue True energy (MeV).
 * @param[in] ereco Reconstructed energy (MeV).
 * @return COMPTEL D2 module response.
 *
 * @todo Add Compton components to response
 ***************************************************************************/
double GCOMD2Response::operator()(const double& etrue, const double& ereco) const
{
    // Initialise response and area with zero
    double response = 0.0;

    // Continue only if a response was loaded
    if (!m_energies.is_empty()) {

        // Update response vector
        update_response_vector(etrue);

        // Get response value from response vector
        response = m_rsp_energies.interpolate(ereco, m_rsp_values);
        if (response < 0.0) {
            response = 0.0;
        }

        // Apply threshold. We use a width of 0.1 MeV for the high-energy
        // threshold width to assure a smooth transition to zero.
        if (ereco < m_emin) {
            double arg  = (m_emin-ereco) / m_ewidth;
            response   *= std::exp(-0.5 * arg * arg);
        }
        else if (ereco > m_emax) {
            double arg  = (m_emax-ereco) / 0.1;
            response   *= std::exp(-0.5 * arg * arg);
        }

    } // endif: response was loaded

    // Return response
    return response;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * Clears COMPTEL D2 module response object by resetting all members to an
 * initial state. Any information that was present in the object before will
 * be lost.
 ***************************************************************************/
void GCOMD2Response::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of COMPTEL D2 module response.
 ***************************************************************************/
GCOMD2Response* GCOMD2Response::clone(void) const
{
    return new GCOMD2Response(*this);
}


/***********************************************************************//**
 * @brief Load COMPTEL D2 module response.
 *
 * @param[in] sdbname COMPTEL D2 module response name.
 *
 * Loads the COMPTEL D2 module response with specified name @p sdbname. The
 * method first searchs for an appropriate response in the calibration
 * database. If no appropriate response is found, the method takes the
 * database root path and response name to build the full path to the
 * response file, and tries to load the response from these paths.
 ***************************************************************************/
void GCOMD2Response::load(const std::string& sdbname)
{
    // Clear instance but conserve calibration database
    GCaldb caldb = m_caldb;
    clear();
    m_caldb = caldb;

    // First attempt reading the response using the GCaldb interface
    GFilename filename = m_caldb.filename("","","SDB","","",sdbname);

    // If filename is empty then build filename from CALDB root path and
    // response name
    if (filename.is_empty()) {
        filename = gammalib::filepath(m_caldb.rootdir(), sdbname);
        if (!filename.exists()) {
            GFilename testname = filename + ".fits";
            if (testname.exists()) {
                filename = testname;
            }
        }
    }

    // Open FITS file
    GFits fits(filename);

    // Get SDB table
    const GFitsTable& sdb = *fits.table(1);

    // Read SDB
    read(sdb);

    // Close SDB FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL D2 module response.
 *
 * @param[in] table FITS table.
 *
 * Read the COMPTEL D2 module response from a SDB FITS table.
 ***************************************************************************/
void GCOMD2Response::read(const GFitsTable& table)
{
    // Initialise COMPTEL D2 module response vectors
    m_energies.clear();
    m_positions.clear();
    m_sigmas.clear();
    m_amplitudes.clear();
    m_escapes1.clear();
    m_escapes2.clear();
    m_tails.clear();
    m_backgrounds.clear();
    m_emins.clear();
    m_ewidths.clear();
    m_emaxs.clear();

    // Extract number of entries in table
    int num = table.nrows();

    // If there are entries then read them
    if (num > 0) {

        // Get column pointers
        const GFitsTableCol* ptr_energy     = table["ENERGY"];
        const GFitsTableCol* ptr_position   = table["POSITION"];
        const GFitsTableCol* ptr_sigma      = table["WIDTH"];
        const GFitsTableCol* ptr_amplitude  = table["AMPLITUDE"];
        const GFitsTableCol* ptr_escape1    = table["ESCAPE1"];
        const GFitsTableCol* ptr_escape2    = table["ESCAPE2"];
        const GFitsTableCol* ptr_tail       = table["TAIL"];
        const GFitsTableCol* ptr_background = table["BACKGROUND"];
        const GFitsTableCol* ptr_emin       = table["EMIN"];
        const GFitsTableCol* ptr_ewidth     = table["EWIDTH"];
        const GFitsTableCol* ptr_emax       = table["EMAX"];

        // Copy data from table into vectors
        for (int i = 0; i < num; ++i) {
            m_energies.append(ptr_energy->real(i));
            m_positions.push_back(ptr_position->real(i));
            m_sigmas.push_back(ptr_sigma->real(i));
            m_amplitudes.push_back(ptr_amplitude->real(i));
            m_escapes1.push_back(ptr_escape1->real(i));
            m_escapes2.push_back(ptr_escape2->real(i));
            m_tails.push_back(ptr_tail->real(i));
            m_backgrounds.push_back(ptr_background->real(i));
            m_emins.push_back(ptr_emin->real(i));
            m_ewidths.push_back(ptr_ewidth->real(i));
            m_emaxs.push_back(ptr_emax->real(i));
        }

    } // endif: there were entries

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL D2 module response information
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL D2 module response information.
 ***************************************************************************/
std::string GCOMD2Response::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMD2Response ===");

        // Append D2 module response information
        result.append("\n"+gammalib::parformat("Energy range"));
        if (m_energies.size() > 0) {
            result.append(gammalib::str(m_energies[0])+ " - ");
            result.append(gammalib::str(m_energies[m_energies.size()-1])+ " MeV");
        }
        else {
            result.append("not defined");
        }
        result.append("\n"+gammalib::parformat("Entries"));
        result.append(gammalib::str(m_energies.size()));

        // Append calibration database
        result.append("\n"+m_caldb.print(chatter));

        // Append information

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCOMD2Response::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_energies.clear();
    m_positions.clear();
    m_sigmas.clear();
    m_amplitudes.clear();
    m_escapes1.clear();
    m_escapes2.clear();
    m_tails.clear();
    m_backgrounds.clear();
    m_emins.clear();
    m_ewidths.clear();
    m_emaxs.clear();

    // Initialise pre-computation cache
    m_energy       = 0.0;
    m_position     = 0.0;
    m_sigma        = 0.0;
    m_amplitude    = 0.0;
    m_escape1      = 0.0;
    m_escape2      = 0.0;
    m_tail         = 0.0;
    m_background   = 0.0;
    m_emin         = 0.0;
    m_ewidth       = 0.0;
    m_emax         = 0.0;
    m_pos_escape1  = 0.0;
    m_pos_escape2  = 0.0;
    m_wgt_photo    = 0.0;
    m_wgt_escape1  = 0.0;
    m_wgt_escape2  = 0.0;
    m_compton_edge = 0.0;

    // Pre-computed response vector
    m_rsp_etrue = 0.0;
    m_rsp_energies.clear();
    m_rsp_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp COMPTEL response.
 ***************************************************************************/
void GCOMD2Response::copy_members(const GCOMD2Response& rsp)
{
    // Copy attributes
    m_caldb       = rsp.m_caldb;
    m_energies    = rsp.m_energies;
    m_positions   = rsp.m_positions;
    m_sigmas      = rsp.m_sigmas;
    m_amplitudes  = rsp.m_amplitudes;
    m_escapes1    = rsp.m_escapes1;
    m_escapes2    = rsp.m_escapes2;
    m_tails       = rsp.m_tails;
    m_backgrounds = rsp.m_backgrounds;
    m_emins       = rsp.m_emins;
    m_ewidths     = rsp.m_ewidths;
    m_emaxs       = rsp.m_emaxs;

    // Copy pre-computation cache
    m_energy       = rsp.m_energy;
    m_position     = rsp.m_position;
    m_sigma        = rsp.m_sigma;
    m_amplitude    = rsp.m_amplitude;
    m_escape1      = rsp.m_escape1;
    m_escape2      = rsp.m_escape2;
    m_tail         = rsp.m_tail;
    m_background   = rsp.m_background;
    m_emin         = rsp.m_emin;
    m_ewidth       = rsp.m_ewidth;
    m_emax         = rsp.m_emax;
    m_pos_escape1  = rsp.m_pos_escape1;
    m_pos_escape2  = rsp.m_pos_escape2;
    m_wgt_photo    = rsp.m_wgt_photo;
    m_wgt_escape1  = rsp.m_wgt_escape1;
    m_wgt_escape2  = rsp.m_wgt_escape2;
    m_compton_edge = rsp.m_compton_edge;

    // Copy pre-computed response vector
    m_rsp_etrue    = rsp.m_rsp_etrue;
    m_rsp_energies = rsp.m_rsp_energies;
    m_rsp_values   = rsp.m_rsp_values;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMD2Response::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update computation cache
 *
 * @param[in] etrue True energy (MeV).
 *
 * The method assumes that there is a valid D2 module response.
 ***************************************************************************/
void GCOMD2Response::update_cache(const double& etrue) const
{
    // Update only if the true energy has changed
    if (etrue != m_energy) {

        // Set true energy
        m_energy = etrue;

        // If true energy is below lowest energy or above largest energy
        // then set response to zero
        if ((etrue < m_energies[0]) ||
            (etrue > m_energies[m_energies.size()-1])) {
            m_position     = 0.0;
            m_sigma        = 0.0;
            m_amplitude    = 0.0;
            m_escape1      = 0.0;
            m_escape2      = 0.0;
            m_tail         = 0.0;
            m_background   = 0.0;
            m_emin         = 0.0;
            m_ewidth       = 0.0;
            m_emax         = 0.0;
            m_pos_escape1  = 0.0;
            m_pos_escape2  = 0.0;
            m_wgt_escape1  = 0.0;
            m_wgt_escape2  = 0.0;
            m_compton_edge = 0.0;
        }

        // ... otherwise interpolate response parameters
        else {

            // Interpolate response parameters
            m_position   = m_energies.interpolate(etrue, m_positions);
            m_sigma      = m_energies.interpolate(etrue, m_sigmas);
            m_amplitude  = m_energies.interpolate(etrue, m_amplitudes);
            m_escape1    = m_energies.interpolate(etrue, m_escapes1);
            m_escape2    = m_energies.interpolate(etrue, m_escapes2);
            m_tail       = m_energies.interpolate(etrue, m_tails);
            m_background = m_energies.interpolate(etrue, m_backgrounds);
            m_emin       = m_energies.interpolate(etrue, m_emins);
            m_ewidth     = m_energies.interpolate(etrue, m_ewidths);
            m_emax       = m_energies.interpolate(etrue, m_emaxs);

            // Derive escape peak positions
            m_pos_escape1 = m_position - gammalib::mec2;
            m_pos_escape2 = m_position - 2.0 * gammalib::mec2;

            // Derive inverse standard deviations of photo and escape peaks
            m_wgt_photo = m_sigma;
            if (m_wgt_photo > 0.0) {
                m_wgt_photo = 1.0 / m_wgt_photo;
            }
            m_wgt_escape1 = m_energies.interpolate(m_pos_escape1, m_sigmas);
            if (m_wgt_escape1 > 0.0) {
                m_wgt_escape1 = 1.0 / m_wgt_escape1;
            }
            m_wgt_escape2 = m_energies.interpolate(m_pos_escape2, m_sigmas);
            if (m_wgt_escape2 > 0.0) {
                m_wgt_escape2 = 1.0 / m_wgt_escape2;
            }

            // Derive Compton edge
            m_compton_edge = m_position / (1.0 + 0.5 * gammalib::mec2 / m_position);

        } // endif: true energy is in valid range

    } // endif: true energy has changed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update response vector
 *
 * @param[in] etrue True energy (MeV).
 *
 * Updates the vector that stores the normalised D2 module response as
 * function of energy. The vector is computed using
 *
 * \f[
 *    R_{\rm D2}(E|E_0) = n \times \left[
 *    B_1 \exp \left( -\frac{1}{2} \frac{(E_0-E)^2}{\sigma^2(E_0)} \right) +
 *    B_2 \exp \left( -\frac{1}{2} \frac{(E_0-m_e c^2-E)^2}{\sigma^2(E_0-m_e c^2)} \right) +
 *    B_3 \exp \left( -\frac{1}{2} \frac{(E_0-2m_e c^2-E)^2}{\sigma^2(E_0-2m_e c^2)} \right) +
 *    B_4 \int_{E'} KN_{\rm mod}(E'|E,E_0) dE' +
 *    B_5 \int_{E'} B_{\rm c}(E'|E,E_0) dE'
 *    \right]
 * \f]
 *
 * where
 * \f$B1\f$ is the amplitude of the photo peak,
 * \f$B2\f$ is the amplitude of the first escape peak,
 * \f$B3\f$ is the amplitude of the second escape peak,
 * \f$B4\f$ is the amplitude of the Compton tail,
 * \f$B5\f$ is the amplitude of the Compton background,
 * \f$E_0\f$ is the position of the photo peak, and
 * \f$\sigma(E)\f$ is the energy dependent width of the photo peak. The
 * constant \f$n\f$ is chosen so that
 *
 * \f[
 *    \int_{E} R_{\rm D2}(E|E_0) dE = 1
 * \f]
 *
 * The method assumes that there is a valid D2 module response.
 ***************************************************************************/
void GCOMD2Response::update_response_vector(const double& etrue) const
{
    // Update only if the true energy has changed
    if (etrue != m_rsp_etrue) {

        // Set true energy
        m_rsp_etrue = etrue;

        // Update cache to determine the spectral parameters at the true
        // energy
        update_cache(etrue);
        
        // Continue only if position is positive
        if (m_position > 0.0) {
        
            // Clear response vectors (this is very very important, otherwise
            // some old elements reside and we just append to them, and things
            // get out of order!!!)
            m_rsp_energies.clear();
            m_rsp_values.clear();

            // Initialise response vector
            double ebin = 0.001 * m_position;
            int    nbin = int(1.35  * m_position / ebin);
            if (nbin < 1) {
                nbin = 1;
            }
            m_rsp_energies.reserve(nbin);
            m_rsp_values.reserve(nbin);

            // Initialise response vector normalisation
            double norm = 0.0;

            // Fill response vector
            for (int i = 0; i < nbin; ++i) {

                // Compute bin energy
                double ereco = (i + 0.5) * ebin;

                // Store bin energy
                m_rsp_energies.append(ereco);

                // Initialise response value
                double value = 0.0;

                // Add photo peak if amplitude is positive and reconstructed
                // energy is within 5 sigma of the Gaussian position
                if (m_amplitude > 0.0) {
                    double arg  = (m_position-ereco) * m_wgt_photo;
                    if (std::abs(arg) < 5.0) {
                        value += m_amplitude * std::exp(-0.5 * arg * arg);
                    }
                }

                // Add first escape peak if amplitude is positive and
                // reconstructed energy is within 5 sigma of the Gaussian
                // position
                if (m_escape1 > 0.0) {
                    double arg  = (m_pos_escape1-ereco) * m_wgt_escape1;
                    if (std::abs(arg) < 5.0) {
                        value += m_escape1 * std::exp(-0.5 * arg * arg);
                    }
                }

                // Add second escape peak if amplitude is positive and
                // reconstructed energy is within 5 sigma of the Gaussian
                // position
                if (m_escape2 > 0.0) {
                    double arg  = (m_pos_escape2-ereco) * m_wgt_escape2;
                    if (std::abs(arg) < 5.0) {
                        value += m_escape2 * std::exp(-0.5 * arg * arg);
                    }
                }

                // Compute sigma for continua. Since this may not work
                double sigma = 0.0;
                if ((m_tail > 0.0) || (m_background > 0.0)) {
                    double e = (ereco > m_emin) ? ereco : m_emin;
                    sigma    = m_energies.interpolate(e, m_sigmas);
                }

                // Add Compton tail if amplitude is positive
                if (m_tail > 0.0) {

                    // Setup integration kernel for Compton tail
                    kn_gauss_kernel integrand(ereco, m_position, m_compton_edge, sigma);

                    // Setup integral
                    GIntegral integral(&integrand);

                    // Set precision
                    integral.eps(1.0e-5);

                    // No warnings
                    #if defined(G_UPDATE_RESPONSE_VECTOR_NO_WARNINGS)
                    integral.silent(true);
                    #endif

                    // Perform integration over Gaussian width
                    double emin = ereco - 3.0 * sigma;
                    double emax = ereco + 3.0 * sigma;
                    value += m_tail * integral.romberg(emin, emax);

                } // endif: added Compton tail

                // Add Compton background if amplitude is positive
                if (m_background > 0.0) {

                    // Setup integration kernel for Compton background
                    bkg_gauss_kernel integrand(ereco, m_position, sigma);

                    // Setup integral
                    GIntegral integral(&integrand);

                    // Set precision
                    integral.eps(1.0e-5);

                    // No warnings
                    #if defined(G_UPDATE_RESPONSE_VECTOR_NO_WARNINGS)
                    integral.silent(true);
                    #endif

                    // Perform integration over Gaussian width
                    double emin = ereco - 3.0 * sigma;
                    double emax = ereco + 3.0 * sigma;
                    value += m_background * integral.romberg(emin, emax);

                } // endif: added Compton background

                // Store response value
                m_rsp_values.push_back(value);

                // Increment normalisation
                norm += value;

            } // endfor: looped over response vector bins

            // Normalise response vector
            if (norm > 0.0) {
                norm = 1.0 / (norm * ebin);
                for (int i = 0; i < nbin; ++i) {
                    m_rsp_values[i] *= norm;
                }
            }

            // Debug
            #if defined(G_DEBUG_UPDATE_RESPONSE_VECTOR)
            std::cout << "etrue=" << etrue;
            std::cout << " ebin=" << ebin;
            std::cout << " nbin=" << nbin;
            std::cout << " norm=" << norm << std::endl;
            #endif

        } // endif: photo-peak position was positive

    } // endif: true energy has changed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes modified Klein-Nishina cross section multiplied with
 *        Gaussian
 *
 * @param[in] e Energy (MeV).
 * @return Modified Klein-Nishina cross section multiplied with Gaussian
 *
 * Computes the modified Klein-Nishina cross section multiplied with a
 * Gaussian using
 *
 * \f[
 *    KN_{\rm mod}(E|E',E_0) = \sigma_{\rm KN}(E|E_0) \times (1-C(E|E_0))
 *                             \times \exp \left( -\frac{1}{2}
 *                                  \frac{(E-E')^2}{\sigma^2(E')}
 *                                  \right)
 * \f]
 *
 * where
 *
 * \f[
 *    \sigma_{\rm KN}(E|E_0) = \theta(E_c-E)
 *                             \left[
 *                             \left( \frac{E/E_0}{1-E/E_0}
 *                                    \frac{m_e c^2}{E_0} \right)^2 -
 *                             \frac{E}{E_0} + \frac{1}{1-E/E_0} \right]
 * \f]
 *
 * is the Klein-Nishina cross section, and
 *
 * \f[
 *    C(E|E_0) = \theta(E-12.14) \times (1 - e^{-\mu(E|E_0) \, l(E_0)})
 * \f]
 *
 * is a modification, where \f$\theta(x)\f$ is a step function that is 1 for
 * positive \f$x\f$ and 0 for non-positive \f$x\f$,
 *
 * \f[
 *    \mu(E|E_0) = 0.72 \, e^{-1.28 (E_0 - E)^{0.35}} +
 *                 0.01 \, (E_0 - E) +
 *                 0.014 \, (E_0 - E)^{-2.5}
 * \f]
 *
 * is the total linear attenuation coefficient in NaI for all processes
 * and
 *
 * \f[
 *    l(E_0) = 2.9 \log( E_0 - 11.14)
 * \f]
 *
 * is the pathlength in the D2 module (energies are in MeV). \f$E_0\f$ is
 * the position of the photo peak,
 *
 * \f[
 *    E_c = \frac{E_0}{1 + \frac{m_e c^2}{2 E_0}}
 * \f]
 *
 * is the Compton edge, \f$E'\f$ is the reconstructed energy,
 * \f$\sigma(E')\f$ is the standard deviation at the reconstructed energy,
 * and \f$E\f$ is the energy over which the convolution is performed.
 ***************************************************************************/
double GCOMD2Response::kn_gauss_kernel::eval(const double& e)
{
    // Initialise result
    double value = 0.0;

    // Continue only if the energy is below the Compton edge
    if (e < m_ec) {

        // Compute terms
        double a    = e/m_e0;
        double term = a/(1.0-a) * gammalib::mec2/m_e0 - 1.0;

        // Compute Klein-Nishina cross section
        value = term * term - a + 1.0/(1.0-a);

        // If the incident energy is above 12.14 MeV then multiply-in the
        // high-energy correction that was introduced by Rob van Dijk
        if (m_e0 > 12.14) {
        
            // Compute linear attenuation coefficient
            double d  = m_e0 - e;
            double mu = 0.72 * std::exp(-1.28 * std::pow(d,0.35)) +
                        0.01 * d + 0.014 * std::pow(d,-2.5);

            // Compute path length through module
            double l = 2.9 * std::log(m_e0 - 11.14);

            // Multiply-in correction factor
            value *= std::exp(-mu * l);

        } // endif: Multiplied-in high-energy correction

        // Compute Gaussian
        double arg   = (e - m_ereco) * m_wgt;
        double gauss = std::exp(-0.5*arg*arg);

        // Multiply-in Gaussian
        value *= gauss;

    } // endif: energy was below the Compton edge

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Computes Compton background multiplied with Gaussian
 *
 * @param[in] e Energy (MeV).
 * @return Compton background multiplied with Gaussian
 *
 * Computes the Compton background multiplied with a Gaussian using
 *
 * \f[
 *    B_{\rm c}(E|E',E_0) = \theta(E_0-E)
 *                          \exp \left( -\frac{1}{2}
 *                               \frac{(E-E')^2}{\sigma^2(E')}
 *                               \right)
 * \f]
 *
 * where \f$\theta(x)\f$ is a step function that is 1 for positive \f$x\f$
 * and 0 for non-positive \f$x\f$, \f$E_0\f$ is the position of the
 * photo peak, \f$E'\f$ is the reconstructed energy, \f$\sigma(E')\f$ is
 * the standard deviation at the reconstructed energy, and \f$E\f$
 * is the energy over which the convolution is performed.
 ***************************************************************************/
double GCOMD2Response::bkg_gauss_kernel::eval(const double& e)
{
    // Initialise background
    double value = (e < m_e0) ? 1.0 : 0.0;

    // Compute Gaussian
    double arg   = (e - m_ereco) * m_wgt;
    double gauss = std::exp(-0.5*arg*arg);

    // Multiply-in Gaussian
    value *= gauss;

    // Return value
    return value;
}
