/***************************************************************************
 *          GCOMD1Response.cpp - COMPTEL D1 module response class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMD1Response.cpp
 * @brief COMPTEL D1 module response class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GMath.hpp"
#include "GCOMD1Response.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_RENORMALIZE                   //!< Renormalise Gaussian amplitude

/* __ Debug definitions __________________________________________________ */
#define G_COMPASS_THRESHOLD        //!< Use threshold as defined in COMPASS

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates an empty COMPTEL D1 module response.
 ***************************************************************************/
GCOMD1Response::GCOMD1Response(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp COMPTEL D1 module response.
 **************************************************************************/
GCOMD1Response::GCOMD1Response(const GCOMD1Response& rsp)
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
 * @param[in] sdaname SDA response name.
 *
 * Create COMPTEL D1 module response by loading an SDA file from a
 * calibration database.
 ***************************************************************************/
GCOMD1Response::GCOMD1Response(const GCaldb& caldb, const std::string& sdaname)
{
    // Initialise members
    init_members();

    // Set calibration database
    this->caldb(caldb);

    // Load D1 module response
    this->load(sdaname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys instance of COMPTEL response object.
 ***************************************************************************/
GCOMD1Response::~GCOMD1Response(void)
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
 * @param[in] rsp COMPTEL D1 module response.
 * @return COMPTEL D1 module response.
 ***************************************************************************/
GCOMD1Response& GCOMD1Response::operator=(const GCOMD1Response& rsp)
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
 * @brief D1 module response evaluation operator
 *
 * @param[in] etrue True energy (MeV).
 * @param[in] ereco Reconstructed energy (MeV).
 * @return COMPTEL D1 module response.
 *
 * Computes
 *
 * \f[
 *    R_{\rm D1}(E|E_0) = n \times \left[
 *    A_1 \exp \left( -\frac{1}{2} \frac{(E_0-E)^2}{\sigma^2(E_0)} \right)
 *    \right]
 * \f]
 *
 * where
 * \f$B1\f$ is the amplitude of the photo peak,
 * \f$\sigma(E_0)\f$ is width of the photo peak, and
 * \f$E_0\f$ is the position of the photo peak.
 * The constant \f$n\f$ is chosen so that
 *
 * \f[
 *    \int_{E} R_{\rm D1}(E|E_0) dE = 1
 * \f]
 *
 * The code implementation is based on the COMPASS RESD1 function
 * SPDERC01.RESD1.F (release ?, date ?).
 ***************************************************************************/
double GCOMD1Response::operator()(const double& etrue, const double& ereco) const
{
    // Initialise response with zero
    double response = 0.0;

    // Continue only if a response was loaded
    if (!m_energies.is_empty()) {

        // Update response evaluation cache
        update_cache(etrue);

        // Continue only if amplitude is positive
        if (m_amplitude > 0.0) {

            // Compute D1 module response. Here is where the real magic
            // happens. Only consider the response within 5 sigma.
            double arg = (m_position-ereco) / m_sigma;
            if (std::abs(arg) < 5.0) {
                response = m_amplitude * std::exp(-0.5 * arg * arg);
            }

        }

        #if defined(G_COMPASS_THRESHOLD)
        // Apply threshold according to COMPASS function resd1
        double thresl = m_emin - 0.5 * m_ewidth;
        if (ereco <= thresl) {
            response = 0.0;
        }
        else {
            if (ereco <= thresl + m_ewidth) {
                response *= (ereco - thresl) / m_ewidth;
            }
            else {
                if (ereco > m_emax) {
                    response = 0.0;
                }
            }
        }

        #else
        // Apply smooth threshold. We use a width of 0.1 MeV for the
        // high-energy threshold width to assure a smooth transition to zero.
        if (ereco < m_emin) {
            double arg  = (m_emin-ereco) / m_ewidth;
            response   *= std::exp(-0.5 * arg * arg);
        }
        else if (ereco > m_emax) {
            double arg  = (m_emax-ereco) / 0.1;
            response   *= std::exp(-0.5 * arg * arg);
        }
        #endif

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
 * Clears COMPTEL D1 module response object by resetting all members to an
 * initial state. Any information that was present in the object before will
 * be lost.
 ***************************************************************************/
void GCOMD1Response::clear(void)
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
 * @return Pointer to deep copy of COMPTEL D1 module response.
 ***************************************************************************/
GCOMD1Response* GCOMD1Response::clone(void) const
{
    return new GCOMD1Response(*this);
}


/***********************************************************************//**
 * @brief Load COMPTEL D1 module response.
 *
 * @param[in] sdaname COMPTEL D1 module response name.
 *
 * Loads the COMPTEL D1 module response with specified name @p sdaname. The
 * method first searchs for an appropriate response in the calibration
 * database. If no appropriate response is found, the method takes the
 * database root path and response name to build the full path to the
 * response file, and tries to load the response from these paths.
 ***************************************************************************/
void GCOMD1Response::load(const std::string& sdaname)
{
    // Clear instance but conserve calibration database
    GCaldb caldb = m_caldb;
    clear();
    m_caldb = caldb;

    // First attempt reading the response using the GCaldb interface
    GFilename filename = m_caldb.filename("","","SDA","","",sdaname);

    // If filename is empty then build filename from CALDB root path and
    // response name
    if (filename.is_empty()) {
        filename = gammalib::filepath(m_caldb.rootdir(), sdaname);
        if (!filename.exists()) {
            GFilename testname = filename + ".fits";
            if (testname.exists()) {
                filename = testname;
            }
        }
    }

    // Open FITS file
    GFits fits(filename);

    // Get SDA table
    const GFitsTable& sda = *fits.table(1);

    // Read SDA
    read(sda);

    // Close SDA FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL D1 module response.
 *
 * @param[in] table FITS table.
 *
 * Read the COMPTEL D1 module response from a SDA FITS table.
 ***************************************************************************/
void GCOMD1Response::read(const GFitsTable& table)
{
    // Initialise COMPTEL D1 module response vectors
    m_energies.clear();
    m_positions.clear();
    m_sigmas.clear();
    m_amplitudes.clear();
    m_emins.clear();
    m_ewidths.clear();
    m_emaxs.clear();

    // Extract number of entries in table
    int num = table.nrows();

    // If there are entries then read them
    if (num > 0) {

        // Get column pointers
        const GFitsTableCol* ptr_energy    = table["ENERGY"];
        const GFitsTableCol* ptr_position  = table["POSITION"];  // PD1(1)
        const GFitsTableCol* ptr_sigma     = table["WIDTH"];     // PD1(2)
        const GFitsTableCol* ptr_amplitude = table["AMPLITUDE"]; // PD1(3)
        const GFitsTableCol* ptr_emin      = table["EMIN"];      // PD1(4)
        const GFitsTableCol* ptr_ewidth    = table["EWIDTH"];    // PD1(5)
        const GFitsTableCol* ptr_emax      = table["EMAX"];      // PD1(6)

        // Copy data from table into vectors
        for (int i = 0; i < num; ++i) {
            m_energies.append(ptr_energy->real(i));
            m_positions.push_back(ptr_position->real(i));
            m_sigmas.push_back(ptr_sigma->real(i));
            m_amplitudes.push_back(ptr_amplitude->real(i));
            m_emins.push_back(ptr_emin->real(i));
            m_ewidths.push_back(ptr_ewidth->real(i));
            m_emaxs.push_back(ptr_emax->real(i));
        }

    } // endif: there were entries

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write COMPTEL D1 module response.
 *
 * @param[in] table FITS table.
 *
 * Write the COMPTEL D1 module response into a SDA FITS table.
 ***************************************************************************/
void GCOMD1Response::write(GFitsBinTable& table)
{
    // Set extension name
    table.extname("SDA");

    // Set number of entries
    int num = m_energies.size();

    // Allocate columns
    GFitsTableFloatCol  col_energy("ENERGY", num);
    GFitsTableDoubleCol col_position("POSITION", num);
    GFitsTableDoubleCol col_width("WIDTH", num);
    GFitsTableDoubleCol col_amplitude("AMPLITUDE", num);
    GFitsTableDoubleCol col_emin("EMIN", num);
    GFitsTableDoubleCol col_ewidth("EWIDTH", num);
    GFitsTableDoubleCol col_emax("EMAX", num);

    // Set column units
    col_energy.unit("MeV");
    col_position.unit("MeV");
    col_width.unit("MeV");
    col_emin.unit("MeV");
    col_ewidth.unit("MeV");
    col_emax.unit("MeV");

    // Copy data from vectors into table
    for (int i = 0; i < num; ++i) {
        col_energy(i)    = m_energies[i];
        col_position(i)  = m_positions[i];
        col_width(i)     = m_sigmas[i];
        col_amplitude(i) = m_amplitudes[i];
        col_emin(i)      = m_emins[i];
        col_ewidth(i)    = m_ewidths[i];
        col_emax(i)      = m_emaxs[i];
    }

    // Append columns to table
    table.append(col_energy);
    table.append(col_position);
    table.append(col_width);
    table.append(col_amplitude);
    table.append(col_emin);
    table.append(col_ewidth);
    table.append(col_emax);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL D1 module response information
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL D1 module response information.
 ***************************************************************************/
std::string GCOMD1Response::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMD1Response ===");

        // Append D1 module response information
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
void GCOMD1Response::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_energies.clear();
    m_positions.clear();
    m_sigmas.clear();
    m_amplitudes.clear();
    m_emins.clear();
    m_ewidths.clear();
    m_emaxs.clear();

    // Initialise pre-computation cache
    m_energy    = 0.0;
    m_position  = 0.0;
    m_sigma     = 0.0;
    m_amplitude = 0.0;
    m_emin      = 0.0;
    m_ewidth    = 0.0;
    m_emax      = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp COMPTEL response.
 ***************************************************************************/
void GCOMD1Response::copy_members(const GCOMD1Response& rsp)
{
    // Copy attributes
    m_caldb      = rsp.m_caldb;
    m_energies   = rsp.m_energies;
    m_positions  = rsp.m_positions;
    m_sigmas     = rsp.m_sigmas;
    m_amplitudes = rsp.m_amplitudes;
    m_emins      = rsp.m_emins;
    m_ewidths    = rsp.m_ewidths;
    m_emaxs      = rsp.m_emaxs;

    // Copy pre-computation cache
    m_energy    = rsp.m_energy;
    m_position  = rsp.m_position;
    m_sigma     = rsp.m_sigma;
    m_amplitude = rsp.m_amplitude;
    m_emin      = rsp.m_emin;
    m_ewidth    = rsp.m_ewidth;
    m_emax      = rsp.m_emax;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMD1Response::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update computation cache
 *
 * @param[in] etrue True energy (MeV).
 *
 * The method assumes that there is a valid D1 module response.
 ***************************************************************************/
void GCOMD1Response::update_cache(const double& etrue) const
{
    // Update only if the true energy has changed
    if (etrue != m_energy) {

        // Set true energy
        m_energy = etrue;

        // If true energy is below lowest energy or above largest energy
        // then set response to zero
        if ((etrue < m_energies[0]) ||
            (etrue > m_energies[m_energies.size()-1])) {
            m_position  = 0.0;
            m_sigma     = 0.0;
            m_amplitude = 0.0;
            m_emin      = 0.0;
            m_ewidth    = 0.0;
            m_emax      = 0.0;
        }

        // ... otherwise interpolate response parameters
        else {

            // Interpolate response parameters
            m_position  = m_energies.interpolate(etrue, m_positions);
            m_sigma     = m_energies.interpolate(etrue, m_sigmas);
            m_amplitude = m_energies.interpolate(etrue, m_amplitudes);
            m_emin      = m_energies.interpolate(etrue, m_emins);
            m_ewidth    = m_energies.interpolate(etrue, m_ewidths);
            m_emax      = m_energies.interpolate(etrue, m_emaxs);

            // Re-compute Gaussian amplitude to assure normalization
            #if defined(G_RENORMALIZE)
            m_amplitude = 1.0 / (m_sigma * gammalib::sqrt_twopi);
            #endif

        }

    } // endif: true energy has changed

    // Return
    return;
}
