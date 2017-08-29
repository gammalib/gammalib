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
#include "GCOMD2Response.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

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
 * @param[in] sdbname SDA response name.
 * @param[in] caldb Calibration database.
 *
 * Create COMPTEL D2 module response response by loading an SDB file from a
 * calibration database.
 ***************************************************************************/
GCOMD2Response::GCOMD2Response(const std::string& sdbname, const GCaldb& caldb)
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
 * @param[in] etrue True energy.
 * @param[in] ereco Reconstructed energy.
 * @return COMPTEL D2 module response.
 ***************************************************************************/
double GCOMD2Response::operator()(const GEnergy& etrue, const GEnergy& ereco) const
{
    // Initialise response with zero
    double response = 0.0;

    // Continue only if a response was loaded
    if (!m_energies.is_empty()) {

        // Update response evaluation cache
        update_cache(etrue);

        // Continue only if amplitude is positive
        if (m_amplitude > 0.0) {

            // Compute D2 module response (here is where the real magic
            // happens)
            //TODO: Implement the correct formula
            double arg = (m_energy.MeV()-ereco.MeV()) / m_sigma;
            response   = m_amplitude * std::exp(-0.5 * arg * arg);

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
    m_energy     = GEnergy(-1.0e30, "MeV");  // To assure initialisation
    m_position   = 0.0;
    m_sigma      = 0.0;
    m_amplitude  = 0.0;
    m_escape1    = 0.0;
    m_escape2    = 0.0;
    m_tail       = 0.0;
    m_background = 0.0;
    m_emin       = 0.0;
    m_ewidth     = 0.0;
    m_emax       = 0.0;

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
    m_energy     = rsp.m_energy;
    m_position   = rsp.m_position;
    m_sigma      = rsp.m_sigma;
    m_amplitude  = rsp.m_amplitude;
    m_escape1    = rsp.m_escape1;
    m_escape2    = rsp.m_escape2;
    m_tail       = rsp.m_tail;
    m_background = rsp.m_background;
    m_emin       = rsp.m_emin;
    m_ewidth     = rsp.m_ewidth;
    m_emax       = rsp.m_emax;

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
 * @param[in] etrue True energy.
 *
 * The method assumes that there is a valid D2 module response.
 ***************************************************************************/
void GCOMD2Response::update_cache(const GEnergy& etrue) const
{
    // Update only if the true energy has changed
    if (etrue != m_energy) {

        // Set true energy
        m_energy = etrue;

        // Get true energy in MeV
        double etrue_MeV = etrue.MeV();

        // If true energy is below lowest energy or above largest energy
        // then set response to zero
        if ((etrue_MeV < m_energies[0]) ||
            (etrue_MeV > m_energies[m_energies.size()-1])) {
            m_position   = 0.0;
            m_sigma      = 0.0;
            m_amplitude  = 0.0;
            m_escape1    = 0.0;
            m_escape2    = 0.0;
            m_tail       = 0.0;
            m_background = 0.0;
            m_emin       = 0.0;
            m_ewidth     = 0.0;
            m_emax       = 0.0;
        }

        // ... otherwise interpolate response parameters
        else {

            // Interpolate response parameters
            m_position   = m_energies.interpolate(etrue_MeV, m_positions);
            m_sigma      = m_energies.interpolate(etrue_MeV, m_sigmas);
            m_amplitude  = m_energies.interpolate(etrue_MeV, m_amplitudes);
            m_escape1    = m_energies.interpolate(etrue_MeV, m_escapes1);
            m_escape2    = m_energies.interpolate(etrue_MeV, m_escapes2);
            m_tail       = m_energies.interpolate(etrue_MeV, m_tails);
            m_background = m_energies.interpolate(etrue_MeV, m_backgrounds);
            m_emin       = m_energies.interpolate(etrue_MeV, m_emins);
            m_ewidth     = m_energies.interpolate(etrue_MeV, m_ewidths);
            m_emax       = m_energies.interpolate(etrue_MeV, m_emaxs);

            //TODO: Assure normalization

        }

    } // endif: true energy has changed

    // Return
    return;
}
