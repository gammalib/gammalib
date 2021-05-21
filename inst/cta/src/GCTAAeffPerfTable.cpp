/***************************************************************************
 *    GCTAAeffPerfTable.hpp - CTA performance table effective area class   *
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
 * @file GCTAAeffPerfTable.hpp
 * @brief CTA performance table effective area class definition
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>             // std::fopen, std::fgets, and std::fclose
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GFitsTable.hpp"
#include "GCTAAeffPerfTable.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                        "GCTAAeffPerfTable::load(std::string&)"

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
 ***************************************************************************/
GCTAAeffPerfTable::GCTAAeffPerfTable(void) : GCTAAeff()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Performance table file name.
 *
 * Construct instance by loading the effective area information from an
 * ASCII performance table.
 ***************************************************************************/
GCTAAeffPerfTable::GCTAAeffPerfTable(const GFilename& filename) : GCTAAeff()
{
    // Initialise class members
    init_members();

    // Load effective area from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] aeff Effective area.
 ***************************************************************************/
GCTAAeffPerfTable::GCTAAeffPerfTable(const GCTAAeffPerfTable& aeff) : GCTAAeff(aeff)
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
GCTAAeffPerfTable::~GCTAAeffPerfTable(void)
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
GCTAAeffPerfTable& GCTAAeffPerfTable::operator=(const GCTAAeffPerfTable& aeff)
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
 *
 * Returns the effective area in units of cm2 for a given energy and
 * offset angle. The effective area is linearily interpolated in
 * log10(energy). The method assures that the effective area value never
 * becomes negative.
 *
 * Outside the energy range that is covered by the performance table the
 * effective area will be set to zero.
 ***************************************************************************/
double GCTAAeffPerfTable::operator()(const double& logE, 
                                     const double& theta, 
                                     const double& phi,
                                     const double& zenith,
                                     const double& azimuth,
                                     const bool&   etrue) const
{
    // Initialise effective area
    double aeff = 0.0;

    // Continue only if logE is in validity range
    if ((logE >= m_logE_min)  && (logE <= m_logE_max)) {

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
            double scale  = std::exp(-0.5 * arg * arg);
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
void GCTAAeffPerfTable::clear(void)
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
GCTAAeffPerfTable* GCTAAeffPerfTable::clone(void) const
{
    return new GCTAAeffPerfTable(*this);
}


/***********************************************************************//**
 * @brief Load effective area from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GException::file_error
 *            File could not be opened for read access.
 *
 * This method loads the effective area information from an ASCII
 * performance table.
 ***************************************************************************/
void GCTAAeffPerfTable::load(const GFilename& filename)
{
    // Clear arrays
    m_logE.clear();
    m_aeff.clear();

    // Allocate line buffer
    const int n = 1000;
    char  line[n];

    // Open performance table readonly
    FILE* fptr = std::fopen(filename.url().c_str(), "r");
    if (fptr == NULL) {
        std::string msg = "Effective area file \""+filename.url()+"\" not "
                          "found or readable. Please specify a valid and "
                          "readable effective area file.";
        throw GException::file_error(G_LOAD, msg);
    }

    // Read lines
    while (std::fgets(line, n, fptr) != NULL) {

        // Split line in elements. Strip empty elements from vector.
        std::vector<std::string> elements = gammalib::split(line, " ");
        for (int i = elements.size()-1; i >= 0; i--) {
            if (gammalib::strip_whitespace(elements[i]).length() == 0) {
                elements.erase(elements.begin()+i);
            }
        }

        // Skip header
        if (elements[0].find("log(E)") != std::string::npos) {
            continue;
        }

        // Break loop if end of data table has been reached
        if (elements[0].find("----------") != std::string::npos) {
            break;
        }

        // Push elements in node array and vector
        m_logE.append(gammalib::todouble(elements[0]));
        m_aeff.push_back(gammalib::todouble(elements[1])*10000.0);

    } // endwhile: looped over lines

    // Close file
    std::fclose(fptr);

    // Store filename
    m_filename = filename;

    // Set the energy boundaries
    set_boundaries();

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
double GCTAAeffPerfTable::max(const double& logE,
                              const double& zenith,
                              const double& azimuth,
                              const bool& etrue) const
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
std::string GCTAAeffPerfTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAAeffPerfTable ===");

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
void GCTAAeffPerfTable::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_logE.clear();
    m_aeff.clear();
    m_ebounds.clear();
    m_sigma    = 3.0;
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
void GCTAAeffPerfTable::copy_members(const GCTAAeffPerfTable& aeff)
{
    // Copy members
    m_filename = aeff.m_filename;
    m_logE     = aeff.m_logE;
    m_aeff     = aeff.m_aeff;
    m_ebounds  = aeff.m_ebounds;
    m_sigma    = aeff.m_sigma;
    m_logE_min = aeff.m_logE_min;
    m_logE_max = aeff.m_logE_max;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAAeffPerfTable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set effective area boundaries
 *
 * Sets the data members m_ebounds, m_logE_min and m_logE_max that define
 * the validity range of the effective area.
 ***************************************************************************/
void GCTAAeffPerfTable::set_boundaries(void)
{
    // Clear energy boundaries
    m_ebounds.clear();

    // Set log10 of minimum and maximum energies. Since the energy values are
    // given at the bin centre we subtract half of the distance to the second
    // bin from the minimum energy and we add half of the distance to the before
    // last bin to the maximum energy
    m_logE_min = m_logE[0]               - 0.5*(m_logE[1] - m_logE[0]);
    m_logE_max = m_logE[m_logE.size()-1] + 0.5*(m_logE[m_logE.size()-1] -
                                                m_logE[m_logE.size()-2]);

    // Set energy boundaries
    GEnergy emin;
    GEnergy emax;
    emin.log10TeV(m_logE_min);
    emax.log10TeV(m_logE_max);
    m_ebounds.append(emin, emax);

    // Return
    return;
}
