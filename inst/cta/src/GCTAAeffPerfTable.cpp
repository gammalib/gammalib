/***************************************************************************
 *    GCTAAeffPerfTable.hpp - CTA performance table effective area class   *
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
 * @file GCTAAeffPerfTable.hpp
 * @brief CTA performance table effective area class definition
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>             // std::fopen, std::fgets, and std::fclose
#include "GTools.hpp"
#include "GFitsTable.hpp"
#include "GCTAAeffPerfTable.hpp"
#include "GCTAException.hpp"

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
GCTAAeffPerfTable::GCTAAeffPerfTable(const std::string& filename) : GCTAAeff()
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
GCTAAeffPerfTable& GCTAAeffPerfTable::operator= (const GCTAAeffPerfTable& aeff)
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
double GCTAAeffPerfTable::operator()(const double& logE, 
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
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the effective area information from an ASCII
 * performance table.
 ***************************************************************************/
void GCTAAeffPerfTable::load(const std::string& filename)
{
    // Clear arrays
    m_logE.clear();
    m_aeff.clear();

    // Allocate line buffer
    const int n = 1000;
    char  line[n];

    // Expand environment variables
    std::string fname = expand_env(filename);

    // Open performance table readonly
    FILE* fptr = std::fopen(fname.c_str(), "r");
    if (fptr == NULL) {
        throw GCTAException::file_open_error(G_LOAD, fname);
    }

    // Read lines
    while (std::fgets(line, n, fptr) != NULL) {

        // Split line in elements. Strip empty elements from vector.
        std::vector<std::string> elements = split(line, " ");
        for (int i = elements.size()-1; i >= 0; i--) {
            if (strip_whitespace(elements[i]).length() == 0) {
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
        m_logE.append(todouble(elements[0]));
        m_aeff.push_back(todouble(elements[1])*10000.0);

    } // endwhile: looped over lines

    // Close file
    std::fclose(fptr);

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
std::string GCTAAeffPerfTable::filename(void) const
{
    // Return filename
    return m_filename;
}


/***********************************************************************//**
 * @brief Print effective area information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing effective area information.
 ***************************************************************************/
std::string GCTAAeffPerfTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        double emin = std::pow(10.0, m_logE[0]);
        double emax = std::pow(10.0, m_logE[size()-1]);

        // Append header
        result.append("=== GCTAAeffPerfTable ===");

        // Append information
        result.append("\n"+parformat("Filename")+m_filename);
        result.append("\n"+parformat("Number of energy bins")+str(size()));
        result.append("\n"+parformat("Log10(Energy) range"));
        result.append(str(emin)+" - "+str(emax)+" TeV");

        // Append offset angle dependence
        if (m_sigma == 0) {
            result.append("\n"+parformat("Offset angle dependence")+"none");
        }
        else {
            std::string txt = "Fixed sigma="+str(m_sigma);
            result.append("\n"+parformat("Offset angle dependence")+txt);
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
    m_sigma = 3.0;

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
    m_sigma    = aeff.m_sigma;

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
