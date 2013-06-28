/***************************************************************************
 *               GArf.cpp - XSPEC Auxiliary Response File class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GArf.cpp
 * @brief XSPEC Auxiliary Response File class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GArf.hpp"
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                                 "GArf::at(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GArf::GArf(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name.
 ***************************************************************************/
GArf::GArf(const std::string& filename)
{
    // Initialise members
    init_members();

    // Load ARF file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy boundary constructor
 *
 * @param[in] ebds Energy boundaries.
 ***************************************************************************/
GArf::GArf(const GEbounds& ebds)
{
    // Initialise members
    init_members();

    // Set energy boundaries
    m_ebounds = ebds;

    // Initialize spectral response
    m_specresp.assign(ebds.size(), 0.0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] arf Auxiliary Response File.
 ***************************************************************************/
GArf::GArf(const GArf& arf)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(arf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GArf::~GArf(void)
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
 * @param[in] arf Auxiliary Response File.
 * @return Auxiliary Response File.
 ***************************************************************************/
GArf& GArf::operator=(const GArf& arf)
{
    // Execute only if object is not identical
    if (this != &arf) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(arf);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object.
 *
 * Reset object to a clean initial state.
 ***************************************************************************/
void GArf::clear(void)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 *
 * @return Pointer to Auxiliary Response File.
 ***************************************************************************/
GArf* GArf::clone(void) const
{
    // Clone client
    return new GArf(*this);
}


/***********************************************************************//**
 * @brief Return content of spectral bin
 *
 * @param[in] index Bin index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Spectral bin index is out of range.
 *
 * Returns reference to content of spectral bin with specified @p index.
 ***************************************************************************/
double& GArf::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return (m_specresp[index]);
}


/***********************************************************************//**
 * @brief Return content of spectral bin (const version)
 *
 * @param[in] index Bin index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Spectral bin index is out of range.
 *
 * Returns reference to content of spectral bin with specified @p index.
 ***************************************************************************/
const double& GArf::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return (m_specresp[index]);
}


/***********************************************************************//**
 * @brief Load Auxiliary Response File
 *
 * @param[in] filename File name.
 ***************************************************************************/
void GArf::load(const std::string& filename)
{
    // Clear any existing models
    clear();

    // Open FITS file
    GFits file(filename);

    // Get ARF table
    GFitsTable* table = file.table("SPECRESP");

    // Read ARF data
    read(table);

    // Close FITS file
    file.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Auxiliary Response File
 *
 * @param[in] filename File name.
 * @param[in] clobber Overwrite existing file? (defaults to true)
 ***************************************************************************/
void GArf::save(const std::string& filename, const bool& clobber) const
{
    // Open FITS file
    GFits fits;

    // Write ARF into file
    write(fits);

    // Close FITS file
    fits.saveto(filename, clobber);

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Auxiliary Response File
 *
 * @param[in] hdu ARF FITS table.
 ***************************************************************************/
void GArf::read(const GFitsTable* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Get pointer to data columns
        const GFitsTableCol* energy_lo = &(*hdu)["ENERG_LO"];
        const GFitsTableCol* energy_hi = &(*hdu)["ENERG_HI"];
        const GFitsTableCol* specresp  = &(*hdu)["SPECRESP"];

        // Determine effective area conversion factor. Internal
        // units are cm^2
        std::string u_specresp = gammalib::tolower(gammalib::strip_whitespace(specresp->unit()));
        double      c_specresp = 1.0;
        if (u_specresp == "m^2" || u_specresp == "m2") {
            c_specresp = 10000.0;
        }

        // Extract number of energy bins
        int num = energy_lo->length();

        // Set energy bins
        for (int i = 0; i < num; ++i) {
    
            // Append energy bin
            GEnergy emin(energy_lo->real(i), energy_lo->unit());
            GEnergy emax(energy_hi->real(i), energy_hi->unit());
            m_ebounds.append(emin, emax);

            // Append effective area value
            double aeff = specresp->real(i) * c_specresp;
            m_specresp.push_back(aeff);

        } // endfor: looped over energy bins

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Auxiliary Response File
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GArf::write(GFits& fits) const
{
    // Set column length
    int length = size();

    // Continue only if there are bins
    if (length > 0) {

        // Create new binary table
        GFitsBinTable* hdu = new GFitsBinTable;

        // Allocate floating point vector columns
        GFitsTableFloatCol energy_lo("ENERG_LO", length);
        GFitsTableFloatCol energy_hi("ENERG_HI", length);
        GFitsTableFloatCol specresp("SPECRESP",  length);

        // Fill columns
        for (int i = 0; i < length; ++i) {
            energy_lo(i) = m_ebounds.emin(i).keV();
            energy_hi(i) = m_ebounds.emax(i).keV();
            specresp(i)  = m_specresp[i];
        }

        // Set column units
        energy_lo.unit("keV");
        energy_hi.unit("keV");
        specresp.unit("cm^2");

        // Set table attributes
        hdu->extname("SPECRESP");

        // Append columns to table
        hdu->append_column(energy_lo);
        hdu->append_column(energy_hi);
        hdu->append_column(specresp);

        // Append HDU to FITS file
        fits.append(*hdu);

        // Free binary table
        delete hdu;

    } // endif: there were data to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Auxiliary Response File
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing Auxiliary Response File information.
 ***************************************************************************/
std::string GArf::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GArf ===");

        // Append energy boundary information
        result.append("\n"+gammalib::parformat("Number of bins"));
        result.append(gammalib::str(m_ebounds.size()));
        result.append("\n"+gammalib::parformat("Energy range"));
        result.append(m_ebounds.emin().print());
        result.append(" - ");
        result.append(m_ebounds.emax().print());

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
void GArf::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_ebounds.clear();
    m_specresp.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] arf Auxiliary Response File.
 ***************************************************************************/
void GArf::copy_members(const GArf& arf)
{
    // Copy members
    m_filename = arf.m_filename;
    m_ebounds  = arf.m_ebounds;
    m_specresp = arf.m_specresp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GArf::free_members(void)
{
    // Return
    return;
}
