/***************************************************************************
 *                GPha.cpp - XSPEC Pulse Height Analyzer class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2017 by Juergen Knoedlseder                         *
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
 * @file GPha.cpp
 * @brief XSPEC Pulse Height Analyzer class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GPha.hpp"
#include "GException.hpp"
#include "GTools.hpp"
#include "GEnergy.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                                 "GPha::at(int&)"
#define G_READ                                      "GPha::read(GFitsTable*)"

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
GPha::GPha(void)
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
GPha::GPha(const GFilename& filename)
{
    // Initialise members
    init_members();

    // Load PHA file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy boundary constructor
 *
 * @param[in] ebds Energy boundaries.
 ***************************************************************************/
GPha::GPha(const GEbounds& ebds)
{
    // Initialise members
    init_members();

    // Set energy boundaries
    m_ebounds = ebds;

    // Initialize spectrum
    m_counts.assign(ebds.size(), 0.0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy bins constructor
 *
 * @param[in] bins Number of energy bins.
 ***************************************************************************/
GPha::GPha(const int& bins)
{
    // Initialise members
    init_members();

    // Initialize spectrum
    m_counts.assign(bins, 0.0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pha Pulse Height Analyzer spectrum.
 ***************************************************************************/
GPha::GPha(const GPha& pha)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(pha);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPha::~GPha(void)
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
 * @param[in] pha Pulse Height Analyzer spectrum.
 * @return Pulse Height Analyzer spectrum.
 ***************************************************************************/
GPha& GPha::operator=(const GPha& pha)
{
    // Execute only if object is not identical
    if (this != &pha) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(pha);

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
void GPha::clear(void)
{
    // Free memory
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 *
 * @return Pulse Height Analyzer spectrum.
 ***************************************************************************/
GPha* GPha::clone(void) const
{
    // Clone client
    return new GPha(*this);
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
double& GPha::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return (m_counts[index]);
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
const double& GPha::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return (m_counts[index]);
}


/***********************************************************************//**
 * @brief Number of counts in spectrum
 *
 * Returns the number of counts in the spectrum.
 ***************************************************************************/
double GPha::counts(void) const
{
    // Initialise counts
    double counts = 0.0;
    
    // Compute content
    for (int i = 0; i < m_counts.size(); ++i) {
        counts += m_counts[i];
    }

    // Return counts
    return counts;
}


/***********************************************************************//**
 * @brief Fill spectrum with a value.
 *
 * @param[in] energy Energy.
 * @param[in] value Fill value (defaults to 1.0).
 *
 * Fills the specified @p value at a given @p energy in the spectrum.
 ***************************************************************************/
void GPha::fill(const GEnergy& energy, const double& value)
{
    // Get index
    int index = m_ebounds.index(energy);

    // If index is not valid, then check for underflow, overflow or outflow
    if (index == -1) {
        if (energy < m_ebounds.emin()) {
            m_underflow += value;
        }
        else if (energy >= m_ebounds.emax()) {
            m_overflow += value;
        }
        else {
            m_outflow += value;
        }
    }

    // ... otherwise fill the histogram
    else {
        m_counts[index] += value;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Pulse Height Analyzer spectrum
 *
 * @param[in] filename File name.
 *
 * Loads the Pulse Height Analyzer spectrum from the `SPECTRUM` extension
 * of the FITS file. If the file contains also an `EBOUNDS` extension the
 * energy boundaries of all Pulse Height Analyzer channels are also loaded.
 ***************************************************************************/
void GPha::load(const GFilename& filename)
{
    // Clear spectrum
    clear();

    // Open FITS file (without extension name as the user is not allowed
    // to modify the extension names)
    GFits fits(filename.url());

    // Get PHA table
    const GFitsTable& pha = *fits.table(gammalib::extname_pha);

    // Read PHA data
    read(pha);

    // Optionally read energy boundaries
    if (fits.contains(gammalib::extname_ebounds)) {

        // Get energy boundary table
        const GFitsTable& ebounds = *fits.table(gammalib::extname_ebounds);

        // Read energy boundaries
        m_ebounds.read(ebounds);

    } // endif: had energy boundary table

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename.url();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Pulse Height Analyzer spectrum
 *
 * @param[in] filename File name.
 * @param[in] clobber Overwrite existing file?
 *
 * Saves the Pulse Height Analyzer spectrum and energy boundaries into a
 * FITS file. If a file with the given @p filename does not yet exist it
 * will be created. If the file exists it can be overwritten if the
 * @p clobber flag is set to `true`. Otherwise an exception is thrown.
 *
 * The method will save two binary FITS tables into the FITS file: a
 * `SPECTRUM` extension that contains the channel values of the Pulse
 * Height Analyzer spectrum and an `EBOUNDS` extension that contains the
 * energy boundaries for all channels.
 ***************************************************************************/
void GPha::save(const GFilename& filename, const bool& clobber) const
{
    // Create FITS file
    GFits fits;

    // Write PHA into file
    write(fits);

    // Save to file (without extension name since the requested extension
    // may not yet exist in the file)
    fits.saveto(filename.url(), clobber);

    // Store filename
    m_filename = filename.url();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Pulse Height Analyzer spectrum
 *
 * @param[in] table FITS table.
 *
 * @exception GException::invalid_value
 *            Mismatch between PHA file and energy boundaries.
 *
 * Reads the Pulse Height Analyzer spectrum from a FITS table. The channel
 * values are expected in the `COUNTS` column of the table. All other
 * columns are ignored.
 *
 * See
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html
 * for details about the Pulse Height Analyzer spectrum format.
 ***************************************************************************/
void GPha::read(const GFitsTable& table)
{
    // Clear spectrum
    clear();

    // Get data column
    const GFitsTableCol* col_data = table["COUNTS"];

    // Extract number of channels in FITS file
    int length = col_data->nrows();

    // Check whether column length is consistent with energy boundaries
    if (m_ebounds.size() > 0) {
        if (m_ebounds.size() != length) {
            std::string msg = "Mismatch between the "+gammalib::str(length)+
                              " channels in the PHA file and the "+
                              gammalib::str(m_ebounds.size())+" energy "
                              "boundaries. Please correct either the energy "
                              "boundaris or the PHA file.";
            throw GException::invalid_value(G_READ, msg);
        }
    }

    // Initialize spectrum
    m_counts.assign(length, 0.0);

    // Copy data
    for (int i = 0; i < length; ++i) {
        m_counts[i] = col_data->real(i);
    }

    // Read keywords
    m_exposure = (table.has_card("EXPOSURE")) ? table.real("EXPOSURE") : 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Pulse Height Analyzer spectrum
 *
 * @param[in] fits FITS file.
 *
 * Writes the Pulse Height Analyzer spectrum into `SPECTRUM` and `EBOUNDS`
 * extensions of the FITS file. Extensions with these names will be removed
 * from the FITS file before writing.
 *
 * The columns `CHANNEL`, `COUNTS`, `STAT_ERR`, `SYS_ERR`, `QUALITY`,
 * `GROUPING`, `AREASCAL`, and `BACKSCAL` will be written into the `SPECTRUM`
 * extension, but only the `CHANNEL` and `COUNTS` columns will be filled with
 * values. Note that the channels start from 1 in the Pulse Height Analyzer
 * spectrum.
 *
 * See
 * https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html
 * for details about the PHA file format.
 ***************************************************************************/
void GPha::write(GFits& fits) const
{
    // Remove extensions if they exist already
    if (fits.contains(gammalib::extname_ebounds)) {
        fits.remove(gammalib::extname_ebounds);
    }
    if (fits.contains(gammalib::extname_pha)) {
        fits.remove(gammalib::extname_pha);
    }

    // Set column length
    int length = size();

    // Continue only if there are bins
    if (length > 0) {

        // Create new binary table
        GFitsBinTable hdu;

        // Allocate floating point vector columns
        GFitsTableShortCol col_chan("CHANNEL",  length);
        GFitsTableFloatCol col_data("COUNTS",   length);
        GFitsTableFloatCol col_stat("STAT_ERR", length);
        GFitsTableFloatCol col_syst("SYS_ERR",  length);
        GFitsTableShortCol col_qual("QUALITY",  length);
        GFitsTableShortCol col_grpg("GROUPING", length);
        GFitsTableFloatCol col_area("AREASCAL", length);
        GFitsTableFloatCol col_back("BACKSCAL", length);

        // Fill columns
        for (int i = 0; i < length; ++i) {
            col_chan(i) = i+1; // Channels start at 1
            col_data(i) = float(m_counts[i]);
        }

        // Set table attributes
        hdu.extname(gammalib::extname_pha);

        // Append columns to table
        hdu.append(col_chan);
        hdu.append(col_data);
        hdu.append(col_stat);
        hdu.append(col_syst);
        hdu.append(col_qual);
        hdu.append(col_grpg);
        hdu.append(col_area);
        hdu.append(col_back);

        // Write keywords
        hdu.card("EXPOSURE", m_exposure, "[s] Deadtime corrected exposure time");

        // Append HDU to FITS file
        fits.append(hdu);

        // Optionally append energy boundaries
        if (m_ebounds.size() > 0) {
            m_ebounds.write(fits);
        }

    } // endif: there were data to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Pulse Height Analyzer spectrum
 *
 * @param[in] chatter Chattiness.
 * @return String containing Pulse Height Analyzer spectrum information.
 ***************************************************************************/
std::string GPha::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GPha ===");

        // Append energy boundary information
        result.append("\n"+gammalib::parformat("Number of bins"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Energy range"));
        if (m_ebounds.size() > 0) {
            result.append(m_ebounds.emin().print());
            result.append(" - ");
            result.append(m_ebounds.emax().print());
        }
        else {
            result.append("not specified");
        }

        // Append content information
        result.append("\n"+gammalib::parformat("Total number of counts"));
        result.append(gammalib::str(counts()));
        result.append("\n"+gammalib::parformat("Underflow counts"));
        result.append(gammalib::str(m_underflow));
        result.append("\n"+gammalib::parformat("Overflow counts"));
        result.append(gammalib::str(m_overflow));
        result.append("\n"+gammalib::parformat("Outflow counts"));
        result.append(gammalib::str(m_outflow));

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
void GPha::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_ebounds.clear();
    m_counts.clear();
    m_underflow = 0.0;
    m_overflow  = 0.0;
    m_outflow   = 0.0;
    m_exposure  = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pha Pulse Height Analyzer spectrum.
 ***************************************************************************/
void GPha::copy_members(const GPha& pha)
{
    // Copy members
    m_filename  = pha.m_filename;
    m_counts    = pha.m_counts;
    m_underflow = pha.m_underflow;
    m_overflow  = pha.m_overflow;
    m_outflow   = pha.m_outflow;
    m_ebounds   = pha.m_ebounds;
    m_exposure  = pha.m_exposure;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPha::free_members(void)
{
    // Return
    return;
}
