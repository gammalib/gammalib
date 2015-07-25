/***************************************************************************
 *               GLATEdisp.cpp - Fermi-LAT energy dispersion               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2015 by Juergen Knoedlseder                         *
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
 * @file GLATEdisp.cpp
 * @brief Fermi-LAT energy dispersion class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATEdisp.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                           "GLATEdisp::read(const GFits& file)"
#define G_READ_EDISP                     "GLATEdisp::read_edisp(GFitsTable&)"

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
GLATEdisp::GLATEdisp(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename FITS file name.
 * @param[in] evtype Event type.
 *
 * Construct instance by loading the energy dispersion information from FITS
 * file.
 ***************************************************************************/
GLATEdisp::GLATEdisp(const std::string& filename, const std::string& evtype)
{
    // Initialise class members
    init_members();

    // Load energy dispersion from file
    load(filename, evtype);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] edisp Energy dispersion.
 ***************************************************************************/
GLATEdisp::GLATEdisp(const GLATEdisp& edisp)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(edisp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATEdisp::~GLATEdisp(void)
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
 * @param[in] edisp Energy dispersion.
 * @return Energy dispersion.
 ***************************************************************************/
GLATEdisp& GLATEdisp::operator=(const GLATEdisp& edisp)
{
    // Execute only if object is not identical
    if (this != &edisp) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(edisp);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear energy dispersion response
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATEdisp::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone energy dispersion response
 *
 * @return Pointer to deep copy of energy dispersion response.
 ***************************************************************************/
GLATEdisp* GLATEdisp::clone(void) const
{
    return new GLATEdisp(*this);
}


/***********************************************************************//**
 * @brief Load energy dispersion from FITS file
 *
 * @param[in] filename Energy dispersion file name.
 * @param[in] evtype Event type.
 ***************************************************************************/
void GLATEdisp::load(const std::string& filename, const std::string& evtype)
{
    // Store event type
    m_evtype = evtype;

    // Open FITS file
    GFits fits(filename);

    // Read energy dispersion from file
    read(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy dispersion into FITS file
 *
 * @param[in] filename FITS file.
 * @param[in] clobber Overwrite existing file?.
 ***************************************************************************/
void GLATEdisp::save(const std::string& filename, const bool& clobber)
{
    // Open FITS file
    GFits fits(filename, true);

    // Write energy dispersion into file
    write(fits);

    // Close FITS file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy dispersion from FITS file
 *
 * @param[in] fits FITS file.
 *
 * @todo Implement reading of scaling parameters
 ***************************************************************************/
void GLATEdisp::read(const GFits& fits)
{
    // Clear instance (keep event type)
    std::string evtype = m_evtype;
    clear();
    m_evtype = evtype;

    // Set extension names
    std::string engdisp = "ENERGY DISPERSION";
    std::string escales = "EDISP_SCALING_PARAMS";
    if (!fits.contains(engdisp)) {
        engdisp += "_" + m_evtype;
    }
    if (!fits.contains(escales)) {
        escales += "_" + m_evtype;
    }

    // Get pointer to effective area HDU
    const GFitsTable& hdu_edisp = *fits.table(engdisp);
    const GFitsTable& hdu_scale = *fits.table(escales);

    // Read energy dispersion
    read_edisp(hdu_edisp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energy dispersion into FITS file
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GLATEdisp::write(GFits& fits) const
{
    // Write energy dispersion
    write_edisp(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print energy dispersion information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy dispersion information.
 ***************************************************************************/
std::string GLATEdisp::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

    // Append header
    result.append("=== GLATEdisp ===");

        // Append information
    result.append("\n"+gammalib::parformat("Number of energy bins") +
                  gammalib::str(nenergies()));
    result.append("\n"+gammalib::parformat("Number of cos theta bins") +
                  gammalib::str(ncostheta()));

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
void GLATEdisp::init_members(void)
{
    // Initialise members
    m_evtype.clear();
    m_edisp_bins.clear();
    m_norm.clear();
    m_ls1.clear();
    m_scale.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] edisp Energy dispersion.
 ***************************************************************************/
void GLATEdisp::copy_members(const GLATEdisp& edisp)
{
    // Copy members
    m_evtype     = edisp.m_evtype;
    m_edisp_bins = edisp.m_edisp_bins;
    m_norm       = edisp.m_norm;
    m_ls1        = edisp.m_ls1;
    m_scale      = edisp.m_scale;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEdisp::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy dispersion from FITS table
 *
 * @param[in] hdu FITS table pointer.
 *
 * @exception GLATException::inconsistent_response
 *            Inconsistent response table encountered
 ***************************************************************************/
void GLATEdisp::read_edisp(const GFitsTable& hdu)
{
    // Clear arrays
    m_norm.clear();
    m_ls1.clear();

    // Get energy and cos theta binning
    m_edisp_bins.read(hdu);

    // Continue only if there are bins
    int size = m_edisp_bins.size();
    if (size > 0) {
/*
        // Allocate arrays
        m_norm.reserve(size);
        m_ls1.reserve(size);

        // Get pointer to columns
        const GFitsTableCol* norm = hdu["NORM"];
        const GFitsTableCol* ls1  = hdu["LS1"];

        // Check consistency of columns
        if (norm->number() != size) {
            throw GLATException::inconsistent_response(G_READ_EDISP,
                                                       norm->number(), size);
        }
        if (ls1->number() != size) {
            throw GLATException::inconsistent_response(G_READ_EDISP,
                                                       ls1->number(), size);
        }

        // Copy data
        for (int i = 0; i < size; ++i) {
            m_norm.push_back(norm->real(0,i));
            m_ls1.push_back(ls1->real(0,i));
        }
*/
    } // endif: there were bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energy dispersion into FITS file
 *
 * @param[in] file FITS file.
 *
 * This method does not write anything if the instance is empty.
 ***************************************************************************/
void GLATEdisp::write_edisp(GFits& file) const
{
    // Continue only if there are bins
    int size = m_edisp_bins.size();
    if (size > 0) {

        // Create new binary table
        GFitsBinTable* hdu_edisp = new GFitsBinTable;

        // Set table attributes
        hdu_edisp->extname("ENERGY DISPERSION");

        // Write boundaries into table
        m_edisp_bins.write(*hdu_edisp);

        // Allocate floating point vector columns
        GFitsTableFloatCol col_norm = GFitsTableFloatCol("NORM",  1, size);
        GFitsTableFloatCol col_ls1  = GFitsTableFloatCol("LS1",   1, size);

        // Fill columns
        for (int i = 0; i < size; ++i) {
            col_norm(0,i) = m_norm[i];
            col_ls1(0,i)  = m_ls1[i];
        }

        // Append columns to table
        hdu_edisp->append(col_norm);
        hdu_edisp->append(col_ls1);

        // Append HDU to FITS file
        file.append(*hdu_edisp);

        // Free table
        delete hdu_edisp;

    } // endif: there were data to write

    // Return
    return;
}
