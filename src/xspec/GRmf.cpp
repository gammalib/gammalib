/***************************************************************************
 *            GRmf.cpp - XSPEC Redistribution Matrix File class            *
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
 * @file GRmf.cpp
 * @brief XSPEC Redistribution Matrix File class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GRmf.hpp"
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                           "GRmf::at(int&, int&)"

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
GRmf::GRmf(void)
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
GRmf::GRmf(const std::string& filename)
{
    // Initialise members
    init_members();

    // Load RMF file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy boundary constructor
 *
 * @param[in] etrue True energy boundaries.
 * @param[in] emeasured Measured energy boundaries.
 ***************************************************************************/
GRmf::GRmf(const GEbounds& etrue, const GEbounds& emeasured)
{
    // Initialise members
    init_members();

    // Set energy boundaries
    m_ebds_true     = etrue;
    m_ebds_measured = emeasured;

    // Initialize matrix
    m_matrix = GMatrixSparse(m_ebds_true.size(), m_ebds_measured.size());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rmf Redistribution Matrix File.
 ***************************************************************************/
GRmf::GRmf(const GRmf& rmf)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(rmf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GRmf::~GRmf(void)
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
 * @param[in] rmf Redistribution Matrix File.
 * @return Redistribution Matrix File.
 ***************************************************************************/
GRmf& GRmf::operator=(const GRmf& rmf)
{
    // Execute only if object is not identical
    if (this != &rmf) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(rmf);

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
void GRmf::clear(void)
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
 * @return Redistribution Matrix File.
 ***************************************************************************/
GRmf* GRmf::clone(void) const
{
    // Clone object
    return new GRmf(*this);
}


/***********************************************************************//**
 * @brief Return content of redistribution matrix bin
 *
 * @param[in] itrue True energy index [0,...,ntrue()-1].
 * @param[in] imeasured Measured energy index [0,...,nmeasured()-1].
 *
 * @exception GException::out_of_range
 *            Bin index is out of range.
 *
 * Returns reference to content of redistribution matrix bin bin with true
 * energy index @p itrue and measured energy index @p imeasured.
 ***************************************************************************/
double& GRmf::at(const int& itrue, const int& imeasured)
{
    // Raise exception if indices are out of range
    if (itrue < 0 || itrue >= ntrue()) {
        throw GException::out_of_range(G_AT, itrue, 0, ntrue()-1);
    }
    // Raise exception if indices are out of range
    if (imeasured < 0 || imeasured >= nmeasured()) {
        throw GException::out_of_range(G_AT, imeasured, 0, nmeasured()-1);
    }

    // Return reference
    return (m_matrix(itrue, imeasured));
}


/***********************************************************************//**
 * @brief Return content of redistribution matrix bin (const version)
 *
 * @param[in] itrue True energy index [0,...,ntrue()-1].
 * @param[in] imeasured Measured energy index [0,...,nmeasured()-1].
 *
 * @exception GException::out_of_range
 *            Bin index is out of range.
 *
 * Returns reference to content of redistribution matrix bin bin with true
 * energy index @p itrue and measured energy index @p imeasured.
 ***************************************************************************/
const double& GRmf::at(const int& itrue, const int& imeasured) const
{
    // Raise exception if indices are out of range
    if (itrue < 0 || itrue >= ntrue()) {
        throw GException::out_of_range(G_AT, itrue, 0, ntrue()-1);
    }
    // Raise exception if indices are out of range
    if (imeasured < 0 || imeasured >= nmeasured()) {
        throw GException::out_of_range(G_AT, imeasured, 0, nmeasured()-1);
    }

    // Return reference
    return (m_matrix(itrue, imeasured));
}


/***********************************************************************//**
 * @brief Load Redistribution Matrix File
 *
 * @param[in] filename File name.
 ***************************************************************************/
void GRmf::load(const std::string& filename)
{
    // Clear any existing models
    clear();

    // Open FITS file
    GFits file(filename);

    // Read measured energy boundaries
    m_ebds_measured.load(filename);

    // Get RMF table
    GFitsTable* table = file.table("MATRIX");

    // Read RMF data
    read(table);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Redistribution Matrix File
 *
 * @param[in] filename File name.
 * @param[in] clobber Overwrite existing file? (defaults to true)
 ***************************************************************************/
void GRmf::save(const std::string& filename, const bool& clobber) const
{
    // Open FITS file
    GFits fits;

    // Write RMF into file
    write(fits);

    // Close FITS file
    fits.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Redistribution Matrix File
 *
 * @param[in] hdu RMF FITS table.
 ***************************************************************************/
void GRmf::read(const GFitsTable* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Initialise members
        m_ebds_true.clear();

        // Get pointer to data columns
        const GFitsTableCol* energy_lo = &(*hdu)["ENERG_LO"];
        const GFitsTableCol* energy_hi = &(*hdu)["ENERG_HI"];
        const GFitsTableCol* n_grp     = &(*hdu)["N_GRP"];
        const GFitsTableCol* f_chan    = &(*hdu)["F_CHAN"];
        const GFitsTableCol* n_chan    = &(*hdu)["N_CHAN"];
        const GFitsTableCol* matrix    = &(*hdu)["MATRIX"];

        // Set matrix rows and columns
        int rows    = energy_lo->length();
        int columns = m_ebds_measured.size();

        // Initialize matrix
        m_matrix = GMatrixSparse(rows, columns);

        // Set true energy bins
        for (int itrue = 0; itrue < rows; ++itrue) {
    
            // Append energy bin
            GEnergy emin(energy_lo->real(itrue), energy_lo->unit());
            GEnergy emax(energy_hi->real(itrue), energy_hi->unit());
            m_ebds_true.append(emin, emax);

            // Loop over groups
            int icolumn = 0;
            int ngroups = n_grp->integer(itrue);
            for (int igroup = 0; igroup < ngroups; ++igroup) {

                // Get start column index and number of columns
                int imeasured = f_chan->integer(itrue, igroup);
                int nvalues   = n_chan->integer(itrue, igroup);

                // Get values
                for (int i = 0; i < nvalues; ++i, ++imeasured, ++icolumn) {
                    m_matrix(itrue, imeasured) = matrix->real(itrue, icolumn);
                }

            } // endfor: looped over groups

        } // endfor: looped over true energy bins

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Redistribution Matrix File
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GRmf::write(GFits& fits) const
{
    // Set table length
    int length = ntrue();

    // Continue only if there are true energies
    if (length > 0) {

        // Compute maximum number of matrix elements
        int maxcolumns = 1;
        for (int i = 0; i < length; i++) {
            
            // Search first non-zero element
            int ifirst = 0;
            for (; ifirst < nmeasured(); ++ifirst) {
                if (m_matrix(i, ifirst) != 0.0) {
                    break;
                }
            }

            // Search last non-zero element
            int ilast = nmeasured() - 1;
            for (; ilast >= 0; --ilast) {
                if (m_matrix(i, ilast) != 0.0) {
                    break;
                }
            }

            // Compute number of non-zero matrix elements
            int num = ilast - ifirst + 1;

            // Set maximum
            if (num > maxcolumns) {
                maxcolumns = num;
            }
        
        } // endfor: looped over true energies

        // Create new binary table
        GFitsBinTable* hdu = new GFitsBinTable;

        // Allocate columns
        GFitsTableFloatCol energy_lo("ENERG_LO", length);
        GFitsTableFloatCol energy_hi("ENERG_HI", length);
        GFitsTableShortCol n_grp("N_GRP", length);
        GFitsTableShortCol f_chan("F_CHAN", length);
        GFitsTableShortCol n_chan("N_CHAN", length);
        GFitsTableFloatCol matrix("MATRIX", length, maxcolumns);

        // Fill columns
        for (int itrue = 0; itrue < length; ++itrue) {

            // Search first non-zero element
            int ifirst = 0;
            for (; ifirst < nmeasured(); ++ifirst) {
                if (m_matrix(itrue, ifirst) != 0.0) {
                    break;
                }
            }

            // Search last non-zero element
            int ilast = nmeasured() - 1;
            for (; ilast >= 0; --ilast) {
                if (m_matrix(itrue, ilast) != 0.0) {
                    break;
                }
            }

            // Compute number of non-zero matrix elements
            int num = ilast - ifirst + 1;
            if (num < 0) {
                ifirst = 0;
                num    = 0;
            }
    
            // Set energy bin
            energy_lo(itrue) = m_ebds_true.emin(itrue).keV();
            energy_hi(itrue) = m_ebds_true.emax(itrue).keV();

            // Set group information
            n_grp(itrue)  = 1;
            f_chan(itrue) = ifirst;
            n_chan(itrue) = num;

            // Set matrix elements
            for (int i = 0; i < num; ++i, ++ifirst) {
                matrix(itrue, i) = m_matrix(itrue, ifirst);
            }

        } // endfor: looped over true energy bins

        // Set column units
        energy_lo.unit("keV");
        energy_hi.unit("keV");

        // Set table attributes
        hdu->extname("MATRIX");

        // Append columns to table
        hdu->append_column(energy_lo);
        hdu->append_column(energy_hi);
        hdu->append_column(n_grp);
        hdu->append_column(f_chan);
        hdu->append_column(n_chan);
        hdu->append_column(matrix);

        // Append HDU to FITS file
        fits.append(*hdu);

        // Free binary table
        delete hdu;

        // Append measured energy boundaries
        m_ebds_measured.write(&fits);

    } // endif: there were data to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Redistribution Matrix File
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing Redistribution Matrix File information.
 ***************************************************************************/
std::string GRmf::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GRmf ===");

        // Append energy boundary information
        result.append("\n"+gammalib::parformat("Number of true energy bins"));
        result.append(gammalib::str(m_ebds_true.size()));
        result.append("\n"+gammalib::parformat("Number of measured bins"));
        result.append(gammalib::str(m_ebds_measured.size()));
        result.append("\n"+gammalib::parformat("True energy range"));
        result.append(m_ebds_true.emin().print());
        result.append(" - ");
        result.append(m_ebds_true.emax().print());
        result.append("\n"+gammalib::parformat("Measured energy range"));
        result.append(m_ebds_measured.emin().print());
        result.append(" - ");
        result.append(m_ebds_measured.emax().print());

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
void GRmf::init_members(void)
{
    // Initialise members
    m_ebds_true.clear();
    m_ebds_measured.clear();
    m_matrix.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rmf Redistribution Matrix File.
 ***************************************************************************/
void GRmf::copy_members(const GRmf& rmf)
{
    // Copy members
    m_ebds_true     = rmf.m_ebds_true;
    m_ebds_measured = rmf.m_ebds_measured;
    m_matrix        = rmf.m_matrix;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GRmf::free_members(void)
{
    // Return
    return;
}
