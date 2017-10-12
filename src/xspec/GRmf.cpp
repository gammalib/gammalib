/***************************************************************************
 *            GRmf.cpp - XSPEC Redistribution Matrix File class            *
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
 * @file GRmf.cpp
 * @brief XSPEC Redistribution Matrix File class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GRmf.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR_PLUS                             "GRmf::operator+=(GRmf&)"
#define G_OPERATOR_MINUS                            "GRmf::operator-=(GRmf&)"
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
GRmf::GRmf(const GFilename& filename)
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


/***********************************************************************//**
 * @brief Add Redistribution Matrix File
 *
 * @param[in] rmf Redistribution Matrix File.
 * @return Sum of Redistribution Matrix File.
 *
 * @exception GException::invalid_value
 *            Incompatible Redistribution Matrix Files.
 *
 * Adds the RMF values of an Redistribution Matrix File to the current
 * values.
 *
 * The operator only works if the provide Redistribution Matrix File has the
 * same energy binning than the current Redistribution Matrix File.
 ***************************************************************************/
GRmf& GRmf::operator+=(const GRmf& rmf)
{
    // Throw an exception if the RMF are not compatible
    if (this->etrue() != rmf.etrue()) {
        std::string msg = "Incompatible true energy binning of "
                          "Redistribution Matrix File.";
        throw GException::invalid_value(G_OPERATOR_PLUS, msg);
    }
    if (this->emeasured() != rmf.emeasured()) {
        std::string msg = "Incompatible measured energy binning of "
                          "Redistribution Matrix File.";
        throw GException::invalid_value(G_OPERATOR_PLUS, msg);
    }
    
    // Add RMF values
    for (int itrue = 0; itrue < ntrue(); ++itrue) {
        for (int imeasured = 0; imeasured < nmeasured(); ++imeasured) {
            m_matrix(itrue, imeasured) += rmf.m_matrix(itrue, imeasured);
        }
    }

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Subtract Redistribution Matrix File
 *
 * @param[in] rmf Redistribution Matrix File.
 * @return Difference of Redistribution Matrix File.
 *
 * @exception GException::invalid_value
 *            Incompatible Redistribution Matrix Files.
 *
 * Subtracts the RMF values of an Redistribution Matrix File from the current
 * values.
 *
 * The operator only works if the provide Redistribution Matrix File has the
 * same energy binning than the current Redistribution Matrix File.
 ***************************************************************************/
GRmf& GRmf::operator-=(const GRmf& rmf)
{
    // Throw an exception if the RMF are not compatible
    if (this->etrue() != rmf.etrue()) {
        std::string msg = "Incompatible true energy binning of "
                          "Redistribution Matrix File.";
        throw GException::invalid_value(G_OPERATOR_MINUS, msg);
    }
    if (this->emeasured() != rmf.emeasured()) {
        std::string msg = "Incompatible measured energy binning of "
                          "Redistribution Matrix File.";
        throw GException::invalid_value(G_OPERATOR_MINUS, msg);
    }
    
    // Subtract RMF values
    for (int itrue = 0; itrue < ntrue(); ++itrue) {
        for (int imeasured = 0; imeasured < nmeasured(); ++imeasured) {
            m_matrix(itrue, imeasured) -= rmf.m_matrix(itrue, imeasured);
        }
    }

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Scale Redistribution Matrix File values
 *
 * @param[in] scale Scale factor.
 * @return Scaled Redistribution Matrix File.
 *
 * Multiplies the values of the Redistribution Matrix File with a scale
 * factor.
 ***************************************************************************/
GRmf& GRmf::operator*=(const double& scale)
{
    // Scale RMF values
    for (int itrue = 0; itrue < ntrue(); ++itrue) {
        for (int imeasured = 0; imeasured < nmeasured(); ++imeasured) {
            m_matrix(itrue, imeasured) *= scale;
        }
    }

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Divide Redistribution Matrix File values
 *
 * @param[in] scale Division factor.
 * @return Divided Redistribution Matrix File.
 *
 * Divides the values of the Redistribution Matrix File by a division factor.
 ***************************************************************************/
GRmf& GRmf::operator/=(const double& scale)
{
    // Divide RMF values
    for (int itrue = 0; itrue < ntrue(); ++itrue) {
        for (int imeasured = 0; imeasured < nmeasured(); ++imeasured) {
            m_matrix(itrue, imeasured) /= scale;
        }
    }

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
 * @brief Return true energy boundaries for specified measured energy
 *
 * @param[in] emeasured Measured energy.
 * @return True energy boundaries for specified measured energy.
 *
 * Returns the true energy boundaries for the specified measured energy
 * @p emeasured. If the RMF is not covering the specified measured energy,
 * an empty energy boundary object is returned.
 ***************************************************************************/
GEbounds GRmf::etrue(const GEnergy& emeasured) const
{
    // Initialise empty energy boundaries
    GEbounds ebounds;

    // Determine matrix column that corresponds to specified measured energy
    int column = m_ebds_measured.index(emeasured);

    // If matrix column was found then determine the boundaries for this
    // column
    if (column != -1) {

        // Determine first and last non-zero indices
        int row_start = -1;
        int row_stop  = -1;
        for (int row = 0; row < m_matrix.rows(); ++row) {
            if (m_matrix(row, column) > 0.0) {
                row_start = row;
                break;
            }
        }
        for (int row = m_matrix.rows()-1; row >= 0; --row) {
            if (m_matrix(row, column) > 0.0) {
                row_stop = row;
                break;
            }
        }

        // Set energy boundaries if valid indices have been found
        if (row_start != -1 && row_stop != -1) {
            ebounds = GEbounds(m_ebds_true.emin(row_start),
                               m_ebds_true.emax(row_stop));
        }

    } // endif: row containing true energy was found

    // Return energy boundaries
    return ebounds;
}
//m_matrix(itrue, imeasured) (row, column)


/***********************************************************************//**
 * @brief Return measured energy boundaries for specified true energy
 *
 * @param[in] etrue True energy.
 * @return Measured energy boundaries for specified true energy.
 *
 * Returns the measured energy boundaries for the specified true energy
 * @p etrue. If the RMF is not covering the specified true energy,
 * an empty energy boundary object is returned.
 ***************************************************************************/
GEbounds GRmf::emeasured(const GEnergy& etrue) const
{
    // Initialise empty energy boundaries
    GEbounds ebounds;

    // Determine matrix row that corresponds to specified true energy
    int row = m_ebds_true.index(etrue);

    // If matrix row was found then determine the boundaries for this
    // row
    if (row != -1) {

        // Determine first and last non-zero indices
        int column_start = -1;
        int column_stop  = -1;
        for (int column = 0; column < m_matrix.columns(); ++column) {
            if (m_matrix(row, column) > 0.0) {
                column_start = column;
                break;
            }
        }
        for (int column = m_matrix.columns()-1; column >= 0; --column) {
            if (m_matrix(row, column) > 0.0) {
                column_stop = column;
                break;
            }
        }

        // Set energy boundaries if valid indices have been found
        if (column_start != -1 && column_stop != -1) {
            ebounds = GEbounds(m_ebds_measured.emin(column_start),
                               m_ebds_measured.emax(column_stop));
        }

    } // endif: row containing true energy was found

    // Return energy boundaries
    return ebounds;
}


/***********************************************************************//**
 * @brief Load Redistribution Matrix File
 *
 * @param[in] filename File name.
 *
 * Loads the Redistribution Matrix File from the `MATRIX` extension of the
 * FITS file. If the file contains also an `EBOUNDS` extension the energy
 * boundaries of all measured energies are also loaded.
 ***************************************************************************/
void GRmf::load(const GFilename& filename)
{
    // Clear matrix
    clear();

    // Open FITS file (without extension name as the user is not allowed
    // to modify the extension names)
    GFits fits(filename.url());

    // Read measured energy boundaries
    m_ebds_measured.load(filename.url());

    // Get RMF table
    const GFitsTable& table = *fits.table(gammalib::extname_rmf);

    // Read RMF data
    read(table);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename.url();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Redistribution Matrix File
 *
 * @param[in] filename File name.
 * @param[in] clobber Overwrite existing file?
 * @param[in] unit Energy unit.
 *
 * Saves the Redistribution Matrix File into a FITS file. If a file with the
 * given @p filename does not yet exist it will be created. If the file
 * exists it can be overwritten if the @p clobber flag is set to `true`.
 * Otherwise an exception is thrown.
 *
 * The method will save two binary FITS tables into the FITS file: a
 * `MATRIX` extension that contains the Redistribution Matrix File elements,
 * and an `EBOUNDS` extension that contains the measured energy boundaries
 * for all matrix. The @p unit argument specifies the unit of the true
 * energies that will be saved in the `MATRIX` extension.
 ***************************************************************************/
void GRmf::save(const GFilename&   filename,
                const bool&        clobber,
                const std::string& unit) const
{
    // Create FITS file
    GFits fits;

    // Write RMF into file
    write(fits, unit);

    // Save to file (without extension name since the requested extension
    // may not yet exist in the file)
    fits.saveto(filename.url(), clobber);

    // Store filename
    m_filename = filename.url();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Redistribution Matrix File
 *
 * @param[in] table RMF FITS table.
 *
 * Reads the Redistribution Matrix File from a FITS table. The true energy
 * bins are expected in the `ENERG_LO` and `ENERG_HI` columns. Information
 * for matrix compression are stored in the `N_GRP`, `F_CHAN` and `N_CHAN`
 * columns, and the matrix elements are stored in the `MATRIX` column.
 *
 * The method automatically decompresses the Redistribution Matrix File and
 * stores the matrix elements in a spare matrix.
 *
 * See
 * http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc3
 * for details about the Redistribution Matrix File file format.
 ***************************************************************************/
void GRmf::read(const GFitsTable& table)
{
    // Initialise members
    m_ebds_true.clear();

    // Get pointer to data columns
    const GFitsTableCol* energy_lo = table["ENERG_LO"];
    const GFitsTableCol* energy_hi = table["ENERG_HI"];
    const GFitsTableCol* n_grp     = table["N_GRP"];
    const GFitsTableCol* f_chan    = table["F_CHAN"];
    const GFitsTableCol* n_chan    = table["N_CHAN"];
    const GFitsTableCol* matrix    = table["MATRIX"];

    // Set matrix rows and columns
    int rows    = energy_lo->nrows();
    int columns = m_ebds_measured.size();

    // Initialize matrix
    m_matrix = GMatrixSparse(rows, columns);

    // Initialize matrix maximum
    double max = 0.0;

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

                // Set matrix value
                m_matrix(itrue, imeasured) = matrix->real(itrue, icolumn);

                // Update fmax
                if (max < m_matrix(itrue, imeasured)) {
                    max = m_matrix(itrue, imeasured);
                    m_itruemax = itrue;
                    m_imeasmax = imeasured;
                }

            } // endfor: looped over measured energy bins

        } // endfor: looped over groups

    } // endfor: looped over true energy bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Redistribution Matrix File
 *
 * @param[in] fits FITS file.
 * @param[in] unit Energy unit.
 *
 * Writes the Redistribution Matrix File into `MATRIX` and `EBOUNDS`
 * extensions of the FITS file. Extensions with these names will be removed
 * from the FITS file before writing.
 *
 * The true energy bins are written into the `ENERG_LO` and `ENERG_HI`
 * columns. The units of the energies in these columns is specified by the
 * @p unit argument. Information for matrix compression are stored in the
 * `N_GRP`, `F_CHAN` and `N_CHAN` columns, and the matrix elements are
 * written into the `MATRIX` column.
 *
 * The measured energy bins are written into the `EBOUNDS` extension.
 *
 * See
 * http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc3
 * for details about the Redistribution Matrix File file format.
 ***************************************************************************/
void GRmf::write(GFits& fits, const std::string& unit) const
{
    // Remove extensions if they exist already
    if (fits.contains(gammalib::extname_rmf)) {
        fits.remove(gammalib::extname_rmf);
    }
    if (fits.contains(gammalib::extname_ebounds)) {
        fits.remove(gammalib::extname_ebounds);
    }

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
        GFitsBinTable hdu;

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
            energy_lo(itrue) = m_ebds_true.emin(itrue)(unit);
            energy_hi(itrue) = m_ebds_true.emax(itrue)(unit);

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
        energy_lo.unit(unit);
        energy_hi.unit(unit);

        // Set table attributes
        hdu.extname(gammalib::extname_rmf);

        // Append columns to table
        hdu.append(energy_lo);
        hdu.append(energy_hi);
        hdu.append(n_grp);
        hdu.append(f_chan);
        hdu.append(n_chan);
        hdu.append(matrix);

        // Append HDU to FITS file
        fits.append(hdu);

        // Append measured energy boundaries
        m_ebds_measured.write(fits);

    } // endif: there were data to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Redistribution Matrix File
 *
 * @param[in] chatter Chattiness.
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
    m_filename.clear();
    m_ebds_true.clear();
    m_ebds_measured.clear();
    m_matrix.clear();
    m_itruemax = 0;
    m_imeasmax = 0;

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
    m_filename      = rmf.m_filename;
    m_ebds_true     = rmf.m_ebds_true;
    m_ebds_measured = rmf.m_ebds_measured;
    m_matrix        = rmf.m_matrix;
    m_itruemax      = rmf.m_itruemax;
    m_imeasmax      = rmf.m_imeasmax;

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
