/***************************************************************************
 *               GArf.cpp - XSPEC Auxiliary Response File class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2018 by Juergen Knoedlseder                         *
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
#include "GFilename.hpp"
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR_PLUS                             "GArf::operator+=(GArf&)"
#define G_OPERATOR_MINUS                            "GArf::operator-=(GArf&)"
#define G_OPERATOR                           "GArf::operator[](std::string&)"
#define G_OPERATOR2                "GArf::operator()(std::string&, GEnergy&)"
#define G_AT1                                                "GArf::at(int&)"
#define G_AT2                                          "GArf::at(int&, int&)"
#define G_APPEND           "GArf::append(std::string&, std::vector<double>&)"

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
GArf::GArf(const GFilename& filename)
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

    // Set log true energy node array
    set_logetrue();

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


/***********************************************************************//**
 * @brief Add Auxiliary Response File
 *
 * @param[in] arf Auxiliary Response File.
 * @return Sum of Auxiliary Response Files.
 *
 * @exception GException::invalid_value
 *            Incompatible Auxiliary Response Files.
 *
 * Adds the ARF values of an Auxiliary Response File to the current values.
 *
 * The operator only works if the provide Auxiliary Response File has the
 * same energy binning than the current Auxiliary Response File.
 ***************************************************************************/
GArf& GArf::operator+=(const GArf& arf)
{
    // Throw an exception if the ARF are not compatible
    if (this->ebounds() != arf.ebounds()) {
        std::string msg = "Incompatible energy binning of Auxiliary "
                          "Response File.";
        throw GException::invalid_value(G_OPERATOR_PLUS, msg);
    }
    
    // Add ARF values
    for (int i = 0; i < this->size(); ++i) {
        m_specresp[i] += arf.m_specresp[i];
    }

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Subtract Auxiliary Response File
 *
 * @param[in] arf Auxiliary Response File.
 * @return Difference of Auxiliary Response File.
 *
 * @exception GException::invalid_value
 *            Incompatible Auxiliary Response Files.
 *
 * Subtracts the ARF values of an Auxiliary Response File from the current
 * values.
 *
 * The operator only works if the provide Auxiliary Response File has the
 * same energy binning than the current Auxiliary Response File.
 ***************************************************************************/
GArf& GArf::operator-=(const GArf& arf)
{
    // Throw an exception if the ARF are not compatible
    if (this->ebounds() != arf.ebounds()) {
        std::string msg = "Incompatible energy binning of Auxiliary "
                          "Response File.";
        throw GException::invalid_value(G_OPERATOR_MINUS, msg);
    }
    
    // Subtract ARF values
    for (int i = 0; i < this->size(); ++i) {
        m_specresp[i] -= arf.m_specresp[i];
    }

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Scale Auxiliary Response File values
 *
 * @param[in] scale Scale factor.
 * @return Scaled Auxiliary Response File.
 *
 * Multiplies the values of the Auxiliary Response File with a scale factor. 
 ***************************************************************************/
GArf& GArf::operator*=(const double& scale)
{
    // Scale ARF values
    for (int i = 0; i < this->size(); ++i) {
        m_specresp[i] *= scale;
    }

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Divide Auxiliary Response File values
 *
 * @param[in] scale Division factor.
 * @return Divided Auxiliary Response File.
 *
 * Divides the values of the Auxiliary Response File by a division factor.
 ***************************************************************************/
GArf& GArf::operator/=(const double& scale)
{
    // Divide ARF values
    for (int i = 0; i < this->size(); ++i) {
        m_specresp[i] /= scale;
    }

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Return additional vector column
 *
 * @param[in] colname Vector column name.
 * @return Vector column.
 *
 * Returns reference to additional vector column.
 ***************************************************************************/
std::vector<double>& GArf::operator[](const std::string& colname)
{
    // Determine index of additional column (-1 if not found)
    int index = column_index(colname);

    // Throw exception if index not found
    if (index == -1) {
        std::string msg = "Could not find additional column with name \""+
                          colname+"\".";
        throw GException::invalid_value(G_OPERATOR, msg);
    }

    // Return reference to vector column
    return (m_coldata[index]);
}


/***********************************************************************//**
 * @brief Return additional vector column (const version)
 *
 * @param[in] colname Vector column name.
 * @return Vector column.
 *
 * Returns reference to additional vector column.
 ***************************************************************************/
const std::vector<double>& GArf::operator[](const std::string& colname) const
{
    // Determine index of additional column (-1 if not found)
    int index = column_index(colname);

    // Throw exception if index not found
    if (index == -1) {
        std::string msg = "Could not find additional column with name \""+
                          colname+"\".";
        throw GException::invalid_value(G_OPERATOR, msg);
    }

    // Return reference to vector column
    return (m_coldata[index]);
}


/***********************************************************************//**
 * @brief Return vector column content as function of energy
 *
 * @param[in] colname Vector column name.
 * @param[in] energy Energy.
 * @return Vector column content.
 *
 * Returns content of additional vector column as function of @p energy.
 ***************************************************************************/
double GArf::operator()(const std::string& colname,
                        const GEnergy&     energy) const
{
    // Determine index of additional column (-1 if not found)
    int index = column_index(colname);

    // Throw exception if index not found
    if (index == -1) {
        std::string msg = "Could not find additional column with name \""+
                          colname+"\".";
        throw GException::invalid_value(G_OPERATOR2, msg);
    }

    // Get reference to vector column
    const std::vector<double>& column = m_coldata[index];

    // Perform log-log interpolation of vector column content
    m_logetrue.set_value(energy.log10TeV());
    double wgt_left    = m_logetrue.wgt_left();
    double wgt_right   = m_logetrue.wgt_right();
    double value_left  = column[m_logetrue.inx_left()];
    double value_right = column[m_logetrue.inx_right()];
    double value       = 0.0;
    if (value_left > 0.0 && value_right > 0.0) {
        value = std::exp(wgt_left  * std::log(value_left) +
                         wgt_right * std::log(value_right));
    }
    if (value < 0.0) {
        value = 0.0;
    }

    // Return value
    return value;
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
    // Throw an exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT1, index, 0, size()-1);
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
    // Throw an exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT1, index, 0, size()-1);
    }

    // Return reference
    return (m_specresp[index]);
}


/***********************************************************************//**
 * @brief Return content of additional columns
 *
 * @param[in] index Bin index [0,...,size()-1].
 * @param[in] col Columns index [0,...,columns()-1].
 *
 * @exception GException::out_of_range
 *            Spectral bin or column index is out of range.
 *
 * Returns reference to content of additional columns.
 ***************************************************************************/
double& GArf::at(const int& index, const int& col)
{
    // Throw an exception if bin or column index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT2, "Bin index", index, size());
    }
    if (col < 0 || col >= columns()) {
        throw GException::out_of_range(G_AT2, "Column index", col, columns());
    }

    // Return reference
    return (m_coldata[col][index]);
}


/***********************************************************************//**
 * @brief Return content of additional columns (const version)
 *
 * @param[in] index Bin index [0,...,size()-1].
 * @param[in] col Columns index [0,...,columns()-1].
 *
 * @exception GException::out_of_range
 *            Spectral bin or column index is out of range.
 *
 * Returns reference to content of additional columns.
 ***************************************************************************/
const double& GArf::at(const int& index, const int& col) const
{
    // Throw an exception if bin or column index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT2, "Bin index", index, size());
    }
    if (col < 0 || col >= columns()) {
        throw GException::out_of_range(G_AT2, "Column index", col, columns());
    }

    // Return reference
    return (m_coldata[col][index]);
}


/***********************************************************************//**
 * @brief Append additional column to spectrum
 *
 * @param[in] name Additional column name.
 * @param[in] column Additional column data.
 ***************************************************************************/
void GArf::append(const std::string& name, const std::vector<double>& column)
{
    // Throw an exception if the number of elements in the column does not
    // correspond to the size of the spectrum
    if (column.size() != size()) {
        std::string msg = "Size of column "+gammalib::str(column.size())+
                          " is incompatible with size of spectrum "+
                          gammalib::str(size())+".";
        throw GException::invalid_argument(G_APPEND, msg);
    }

    // Append column name
    m_colnames.push_back(name);
    
    // Append column data
    m_coldata.push_back(column);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Auxiliary Response File
 *
 * @param[in] filename File name.
 *
 * Loads the Auxiliary Response File from the `SPECRESP` extension of the
 * FITS file.
 ***************************************************************************/
void GArf::load(const GFilename& filename)
{
    // Clear response
    clear();

    // Open FITS file (without extension name as the user is not allowed
    // to modify the extension names)
    GFits fits(filename.url());

    // Read ARF data
    read(fits);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename.url();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Auxiliary Response File
 *
 * @param[in] filename File name.
 * @param[in] clobber Overwrite existing file?
 *
 * Saves the Auxiliary Response File into a FITS file. If a file with the
 * given @p filename does not yet exist it will be created. If the file
 * exists it can be overwritten if the @p clobber flag is set to `true`.
 * Otherwise an exception is thrown.
 *
 * The method will save the `SPECRESP` binary FITS table into the FITS file
 * that contains the values of the Auxiliary Response File.
 ***************************************************************************/
void GArf::save(const GFilename& filename, const bool& clobber) const
{
    // Create FITS file
    GFits fits;

    // Write ARF into file
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
 * @brief Read Auxiliary Response File
 *
 * @param[in] fits FITS file.
 *
 * Loads the Auxiliary Response File from the `SPECRESP` extension of the
 * FITS file.
 ***************************************************************************/
void GArf::read(const GFits& fits)
{
    // Clear response
    clear();

    // Get ARF table
    const GFitsTable& table = *fits.table(gammalib::extname_arf);

    // Read ARF data
    read(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Auxiliary Response File
 *
 * @param[in] table ARF FITS table.
 *
 * Reads the Auxiliary Response File from a FITS table. The true energy
 * boundaries are expected in the `ENERG_LO` and `ENERG_HI` columns, the
 * response information is expected in the `SPECRESP` column.
 *
 * The method will analyze the unit of the `SPECRESP` column, and if either
 * `m**2`, `m^2` or `m2` are encountered, multiply the values of the column
 * by \f$10^4\f$ to convert the response into units of \f$cm^2\f$. Units of
 * the `ENERG_LO` and `ENERG_HI` columns are also interpreted for conversion.
 *
 * See
 * http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc4
 * for details about the Auxiliary Response File format.
 ***************************************************************************/
void GArf::read(const GFitsTable& table)
{
    // Clear members
    clear();

    // Get pointer to data columns
    const GFitsTableCol* energy_lo = table["ENERG_LO"];
    const GFitsTableCol* energy_hi = table["ENERG_HI"];
    const GFitsTableCol* specresp  = table["SPECRESP"];

    // Determine effective area conversion factor. Internal
    // units are cm^2
    std::string u_specresp =
         gammalib::tolower(gammalib::strip_whitespace(specresp->unit()));
    double      c_specresp = 1.0;
    if (u_specresp == "m**2" || u_specresp == "m^2" || u_specresp == "m2") {
        c_specresp = 10000.0;
    }

    // Extract number of energy bins
    int num = energy_lo->nrows();

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

    // Read any additional columns
    for (int icol = 0; icol < table.ncols(); ++icol) {

        // Fall through if the column is a standard column
        std::string colname(table[icol]->name());
        if ((colname == "ENERG_LO") ||
            (colname == "ENERG_HI") ||
            (colname == "SPECRESP")) {
            continue;
        }

        // Get pointer to column
        const GFitsTableCol* column = table[icol];

        // Set column vector
        std::vector<double> coldata;
        for (int i = 0; i < num; ++i) {
            coldata.push_back(column->real(i));
        }

        // Append column
        append(colname, coldata);

    } // endfor: looped over all additional columns

    // Set log true energy node array
    set_logetrue();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Auxiliary Response File
 *
 * @param[in] fits FITS file.
 *
 * Writes the Auxiliary Response File into the `SPECRESP` extension of the
 * FITS file. An existing extension with the same name will be removed from
 * the FITS file before writing.
 *
 * The method writes the boundaries of the true energy bins into the
 * `ENERG_LO` and `ENERG_HI` columns, response information will be written
 * into the `SPECRESP` column.
 *
 * See
 * http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc4
 * for details about the Auxiliary Response File format.
 ***************************************************************************/
void GArf::write(GFits& fits) const
{
    // If the FITS object contains already an extension with the same
    // name then remove now this extension
    if (fits.contains(gammalib::extname_arf)) {
        fits.remove(gammalib::extname_arf);
    }

    // Set column length
    int length = size();

    // Continue only if there are bins
    if (length > 0) {

        // Create new binary table
        GFitsBinTable hdu;

        // Allocate floating point vector columns
        GFitsTableFloatCol energy_lo("ENERG_LO", length);
        GFitsTableFloatCol energy_hi("ENERG_HI", length);
        GFitsTableFloatCol specresp("SPECRESP",  length);

        // Fill columns
        for (int i = 0; i < length; ++i) {
            energy_lo(i) = (float)m_ebounds.emin(i).keV();
            energy_hi(i) = (float)m_ebounds.emax(i).keV();
            specresp(i)  = (float)m_specresp[i];
        }

        // Set column units
        energy_lo.unit("keV");
        energy_hi.unit("keV");
        specresp.unit("cm**2");

        // Set table attributes
        hdu.extname(gammalib::extname_arf);

        // Append columns to table
        hdu.append(energy_lo);
        hdu.append(energy_hi);
        hdu.append(specresp);

        // Append any additional columns
        for (int icol = 0; icol < columns(); ++icol) {

            // Allocate floating point vector columns
            GFitsTableFloatCol column(m_colnames[icol], length);

            // Fill columns
            for (int i = 0; i < length; ++i) {
                column(i) = (float)m_coldata[icol][i];
            }

            // Append columns to table
            hdu.append(column);

        } // endfor: looped over all additional columns

        // Write mandatory header keywords
        hdu.card("TELESCOP", "unknown",  "Telescope");
        hdu.card("INSTRUME", "unknown",  "Instrument");
        hdu.card("FILTER",   "none",     "Filter");
        hdu.card("HDUCLASS", "OGIP",     "Format conforms to OGIP srandard");
        hdu.card("HDUCLAS1", "RESPONSE", "Extension contains response data");
        hdu.card("HDUCLAS2", "SPECRESP", "Extension contains an ARF");
        hdu.card("HDUVERS",  "1.1.0",    "Version of the file format");

        // Append HDU to FITS file
        fits.append(hdu);

    } // endif: there were data to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Auxiliary Response File
 *
 * @param[in] chatter Chattiness.
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
    m_logetrue.clear();
    m_specresp.clear();
    m_colnames.clear();
    m_coldata.clear();

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
    m_logetrue = arf.m_logetrue;
    m_specresp = arf.m_specresp;
    m_colnames = arf.m_colnames;
    m_coldata  = arf.m_coldata;

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


/***********************************************************************//**
 * @brief Set true energy node array
 ***************************************************************************/
void GArf::set_logetrue(void)
{
    // Clear node array
    m_logetrue.clear();

    // Continue only if there are true energies in Arf
    int netrue = m_ebounds.size();
    if (netrue > 0) {

        // Reserve space in node array
        m_logetrue.reserve(netrue);

        // Append all log mean energies to node array
        for (int i = 0; i < netrue; ++i) {

            // Get log mean of true energy in TeV
            double logE = m_ebounds.elogmean(i).log10TeV();

            // Append energy to node array
            m_logetrue.append(logE);

        } // endfor: appended log mean energies

    } // endif: there were true energies in Arf

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns index of additional vector column
 *
 * @param[in] colname Vector column name.
 * @return Vector column index.
 *
 * Returns index of additional vector column. Returns -1 if the vector column
 * has not found.
 ***************************************************************************/
int GArf::column_index(const std::string& colname) const
{
    // Initialise index
    int index(-1);

    // Search vector column name
    for (int i = 0; i < columns(); ++i) {
        if (m_colnames[i] == colname) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}
