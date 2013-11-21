/***************************************************************************
 *        GFitsTableCol.cpp - FITS table column abstract base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GFitsTableCol.cpp
 * @brief Abstract FITS table column class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ELEMENTS                            "GFitsTableCol::elements(int&)"
#define G_LOAD_COLUMN_FIXED              "GFitsTableCol::load_column_fixed()"
#define G_LOAD_COLUMN_VARIABLE        "GFitsTableCol::load_column_variable()"
#define G_SAVE_COLUMN_FIXED              "GFitsTableCol::save_column_fixed()"
#define G_SAVE_COLUMN_VARIABLE        "GFitsTableCol::save_column_variable()"
#define G_OFFSET                          "GFitsTableCol::offset(int&, int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_CALL_GRAPH                        //!< Dump call graph in console


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GFitsTableCol::GFitsTableCol(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Column constructor
 *
 * @param[in] name Name of column.
 * @param[in] length Length of column (number of rows).
 * @param[in] number Vector size of column.
 * @param[in] width Width of single column element.
 * @param[in] variable Variable-length flag (defaults to false).
 *
 * Construct column instance from name, length, vector size and column width.
 * The repeat value, required for binary tables, is calculated internally.
 * The optional parameter @p variable specifies whether the column is a
 * variable-length (true) or fixed-length (false) column. By default, a
 * fixed-length column will be allocated.
 ***************************************************************************/
GFitsTableCol::GFitsTableCol(const std::string& name,
                             const int&         length,
                             const int&         number,
                             const int&         width,
                             const bool&        variable)
{
    // Initialise class members for clean destruction
    init_members();

    // Store attributes
    m_name   = name;
    m_length = length;
    m_number = number;
    m_width  = width;

    // Calculate repeat value (only used for binary table!)
    m_repeat = m_number * m_width;

    // If column is a variable-length column then initialise row start
    // array and set variable-length flag
    if (variable) {
        m_variable = true;
        m_rowstart.assign(length+1, 0);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] column Column.
 ***************************************************************************/
GFitsTableCol::GFitsTableCol(const GFitsTableCol& column)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(column);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFitsTableCol::~GFitsTableCol(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] column Column.
 * @return Column.
 ***************************************************************************/
GFitsTableCol& GFitsTableCol::operator=(const GFitsTableCol& column)
{
    // Execute only if object is not identical
    if (this != &column) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(column);

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
 * @brief Set number of column elements for specific row
 *
 * @param[in] row Row index.
 * @param[in] elements Number of elements in @p row.
 *
 * Sets the number of elements in column for a specific @p row.
 *
 * @todo Implement method.
 ***************************************************************************/
void GFitsTableCol::elements(const int& row, const int& elements)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns number of elements in column for specific row
 *
 * @param[in] row Row index.
 * @return Number of elements in column at @p row.
 *
 * @exception GException::out_of_range
 *            Row index out of valid range.
 *
 * Returns the number of elements in the column for a specific @p row. For a
 * fixed-length column the returned number is independent of the row. For a
 * variable-length column the returned number is then length of the column
 * for the specified row.
 ***************************************************************************/
int GFitsTableCol::elements(const int& row) const
{
    // Check row value
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_length) {
        throw GException::out_of_range(G_ELEMENTS, row, 0, m_length-1);
    }
    #endif

    // Get number of elements
    int number = (m_variable) ? m_rowstart[row+1] - m_rowstart[row]
                              : m_number; 

    // Return number
    return number;
}


/***********************************************************************//**
 * @brief Returns TFORM code for binary table column
 *
 * Constructs the TFORM code for a binary table column, supporting
 * fixed-length and variable-length column types.
 ***************************************************************************/
std::string GFitsTableCol::tform_binary(void) const
{
    // Set type code string
    std::string typecode = "";
    switch (std::abs(type())) {
    case __TBIT:
        typecode = "X";
        break;
    case __TBYTE:
        typecode = "B";
        break;
    case __TLOGICAL:
        typecode = "L";
        break;
    case __TSTRING:
        // If there are substrings then add width of substring
        typecode = "A";
        if (repeat() > width()) {
            typecode.append(gammalib::str(width()));
        }
        break;
    case __TUSHORT:
        typecode = "U";
        break;
    case __TSHORT:
        typecode = "I";
        break;
    case __TULONG:
        typecode = "V";
        break;
    case __TLONG:
        typecode = "J";
        break;
    case __TFLOAT:
        typecode = "E";
        break;
    case __TLONGLONG:
        typecode = "K";
        break;
    case __TDOUBLE:
        typecode = "D";
        break;
    case __TCOMPLEX:
        typecode = "C";
        break;
    case __TDBLCOMPLEX:
        typecode = "M";
        break;
    default:
        break;
    }

    // Build TFORM code
    std::string tform;
    if (std::abs(type()) == __TSTRING) {
        tform.append(gammalib::str(repeat()));
    }
    else {
        tform.append(gammalib::str(number()));
    }
    if (isvariable()) {
        tform.append("P");
    }
    tform.append(typecode);
    if (isvariable() && m_varlen > 0) {
        tform.append("("+gammalib::str(m_varlen)+")");
    }
    
    // Return TFORM code
    return tform;
}


/***********************************************************************//**
 * @brief Print column information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing column information.
 *
 * @todo Format and cfitsio information is mainly for debugging. This could
 * vanish in a more stable version of the code, or it could be compiled-in
 * conditionally using a debug option. Alternatively, a higher chatter level
 * may be required to see this information.
 ***************************************************************************/
std::string GFitsTableCol::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append formatted column name. Optionally add units
        if (unit().length() > 0) {
            result.append(gammalib::parformat(name()+" ("+unit()+")"));
        }
        else {
            result.append(gammalib::parformat(name()));
        }

        // Append column number. This will be "-" if the column does not exist
        // in the FITS file.
        if (m_colnum > 0) {
            result.append(gammalib::right(gammalib::str(m_colnum),4)+" ");
        }
        else {
            result.append(gammalib::right("[-]",4)+" ");
        }

        // Append loading information
        if (isloaded()) {
            result.append("[loaded]     ");
        }
        else {
            result.append("[not loaded] ");
        }

        // Append format information
        result.append("["+tform_binary()+","+ascii_format()+"]");

        // Append dimensions (if available)
        if (!dim().empty()) {
    
            // Build TDIM string
            std::string value = "("+gammalib::str(dim()[0]);
            for (int k = 1; k < dim().size(); ++k) {
                value += ","+gammalib::str(dim()[k]);
            }
            value += ")";
        
            // Append
            result.append(" "+value);
        }

        // Append cfitsio information
        if (isvariable()) {
            result.append(" repeat=" + gammalib::str(repeat()));
            result.append(" width="  + gammalib::str(width()));
            result.append(" number=" + gammalib::str(number()));
            result.append(" length=" + gammalib::str(length()));
            result.append(" size="   + gammalib::str(m_size));
            result.append(" varlen=");
            if (m_varlen > 0) {
                result.append(gammalib::str(m_varlen));
            }
            else {
                result.append("undetermined");
            }
        }
        else {
            result.append(" repeat=" + gammalib::str(repeat()));
            result.append(" width="  + gammalib::str(width()));
            result.append(" number=" + gammalib::str(number()));
            result.append(" length=" + gammalib::str(length()));
            result.append(" size="   + gammalib::str(m_size));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * Refer to GFitsTableCol::save_column() for more information.
 ***************************************************************************/
void GFitsTableCol::save(void)
{
    // Save column
    save_column();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load table column from FITS file
 *
 * Loads table column from FITS file by calling the load_column_variable()
 * method for variable-length columns and load_column_fixed() for fixed-
 * length columns.
 ***************************************************************************/
void GFitsTableCol::load_column(void)
{
    // Load variable-length or fixed-length column from FITS file
    if (isvariable()) {
        load_column_variable();
    }
    else {
        load_column_fixed();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load fixed-length column from FITS file
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 * @exception GException::fits_error
 *            An error occured while loading column data from FITS file.
 *
 * If a FITS file is attached to the column the data are loaded into memory
 * from the FITS file. If no FITS file is attached, memory is allocated
 * to hold the column data and all cells are set to 0.
 *
 * The method makes use of the virtual methods 
 * GFitsTableCol::alloc_data,
 * GFitsTableCol::init_data,
 * GFitsTableCol::ptr_data, and
 * GFitsTableCol::ptr_nulval.
 * These methods are implemented by the derived column classes which 
 * implement a specific storage class (i.e. float, double, short, ...).
 ***************************************************************************/
void GFitsTableCol::load_column_fixed(void)
{
    // Calculate size of memory
    m_size = m_number * m_length;

    // Load only if the column has a positive size
    if (m_size > 0) {

        // Allocate and initialise fresh memory
        alloc_data();
        init_data();

        // If a FITS file is attached then try loading column data from the
        // FITS file. This may fail in case that no data has yet been written
        // to the FITS file. In that case we just skip loading and return
        // the initalised column ... 
        if (FPTR(m_fitsfile)->Fptr != NULL) {

            // Move to the HDU
            int status = 0;
            status     = __ffmahd(FPTR(m_fitsfile),
                                  (FPTR(m_fitsfile)->HDUposition)+1,
                                  NULL, &status);
            
            // If this failed because:
            // - the primary HDU was not found (status 252)
            // - we moved past the file (status 107)
            // we assume that no data have yet been written to the file and
            // we skip the loading.
            if (status != 252 && status != 107) {
            
                // Break on any other cfitsio error
                if (status != 0) {
                    throw GException::fits_hdu_not_found(G_LOAD_COLUMN_FIXED,
                                  (FPTR(m_fitsfile)->HDUposition)+1,
                                  status);
                }

                // Load data
                status = __ffgcv(FPTR(m_fitsfile), m_type, m_colnum, 
                                 1, 1, m_size, ptr_nulval(), ptr_data(),
                                 &m_anynul, &status);
                if (status != 0) {
                    throw GException::fits_error(G_LOAD_COLUMN_FIXED, status,
                                    "for column '"+m_name+"'.");
                }
        
            } // endif: no primary HDU found

        } // endif: there was a FITS file attached

    } // endif: column has a positive size

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load variable-length column from FITS file
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 * @exception GException::fits_error
 *            An error occured while loading column data from FITS file.
 *
 * If a FITS file is attached to the column the data are loaded into memory
 * from the FITS file. If no FITS file is attached, memory is allocated
 * to hold the column data and all cells are set to 0.
 *
 * The method makes use of the virtual methods 
 * GFitsTableCol::alloc_data,
 * GFitsTableCol::init_data,
 * GFitsTableCol::ptr_data, and
 * GFitsTableCol::ptr_nulval.
 * These methods are implemented by the derived column classes which 
 * implement a specific storage class (i.e. float, double, short, ...).
 ***************************************************************************/
void GFitsTableCol::load_column_variable(void)
{
    // If a FITS file is attached then try loading column data from the
    // FITS file. This may fail in case that no data has yet been written
    // to the FITS file. In that case we just skip loading and return
    // the initalised column ... 
    if (FPTR(m_fitsfile)->Fptr != NULL) {

        // Move to the HDU
        int status = 0;
        status     = __ffmahd(FPTR(m_fitsfile),
                              (FPTR(m_fitsfile)->HDUposition)+1,
                              NULL,
                              &status);
            
        // If this failed because:
        // - the primary HDU was not found (status 252)
        // - we moved past the file (status 107)
        // we assume that no data have yet been written to the file and
        // we skip the loading.
        if (status != 252 && status != 107) {
            
            // Break on any other cfitsio error
            if (status != 0) {
                throw GException::fits_hdu_not_found(G_LOAD_COLUMN_VARIABLE,
                                  (FPTR(m_fitsfile)->HDUposition)+1,
                                  status);
            }

            // Allocate rowstart array
            m_rowstart.assign(m_length+1, 0);

            // Determine the column length for each row by looping over
            // all rows and derive the total memory requirement
            m_size        = 0;
            m_varlen      = 0;
            m_rowstart[0] = 0;
            for (int row = 0; row < m_length; ++row) {

                // Initialise offset and repeat
                long offset(0);
                long repeat(0);

                // Get variable-length of row in repeat
                status = __ffgdes(FPTR(m_fitsfile),
                                  m_colnum,
                                  row+1,
                                  &repeat,
                                  &offset,
                                  &status);
                if (status != 0) {
                    std::string msg = "Unable to get descriptor of row "+
                                      gammalib::str(row+1)+" of column '"+
                                      name()+"' from FITS file.";
                    throw GException::fits_error(G_LOAD_COLUMN_VARIABLE, status,
                                                 msg);
                }

                // Store start of next row
                m_rowstart[row+1] = m_rowstart[row] + repeat;
                m_size           += repeat;
                if (repeat > m_varlen) {
                    m_varlen = repeat;
                }

            } // endfor: looped over all rows

            // Allocate and initialise fresh memory
            alloc_data();
            init_data();

            // Load data for each row
            for (int row = 0; row < m_length; ++row) {

                // Initialise anynul
                int anynul(0);

                // Load data
                status = __ffgcv(FPTR(m_fitsfile),
                                 std::abs(m_type),
                                 m_colnum, 
                                 row+1,
                                 1,
                                 elements(row),
                                 ptr_nulval(),
                                 ptr_data(m_rowstart[row]),
                                 &anynul,
                                 &status);
                if (status != 0) {
                    std::string msg = "Unable to load row "+gammalib::str(row+1)+""
                                      " of column '"+name()+"' from FITS file.";
                    throw GException::fits_error(G_LOAD_COLUMN_VARIABLE, status,
                                                 msg);
                }

                // Increment anynul
                m_anynul += anynul;

            } // endfor: looped over all rows

        } // endif: no primary HDU found

    } // endif: there was a FITS file attached

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * Save table column into FITS file by calling the save_column_variable()
 * method for variable-length columns and save_column_fixed() for fixed-
 * length columns.
 ***************************************************************************/
void GFitsTableCol::save_column(void)
{
    // Save variable-length or fixed-length column into FITS file
    if (isvariable()) {
        save_column_variable();
    }
    else {
        save_column_fixed();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 * @exception GException::fits_error
 *            Error occured during writing of the column data.
 *
 * The table column is only saved if it is linked to a FITS file and if the
 * data are indeed present in the class instance. This avoids saving of data
 * that have not been modified.
 *
 * The method make use of the virtual methods 
 *   GFitsTableCol::ptr_data and
 *   GFitsTableCol::ptr_nulval.
 * These methods are implemented by the derived column classes which 
 * implement a specific storage class (i.e. float, double, short, ...).
 ***************************************************************************/
void GFitsTableCol::save_column_fixed(void)
{
    // Continue only if a FITS file is connected and data have been loaded
    if (FPTR(m_fitsfile)->Fptr != NULL && m_colnum > 0 && ptr_data() != NULL) {

        // Move to the HDU
        int status = 0;
        status     = __ffmahd(FPTR(m_fitsfile),
                              (FPTR(m_fitsfile)->HDUposition)+1, NULL,
                              &status);
        if (status != 0) {
            throw GException::fits_hdu_not_found(G_SAVE_COLUMN_FIXED,
                              (FPTR(m_fitsfile)->HDUposition)+1,
                              status);
        }

        // Save the column data
        status = __ffpcn(FPTR(m_fitsfile), m_type, m_colnum, 1, 1,
                         m_size, ptr_data(), ptr_nulval(), &status);
        if (status != 0) {
            std::string msg = "Unable to save column '"+name()+"' to"
                              " FITS file.";
            throw GException::fits_error(G_SAVE_COLUMN_FIXED, status, msg);
        }

    } // endif: FITS file was connected

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 * @exception GException::fits_error
 *            Error occured during writing of the column data.
 *
 * The table column is only saved if it is linked to a FITS file and if the
 * data are indeed present in the class instance. This avoids saving of data
 * that have not been modified.
 *
 * The method make use of the virtual methods 
 *   GFitsTableCol::ptr_data and
 *   GFitsTableCol::ptr_nulval.
 * These methods are implemented by the derived column classes which 
 * implement a specific storage class (i.e. float, double, short, ...).
 ***************************************************************************/
void GFitsTableCol::save_column_variable(void)
{
    // Continue only if a FITS file is connected and data have been loaded
    if (FPTR(m_fitsfile)->Fptr != NULL && m_colnum > 0 && ptr_data() != NULL) {

        // Move to the HDU
        int status = 0;
        status     = __ffmahd(FPTR(m_fitsfile),
                              (FPTR(m_fitsfile)->HDUposition)+1, NULL,
                              &status);
        if (status != 0) {
            throw GException::fits_hdu_not_found(G_SAVE_COLUMN_VARIABLE,
                              (FPTR(m_fitsfile)->HDUposition)+1,
                              status);
        }

        // Save the column data row-by-row
        for (int row = 0; row < m_length; ++row) {

            // Save row data
            status = __ffpcn(FPTR(m_fitsfile),
                             std::abs(m_type),
                             m_colnum, 
                             row+1,
                             1,
                             elements(row),
                             ptr_data(m_rowstart[row]),
                             ptr_nulval(),
                             &status);
            if (status != 0) {
                std::string msg = "Unable to save row "+gammalib::str(row+1)+""
                                  " of column '"+name()+"' to FITS file.";
                throw GException::fits_error(G_SAVE_COLUMN_VARIABLE, status,
                                             msg);
            }

        } // endfor: looped over rows

    } // endif: FITS file was connected

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute offset of column element in memory
 *
 * @param[in] row Row of column.
 * @param[in] inx Vector index in column row.
 *
 * @exception GException::out_of_range
 *            Table row or vector index are out of valid range.
 *
 * Computes the offset of a column element in the storage array from the
 * @p row number and the vector index. The method also supports both
 * variable-length and fixed-length columns.
 ***************************************************************************/
int GFitsTableCol::offset(const int& row, const int& inx) const
{
    // Check row value
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_length) {
        throw GException::out_of_range(G_OFFSET, row, 0, m_length-1);
    }
    #endif

    // Check inx value
    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= elements(row)) {
        throw GException::out_of_range(G_OFFSET, inx, 0, elements(row));
    }
    #endif

    // Calculate pixel offset
    int offset = (isvariable()) ? m_rowstart[row] + inx : row * m_number + inx;

    // Return offset
    return offset;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsTableCol::init_members(void)
{
    // Optionally print call graph
    #if defined(G_CALL_GRAPH)
    printf("GFitsTableCol::init_members\n");
    #endif

    // Allocate FITS file pointer
    m_fitsfile = new __fitsfile;
    FPTR(m_fitsfile)->HDUposition = 0;
    FPTR(m_fitsfile)->Fptr        = NULL;

    // Initialise members
    m_name.clear();
    m_unit.clear();
    m_dim.clear();
    m_rowstart.clear();
    m_colnum   = 0;
    m_type     = 0;
    m_repeat   = 0;
    m_width    = 0;
    m_number   = 0;
    m_length   = 0;
    m_variable = false;
    m_varlen   = 0;
    m_size     = 0;
    m_anynul   = 0;

    // Optionally print call graph
    #if defined(G_CALL_GRAPH)
    printf("exit GFitsTableCol::init_members\n");
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] column Column to be copied.
 ***************************************************************************/
void GFitsTableCol::copy_members(const GFitsTableCol& column)
{
    // Copy attributes
    m_name     = column.m_name;
    m_unit     = column.m_unit;
    m_dim      = column.m_dim;
    m_colnum   = column.m_colnum;
    m_type     = column.m_type;
    m_repeat   = column.m_repeat;
    m_width    = column.m_width;
    m_number   = column.m_number;
    m_length   = column.m_length;
    m_variable = column.m_variable;
    m_varlen   = column.m_varlen;
    m_rowstart = column.m_rowstart;
    m_size     = column.m_size;
    m_anynul   = column.m_anynul;
    FPTR_COPY(m_fitsfile, column.m_fitsfile);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableCol::free_members(void)
{
    // Optionally print call graph
    #if defined(G_CALL_GRAPH)
    printf("GFitsTableCol::free_members\n");
    #endif

    // Free memory
    if (m_fitsfile != NULL) delete FPTR(m_fitsfile);

    // Mark memory as free
    m_fitsfile = NULL;

    // Optionally print call graph
    #if defined(G_CALL_GRAPH)
    printf("exit GFitsTableCol::free_members\n");
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Connect table column to FITS file
 *
 * @param[in] vptr Column file void pointer.
 *
 * A table column is connected to a FITS file when the m_fitsfile holds a
 * FITS file pointer. In that way, the column knows to which file it belongs.
 ***************************************************************************/
void GFitsTableCol::connect(void* vptr)
{
    // Connect table column by copying the column file pointer
    FPTR_COPY(m_fitsfile, vptr);

    // Return
    return;
}
