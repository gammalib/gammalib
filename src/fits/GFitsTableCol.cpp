/***************************************************************************
 *        GFitsTableCol.cpp  - FITS table column abstract base class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsTableCol.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD_COLUMN                          "GFitsTableCol::load_column()"
#define G_SAVE_COLUMN                          "GFitsTableCol::save_column()"
#define G_OFFSET                           "GFitsTableCol::offset(int&,int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsTableCol::GFitsTableCol(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] name Name of column.
 * @param[in] length Length of column.
 * @param[in] number Vector size of column.
 * @param[in] width Width of single column element.
 *
 * Construct column instance from name, length, vector size and column width.
 * The repeat value, required for binary tables, is calculated internally.
 ***************************************************************************/
GFitsTableCol::GFitsTableCol(const std::string& name,
                             const int&         length,
                             const int&         number,
                             const int&         width)
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

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] column Column from which the class instance should be built
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
 * @param[in] column Column which will be assigned
 ***************************************************************************/
GFitsTableCol& GFitsTableCol::operator= (const GFitsTableCol& column)
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
 * @brief Set column name
 *
 * @param[in] name Name of the column.
 ***************************************************************************/
void GFitsTableCol::name(const std::string& name)
{
    // Set name
    m_name = name;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set column unit
 *
 * @param[in] unit Unit of the column.
 ***************************************************************************/
void GFitsTableCol::unit(const std::string& unit)
{
    // Set name
    m_unit = unit;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns column name
 ***************************************************************************/
std::string GFitsTableCol::name(void) const
{
    // Return Name
    return m_name;
}


/***********************************************************************//**
 * @brief Returns column unit
 ***************************************************************************/
std::string GFitsTableCol::unit(void) const
{
    // Return column length
    return m_unit;
}


/***********************************************************************//**
 * @brief Returns number of column in FITS file (starting from 1)
 ***************************************************************************/
int GFitsTableCol::colnum(void) const
{
    // Return column number
    return m_colnum;
}


/***********************************************************************//**
 * @brief Returns CFITSIO column type
 *
 * Returns one of the following:
 *   1 (TBIT)
 *  11 (TBYTE)
 *  12 (TSBYTE)
 *  14 (TLOGICAL)
 *  16 (TSTRING)
 *  21 (TUSHORT)
 *  21 (TSHORT)
 *  31 (TUINT)
 *  31 (TINT)
 *  41 (TULONG)
 *  41 (TLONG)
 *  42 (TFLOAT)
 *  81 (TLONGLONG)
 *  82 (TDOUBLE)
 *  83 (TCOMPLEX)
 * 163 (TDBLCOMPLEX)
 ***************************************************************************/
int GFitsTableCol::type(void) const
{
    // Return column type
    return m_type;
}


/***********************************************************************//**
 * @brief Returns column repeat value (only used for binary tables)
 ***************************************************************************/
int GFitsTableCol::repeat(void) const
{
    // Return column repeat value
    return m_repeat;
}


/***********************************************************************//**
 * @brief Returns width of one element in column
 ***************************************************************************/
int GFitsTableCol::width(void) const
{
    // Return width of one element in column
    return m_width;
}


/***********************************************************************//**
 * @brief Returns number of elements in a column
 ***************************************************************************/
int GFitsTableCol::number(void) const
{
    // Return number of elements in a column
    return m_number;
}


/***********************************************************************//**
 * @brief Returns column length
 ***************************************************************************/
int GFitsTableCol::length(void) const
{
    // Return column length
    return m_length;
}


/***********************************************************************//**
 * @brief Returns number of NULLs encountered
 ***************************************************************************/
int GFitsTableCol::anynul(void) const
{
    // Return column length
    return m_anynul;
}


/***********************************************************************//**
 * @brief Print column information
 *
 * @todo Format and cfitsio information is mainly for debugging. This could
 * be vanish in a more stable version of the code, or it could be compiled
 * in conditionally using a debug option.
 ***************************************************************************/
std::string GFitsTableCol::print(void) const
{
    // Initialise result string
    std::string result;

    // Append formatted column name
    result.append(parformat(m_name));

    // Append column number
    if (m_colnum > 0)
        result.append(right(str(m_colnum),4)+" ");
    else
        result.append(right("[-]",4)+" ");

    // Append loading information
    if (m_size == 0)
        result.append("[not loaded] ");
    else
        result.append("[loaded]     ");

    // Append format information
    result.append("["+binary_format()+","+ascii_format()+"]");

    // Append cfitsio information
    result.append(" repeat="+str(m_repeat));
    result.append(" width="+str(m_width));
    result.append(" number="+str(m_number));
    result.append(" length="+str(m_length));
    result.append(" size="+str(m_size));

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
 * @brief Fetch column data
 *
 * This method fetches column data when needed. It is declared const, so
 * that const data access methods can be implemented.
 *
 * If a FITS file is attached to the column the data are loaded into memory
 * from the FITS file. If no FITS file is attached, memory is allocated
 * to hold the column data and all cells are initialised.
 *
 * This method calls GFitsTableCol::load_column to do the job.
 ***************************************************************************/
void GFitsTableCol::fetch_data(void) const
{
    // Save column (circumvent const correctness)
    const_cast<GFitsTableCol*>(this)->load_column();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load table column from FITS file
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
void GFitsTableCol::load_column(void)
{
    // Calculate size of memory
    m_size = m_number * m_length;

    // Load only if the column has a positive size
    if (m_size > 0) {

        // Allocate and initialise fresh memory
        alloc_data();
        init_data();

        // If a FITS file is attached then load column data from the FITS
        // file
        if (FPTR(m_fitsfile)->Fptr != NULL) {

            // Move to the HDU
            int status = 0;
            status     = __ffmahd(FPTR(m_fitsfile),
                                  (FPTR(m_fitsfile)->HDUposition)+1,
                                  NULL, &status);
            if (status != 0)
                throw GException::fits_hdu_not_found(G_LOAD_COLUMN,
                                  (FPTR(m_fitsfile)->HDUposition)+1,
                                  status);

            // Load data
            status = __ffgcv(FPTR(m_fitsfile), m_type, m_colnum, 1, 1, m_size,
                             ptr_nulval(), ptr_data(), &m_anynul, &status);
            if (status != 0)
                throw GException::fits_error(G_LOAD_COLUMN, status,
                                  "for column \""+m_name+"\".");
        }

    } // endif: column has a positive size

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
void GFitsTableCol::save_column(void)
{
    // Continue only if a FITS file is connected and data have been loaded
    if (FPTR(m_fitsfile)->Fptr != NULL && m_colnum > 0 && ptr_data() != NULL) {

        // Move to the HDU
        int status = 0;
        status     = __ffmahd(FPTR(m_fitsfile),
                              (FPTR(m_fitsfile)->HDUposition)+1, NULL,
                              &status);
        if (status != 0)
            throw GException::fits_hdu_not_found(G_SAVE_COLUMN,
                              (FPTR(m_fitsfile)->HDUposition)+1,
                              status);

        // Save the column data
        status = __ffpcn(FPTR(m_fitsfile), m_type, m_colnum, 1, 1,
                         m_size, ptr_data(), ptr_nulval(), &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE_COLUMN, status);

    } // endif: FITS file was connected

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write column in output stream
 *
 * @param[in] os Output stream.
 ***************************************************************************/
std::ostream& GFitsTableCol::dump_column(std::ostream& os) const
{
    // Write column in output stream
    os << print();

    // Return stream
    return os;
}


/***********************************************************************//**
 * @brief Write column in logger
 *
 * @param[in] log Logger.
 ***************************************************************************/
GLog& GFitsTableCol::dump_column(GLog& log) const
{
    // Write column in logger
    log << print();

    // Return logger
    return log;
}


/***********************************************************************//**
 * @brief Convert row and vector index into column offset
 *
 * @param[in] row Row of column.
 * @param[in] inx Vector index in column row.
 *
 * @exception GException::out_of_range
 *            Table row or vector index are out of valid range.
 *
 * Converts the row and vector index of a column into a linear offset from
 * the column start.
 ***************************************************************************/
int GFitsTableCol::offset(const int& row, const int& inx) const
{
    // Check row value
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_OFFSET, row, 0, m_length-1);
    #endif

    // Check inx value
    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= m_number)
        throw GException::out_of_range(G_OFFSET, inx, 0, m_number-1);
    #endif

    // Calculate pixel offset
    int offset = row * m_number + inx;

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
    // Allocate FITS file pointer
    m_fitsfile = new __fitsfile;
    FPTR(m_fitsfile)->HDUposition = 0;
    FPTR(m_fitsfile)->Fptr        = NULL;

    // Initialise members
    m_name.clear();
    m_unit.clear();
    m_colnum = 0;
    m_type   = 0;
    m_repeat = 0;
    m_width  = 0;
    m_number = 0;
    m_length = 0;
    m_size   = 0;
    m_anynul = 0;

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
    m_colnum   = column.m_colnum;
    m_type     = column.m_type;
    m_repeat   = column.m_repeat;
    m_width    = column.m_width;
    m_number   = column.m_number;
    m_length   = column.m_length;
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
    // Free memory
    if (m_fitsfile != NULL) delete FPTR(m_fitsfile);

    // Mark memory as free
    m_fitsfile = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Connect table column to FITS file
 *
 * @param[in] vptr Column file void pointer.
 ***************************************************************************/
void GFitsTableCol::connect(void* vptr)
{
    // Connect table column by copying the column file pointer
    FPTR_COPY(m_fitsfile, vptr);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] column FITS table column.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GFitsTableCol& column)
{
    // Return output stream
    return (column.dump_column(os));
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] column FITS table column.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GFitsTableCol& column)
{
    // Return logger
    return (column.dump_column(log));
}
