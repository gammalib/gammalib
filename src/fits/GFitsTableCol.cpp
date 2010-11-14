/***************************************************************************
 *        GFitsTableCol.cpp  - FITS table column abstract base class       *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GFitsTableCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD_COLUMN                          "GFitsTableCol::load_column()"
#define G_SAVE_COLUMN                          "GFitsTableCol::save_column()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                   GFitsTableCol constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsTableCol::GFitsTableCol()
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
GFitsTableCol::~GFitsTableCol()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GFitsTableCol operators                         =
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
 =                        GFitsTableCol public methods                     =
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
 * @brief Get column name
 ***************************************************************************/
std::string GFitsTableCol::name(void)
{
    // Return Name
    return m_name;
}


/***********************************************************************//**
 * @brief Get number of column in FITS file (starting from 1)
 ***************************************************************************/
int GFitsTableCol::colnum(void)
{
    // Return column number
    return m_colnum;
}


/***********************************************************************//**
 * @brief Get CFITSIO column type
 *
 * Returns one of the following:
 *   1 (TBIT)
 *  11 (TBYTE)
 *  14 (TLOGICAL)
 *  16 (TSTRING)
 *  21 (TSHORT)
 *  31 (TINT)
 *  41 (TLONG)
 *  42 (TFLOAT)
 *  81 (TLONGLONG)
 *  82 (TDOUBLE)
 *  83 (TCOMPLEX)
 * 163 (TDBLCOMPLEX)
 ***************************************************************************/
int GFitsTableCol::type(void)
{
    // Return column type
    return m_type;
}


/***********************************************************************//**
 * @brief Get column repeat value (only used for binary tables)
 ***************************************************************************/
int GFitsTableCol::repeat(void)
{
    // Return column repeat value
    return m_repeat;
}


/***********************************************************************//**
 * @brief Get width of one element in column
 ***************************************************************************/
int GFitsTableCol::width(void)
{
    // Return width of one element in column
    return m_width;
}


/***********************************************************************//**
 * @brief Get number of elements in a column
 ***************************************************************************/
int GFitsTableCol::number(void)
{
    // Return number of elements in a column
    return m_number;
}


/***********************************************************************//**
 * @brief Get column length
 ***************************************************************************/
int GFitsTableCol::length(void)
{
    // Return column length
    return m_length;
}


/*==========================================================================
 =                                                                         =
 =                        GFitsTableCol private methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsTableCol::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_unit.clear();
    m_colnum               = 0;
    m_type                 = 0;
    m_repeat               = 0;
    m_width                = 0;
    m_number               = 0;
    m_length               = 0;
    m_size                 = 0;
    m_anynul               = 0;
    m_fitsfile.HDUposition = 0;
    m_fitsfile.Fptr        = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] column Column to be copied
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
    m_fitsfile = column.m_fitsfile;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableCol::free_members(void)
{
    // Free memory

    // Mark memory as freed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Connect table column to FITS file
 *
 * @param[in] fptr FITS file pointer to which the table column should be
 *                 connected
 ***************************************************************************/
void GFitsTableCol::connect(__fitsfile* fptr)
{
    // Connect Image
    m_fitsfile = *fptr;

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
 * The method make use of the virtual methods 
 *   GFitsTableCol::alloc_data,
 *   GFitsTableCol::init_data,
 *   GFitsTableCol::ptr_data, and
 *   GFitsTableCol::ptr_nulval.
 * These methods are implemented by the derived column classes which 
 * implement a specific storage class (i.e. float, double, short, ...).
 ***************************************************************************/
void GFitsTableCol::load_column(void)
{
    // Calculate size of memory
    m_size = m_number * m_length;

    // Load only if the column has a positive size
    if (m_size > 0) {

        // Allocate fresh memory
        alloc_data();

        // If a FITS file is attached then load column data from the FITS
        // file
        if (m_fitsfile.Fptr != NULL) {

            // Move to the HDU
            int status = 0;
            status     = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1,
                                  NULL, &status);
            if (status != 0)
                throw GException::fits_hdu_not_found(G_LOAD_COLUMN,
                                                 (m_fitsfile.HDUposition)+1,
                                                 status);

            // Load data
            status = __ffgcv(&m_fitsfile, m_type, m_colnum, 1, 1, m_size,
                             ptr_nulval(), ptr_data(), &m_anynul, &status);
            if (status != 0)
                throw GException::fits_error(G_LOAD_COLUMN, status);
        }

        // ... otherwise initialise all column values
        else
            init_data();

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
    if (m_fitsfile.Fptr != NULL && m_colnum > 0 && ptr_data() != NULL) {

        // Move to the HDU
        int status = 0;
        status     = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL,
                              &status);
        if (status != 0)
            throw GException::fits_hdu_not_found(G_SAVE_COLUMN,
                                                 (m_fitsfile.HDUposition)+1,
                                                 status);
        // Save the column data
        status = __ffpcn(&m_fitsfile, m_type, m_colnum, 1, 1,
                         (long)m_size, ptr_data(), ptr_nulval(), &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE_COLUMN, status);

    } // endif: FITS file was connected

    // Return
    return;
}


/***********************************************************************//**
 * @brief Put column information in output stream
 *
 * @param[in] os Output stream
 * @param[in] column Column to put in output stream
 ***************************************************************************/
void GFitsTableCol::dump_column(std::ostream& os, void* data) const
{
    // Put column name in stream
    os << "'" << m_name << "'";

    // Put FITS column number in stream
    if (m_colnum > 0)
        os << " [fits_colnum=" << m_colnum << "]";
    else
        os << " [not linked to FITS file]";

    // Put column type in stream
    os << " " << ascii_format();
    os << " " << binary_format();

    // Put data loading in stream
    if (data == NULL)
        os << " (not loaded)";
    else
        os << " (loaded in memory)";

    // Set data area size
    os << " size=" << m_size;

    // Set vector length
    os << " repeat=" << m_repeat;

    // Set width
    os <<  " width=" << m_width;

    // Set number
    os <<  " number=" << m_number;

    // Set length
    os << " length=" << m_length;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GFitsTableCol friends                         =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                    Other functions used by GFitsTableCol                =
 =                                                                         =
 ==========================================================================*/
