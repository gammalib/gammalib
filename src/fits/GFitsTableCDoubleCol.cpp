/***************************************************************************
 *  GFitsTableCDoubleCol.cpp  - FITS table double precision complex column *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsTableCDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */

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
GFitsTableCDoubleCol::GFitsTableCDoubleCol(void) : GFitsTableCol()
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
 * @param[in] size Vector size of column.
 ***************************************************************************/
GFitsTableCDoubleCol::GFitsTableCDoubleCol(const std::string& name,
                                           const int& length, const int& size) :
                      GFitsTableCol(name, length, size, 4)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] column Column from which class instance should be built.
 ***************************************************************************/
GFitsTableCDoubleCol::GFitsTableCDoubleCol(const GFitsTableCDoubleCol& column) :
                      GFitsTableCol(column)
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
GFitsTableCDoubleCol::~GFitsTableCDoubleCol(void)
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
 * @param[in] column Column which should be assigned
 ***************************************************************************/
GFitsTableCDoubleCol& GFitsTableCDoubleCol::operator= (const GFitsTableCDoubleCol& column)
{
    // Execute only if object is not identical
    if (this != &column) {

        // Copy base class members
        this->GFitsTableCol::operator=(column);

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


/***********************************************************************//**
 * @brief Column data access operator
 *
 * @param[in] row Row of column to access.
 * @param[in] inx Vector index in column row to access
 *
 * Provides access to data in a column.
 ***************************************************************************/
GFits::cdouble& GFitsTableCDoubleCol::operator() (const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Return data bin
    return m_data[offset(row, inx)];
}


/***********************************************************************//**
 * @brief Column data access operator (const variant)
 *
 * @param[in] row Row of column to access.
 * @param[in] inx Vector index in column row to access
 *
 * Provides access to data in a column.
 ***************************************************************************/
const GFits::cdouble& GFitsTableCDoubleCol::operator() (const int& row,
                                                        const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) ((GFitsTableCDoubleCol*)this)->fetch_data();

    // Return data bin
    return m_data[offset(row, inx)];
}


/*==========================================================================
 =                                                                         =
 =                               Public methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as string.
 ***************************************************************************/
std::string GFitsTableCDoubleCol::string(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Convert complex into string
    std::string result = str(m_data[offset(row,inx)].re) + ", " +
                         str(m_data[offset(row,inx)].im);

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Get double precision value (real part)
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as double precision.
 ***************************************************************************/
double GFitsTableCDoubleCol::real(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Convert complex into double
    double value = (double)m_data[offset(row,inx)].re;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Get integer value (real part)
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as integer.
 ***************************************************************************/
int GFitsTableCDoubleCol::integer(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Convert complex into int
    int value = (int)m_data[offset(row,inx)].re;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Set nul value
 *
 * @param[in] value Nul value.
 *
 * @todo To correctly reflect the nul value in the data, the column should
 * be reloaded. However, the column may have been changed, so in principle
 * saving is needed. However, we may not want to store the data, hence saving
 * is also not desired. We thus have to develop a method to update the
 * column information for a new nul value in place ...
 ***************************************************************************/
void GFitsTableCDoubleCol::nulval(const GFits::cdouble* value)
{
    // Allocate nul value
    alloc_nulval(value);

    // Update column
//    if (m_data != NULL) {
//        save();
//        load();
//    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsTableCDoubleCol::init_members(void)
{
    // Initialise members
    m_type   = __TDBLCOMPLEX;
    m_data   = NULL;
    m_nulval = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] column Column for which members should be copied.
 *
 * Sets the content of the vector column by copying from another column.
 * If the code is compiled with the small memory option, and if the source
 * column has not yet been loaded, then we only load the column temporarily
 * for copying purposes and release it again once copying is finished.
 ***************************************************************************/
void GFitsTableCDoubleCol::copy_members(const GFitsTableCDoubleCol& column)
{
    // Fetch column data if not yet fetched. The casting circumvents the
    // const correctness
    bool not_loaded = (column.m_data == NULL);
    if (not_loaded) ((GFitsTableCDoubleCol*)(&column))->fetch_data();

    // Copy attributes
    m_type = column.m_type;
    m_size = column.m_size;

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        alloc_data();
        for (int i = 0; i < m_size; ++i)
            m_data[i] = column.m_data[i];
    }

    // Copy NULL value
    alloc_nulval(column.m_nulval);

    // Small memory option: release column if it was fetch above
    #if defined(G_SMALL_MEMORY)
    if (not_loaded) ((GFitsTableCDoubleCol*)(&column))->release_data();
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableCDoubleCol::free_members(void)
{
    // Free memory
    if (m_data   != NULL) delete [] m_data;
    if (m_nulval != NULL) delete m_nulval;

    // Mark memory as freed
    m_data   = NULL;
    m_nulval = NULL;

    // Reset load flag
    m_size = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone column
 ***************************************************************************/
GFitsTableCDoubleCol* GFitsTableCDoubleCol::clone(void) const
{
    return new GFitsTableCDoubleCol(*this);
}


/***********************************************************************//**
 * @brief Returns format string of ASCII table
 ***************************************************************************/
std::string GFitsTableCDoubleCol::ascii_format(void) const
{
    // Initialize format string
    std::string format;

    // Set type code
    format.append("M");

    // Set width
    format.append(str(m_width));

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Returns format string of binary table
 ***************************************************************************/
std::string GFitsTableCDoubleCol::binary_format(void) const
{
    // Initialize format string
    std::string format;

    // Set number of elements
    format.append(str(m_number));

    // Set type code
    format.append("M");

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Allocates column data
 ***************************************************************************/
void GFitsTableCDoubleCol::alloc_data(void)
{
    // Free any existing memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free
    m_data = NULL;

    // Allocate new data
    if (m_size > 0)
        m_data = new GFits::cdouble[m_size];

    // Return
    return;
}


/***********************************************************************//**
 * @brief Release column data
 ***************************************************************************/
void GFitsTableCDoubleCol::release_data(void)
{
    // Free any existing memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free and reset loaded vector size
    m_data = NULL;
    m_size = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocates null value
 ***************************************************************************/
void GFitsTableCDoubleCol::alloc_nulval(const GFits::cdouble* value)
{
    // Free any existing memory
    if (m_nulval != NULL) delete m_nulval;

    // Mark pointer as free
    m_nulval = NULL;

    // If we have valid value, allocate and set nul value
    if (value != NULL) {
        m_nulval  = new GFits::cdouble;
        *m_nulval = *value;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise column data
 ***************************************************************************/
void GFitsTableCDoubleCol::init_data(void)
{
    // Initialise data if they exist
    if (m_data != NULL) {
        for (int i = 0; i < m_size; ++i) {
            m_data[i].re = 0.0;
            m_data[i].im = 0.0;
        }
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
