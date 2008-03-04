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
#include <iostream>
#include "GException.hpp"
#include "GFitsTableCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

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
 * @param[in] length Length of column.
 * @param[in] size Vector size of column.
 ***************************************************************************/
GFitsTableCol::GFitsTableCol(const int& length, const int& size)
{
    // Initialise class members for clean destruction
    init_members();

    // Set column length and vector size
    m_length = length;
    m_repeat = size;

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
 * @brief Set column number
 *
 * @param[in] colnum Number of column.
 ***************************************************************************/
/*
void GFitsTableCol::colnum(const int& colnum)
{
    // Set column number
    m_colnum = colnum;

    // Return
    return;
}
*/

/***********************************************************************//**
 * @brief Set column type
 *
 * @param[in] type Type of column.
 ***************************************************************************/
/*
void GFitsTableCol::type(const int& type)
{
    // Set column type
    m_type = type;

    // Return
    return;
}
*/

/***********************************************************************//**
 * @brief Set column length
 *
 * @param[in] length Length of column.
 ***************************************************************************/
/*
void GFitsTableCol::length(const int& length)
{
    // Set column length
    m_length = length;

    // Return
    return;
}
*/

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
 * @brief Get column repeat value
 ***************************************************************************/
int GFitsTableCol::repeat(void)
{
    // Return column repeat value
    return m_repeat;
}


/***********************************************************************//**
 * @brief Get column width
 ***************************************************************************/
int GFitsTableCol::width(void)
{
    // Return column width
    return m_width;
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
    m_format.clear();
    m_unit.clear();
    m_colnum               = 0;
    m_type                 = 0;
    m_repeat               = 0;
    m_width                = 0;
    m_length               = 0;
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
    m_format   = column.m_format;
    m_unit     = column.m_unit;
    m_colnum   = column.m_colnum;
    m_type     = column.m_type;
    m_repeat   = column.m_repeat;
    m_width    = column.m_width;
    m_length   = column.m_length;
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


/*==========================================================================
 =                                                                         =
 =                           GFitsTableCol friends                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] column Column to put in output stream
 ***************************************************************************/
ostream& operator<< (ostream& os, const GFitsTableCol& column)
{
    // Put header in stream
    os << "'" << column.m_name << "'";
    os << " [" << column.m_colnum << "] ";
    os << column.m_type << " repeat=";
    os << column.m_repeat << " width=";
    os << column.m_width << " length=";
    os << column.m_length << endl;

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                    Other functions used by GFitsTableCol                =
 =                                                                         =
 ==========================================================================*/
