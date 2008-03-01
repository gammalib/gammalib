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

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableCol::GFitsTableCol()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***************************************************************************
 *                              Copy constructor                           *
 * ----------------------------------------------------------------------- *
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


/***************************************************************************
 *                               Destructor                                *
 * ----------------------------------------------------------------------- *
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

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
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


/*==========================================================================
 =                                                                         =
 =                        GFitsTableCol private methods                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableCol::init_members(void)
{
    // Initialise members
    m_name.clear();
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


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableCol::copy_members(const GFitsTableCol& column)
{
    // Copy attributes
    m_name     = column.m_name;
    m_colnum   = column.m_colnum;
    m_type     = column.m_type;
    m_repeat   = column.m_repeat;
    m_width    = column.m_width;
    m_length   = column.m_length;
    m_fitsfile = column.m_fitsfile;
    
    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
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
 * @param fptr FITS file pointer to which the table column should be connected
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

/***************************************************************************
 *                             Output operator                             *
 * ----------------------------------------------------------------------- *
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
