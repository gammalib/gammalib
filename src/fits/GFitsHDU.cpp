/***************************************************************************
 *                   GFitsHDU.cpp  - FITS HDU handling class               *
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
#include "GException.hpp"
#include "GFitsHDU.hpp"
#include <iostream>                           // cout, cerr

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                     GFitsHDU constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHDU::GFitsHDU()
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
GFitsHDU::GFitsHDU(const GFitsHDU& hdu)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(hdu);

    // Return
    return;
}


/***************************************************************************
 *                               Destructor                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHDU::~GFitsHDU()
{
    // Free members
    free_members();

    // Return
    return;
}

/*==========================================================================
 =                                                                         =
 =                           GFitsHDU operators                            =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHDU& GFitsHDU::operator= (const GFitsHDU& hdu)
{
    // Execute only if object is not identical
    if (this != &hdu) {
  
        // Free members
        free_members();
  
        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(hdu);
	
    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsHDU public methods                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                 Open HDU                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::open(__fitsfile* fptr, int hdunum)
{
    // Store information
    m_fitsfile = fptr;
    m_num      = hdunum;
    
    // Move to HDU
    move2hdu();
    
    // Get HDU type
    int status = 0;
    status     = __ffghdt(m_fitsfile, &m_type, &status);
    if (status != 0) {
        throw GException::fits_error("GFitsHDU::open(int)", status);
    }
    
    // Open HDU header
    m_header.open(m_fitsfile);
    
    // Open HDU data area
    //
    
    // Get HDU name from header
    GFitsHeaderCard* extname = m_header.card("EXTNAME");
    if (extname == NULL) {
        if (hdunum == 1)
            m_name = "Primary";
        else
            m_name = "NoName";
    }
    else
        m_name = extname->value();

    cout << "HDU #" << hdunum << ": " << m_type << " : " << m_name << endl;

    // Return
    return;
}


/***************************************************************************
 *                                Close HDU                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::close(void)
{
    // Return
    return;
}


/***************************************************************************
 *                       Return pointer to HDU header                      *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHeader* GFitsHDU::header(void)
{
    // Return header pointer
    return &m_header;
}


/***************************************************************************
 *                        Return pointer to HDU data                       *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsData* GFitsHDU::data(void)
{
    // Return data pointer
    return &m_data;
}


/*==========================================================================
 =                                                                         =
 =                         GFitsHDU private methods                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::init_members(void)
{
    // Initialise members
    m_fitsfile = NULL;
    m_type     = 0;
    //m_header   = NULL;
    //m_data     = NULL;
  
    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::copy_members(const GFitsHDU& hdu)
{
    // Copy membres
    m_fitsfile = hdu.m_fitsfile;
    m_type     = hdu.m_type;
    m_header   = hdu.m_header;
    m_data     = hdu.m_data;
    
    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::free_members(void)
{
    // Return
    return;
}


/***************************************************************************
 *                      Move FITS file pointer to HDU                      *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::move2hdu(void)
{
    // Move FITS file pointer to HDU
    int status = 0;
    status     = __ffmahd(m_fitsfile, m_num, NULL, &status);
    if (status != 0) {
        throw GException::fits_error("GFitsHDU::move2hdu(void)", status);
    }
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GFitsHDU friends                            =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GFitsHDU                    =
 =                                                                         =
 ==========================================================================*/
