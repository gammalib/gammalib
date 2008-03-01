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
#include <iostream>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsHDU.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_OPEN     "GFitsHDU::open(int)"
#define G_SAVE     "GFitsHDU::save()"
#define G_COLUMN   "GFitsHDU::column(const std::string&)"

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
    // Move to HDU
    int status = 0;
    status     = __ffmahd(fptr, hdunum, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Store information. 
    // Note that the HDU number - 1 is stored in m_fitsfile->HDUposition !!!
    m_fitsfile = *fptr;

    // Get HDU type
    status = __ffghdt(&m_fitsfile, &m_type, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Allocate and open HDU header
    if (m_header != NULL) delete m_header;
    m_header = new GFitsHeader();
    m_header->open(&m_fitsfile);

    // Open HDU data area
    if (m_data != NULL) delete m_data;
    switch (m_type) {
    case 0:        // Image HDU
        m_data = new GFitsImage();
        break;
    case 1:        // ASCII Table HDU
        m_data = new GFitsAsciiTable();
        break;
    case 2:        // Binary Table HDU
        m_data = new GFitsBinTable();
        break;
    default:
        throw GException::fits_unknown_HDU_type(G_OPEN, m_type);
        break;
    }
    m_data->open(&m_fitsfile);

    // Get HDU name from header
    m_name = strip_whitespace(m_header->string("EXTNAME"));
    if (m_name.length() == 0) {
        if (hdunum == 1)
            m_name = "Primary";
        else
            m_name = "NoName";
    }

    // Return
    return;
}


/***************************************************************************
 *                           Save HDU to FITS file                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/// NOT YET IMPLEMENTED !!!
void GFitsHDU::save(void)
{
cout << "entry" << endl;
    // Move to HDU
    int status = 0;
    status     = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL, &status);
    cout << status << endl;
    
    // If HDU does not exist then create it now
    
    
    //if (status != 0)
    //    throw GException::fits_error(G_SAVE, status);

    // Save header
    //...

    // Save data
    //...

cout << "exit" << endl;
    // Return
    return;
}


/***************************************************************************
 *                     Return pointer to column of table                   *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableCol* GFitsHDU::column(const std::string& colname) const
{
    // Initialise pointer
    GFitsTableCol* ptr = NULL;

    // Get pointer to table column
    switch (m_type) {
    case 0:        // Image HDU
        throw GException::fits_HDU_not_a_table(G_COLUMN, m_type);
        break;
    case 1:        // ASCII Table HDU
        ptr = ((GFitsAsciiTable*)this->data())->column(colname);
        break;
    case 2:        // Binary Table HDU
        ptr = ((GFitsBinTable*)this->data())->column(colname);
        break;
    default:
        throw GException::fits_unknown_HDU_type(G_COLUMN, m_type);
        break;
    }

    // Return pointer
    return ptr;
}


/***************************************************************************
 *                         Setup minimal primary HDU                       *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::primary(void)
{
    // Free any allocated header and data
    free_members();
    
    // Allocate header
    m_header = new GFitsHeader();
    
    // Add cards
    m_header->update(GFitsHeaderCard("SIMPLE", "T", "/ file does conform to FITS standard"));
    m_header->update(GFitsHeaderCard("BITPIX", "8", "/ number of bits per data pixel"));
    m_header->update(GFitsHeaderCard("NAXIS", "0", "/ number of data axes"));
    m_header->update(GFitsHeaderCard("EXTEND", "T", "/ FITS dataset may contain extensions"));
    m_header->update(GFitsHeaderCard("BZERO", "-128", "/ Make values Signed"));
    m_header->update(GFitsHeaderCard("BSCALE", "1", "/ Make values Signed"));
    
    // Return
    return;
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
    m_fitsfile.HDUposition = 0;
    m_fitsfile.Fptr        = NULL;
    m_name.clear();
    m_type                 = 0;
    m_header               = NULL;
    m_data                 = NULL;

    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::copy_members(const GFitsHDU& hdu)
{
    // Copy members
    m_fitsfile = hdu.m_fitsfile;
    m_name     = hdu.m_name;
    m_type     = hdu.m_type;
    if (hdu.m_header != NULL) m_header = hdu.m_header->clone();
    if (hdu.m_data   != NULL) m_data   = hdu.m_data->clone();

    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsHDU::free_members(void)
{
    // Free memory
    if (m_header != NULL) delete m_header;
    if (m_data   != NULL) delete m_data;

    // Signal free pointers
    m_header = NULL;
    m_data   = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GFitsHDU friends                            =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
ostream& operator<< (ostream& os, const GFitsHDU& hdu)
{
    // Put header in stream
    os << "=== GFitsHDU ===" << endl;
    os << " HDU number ................: " << hdu.m_fitsfile.HDUposition << endl;
    os << " HDU name ..................: " << hdu.m_name << endl;
    os << " HDU type ..................: ";
    switch (hdu.m_type) {
    case 0:        // Image HDU
        os << "Image" << endl;
        break;
    case 1:        // ASCII Table HDU
        os << "ASCII table" << endl;
        break;
    case 2:        // Binary Table HDU
        os << "Binary table" << endl;
        break;
    default:
        os << "Unknown" << endl;
        break;
    }
    
    // Put FITS header in stream
    if (hdu.m_header != NULL)
        os << *hdu.m_header;
        
    // Put FITS data in stream
    if (hdu.m_data != NULL) {
        switch (hdu.m_type) {
        case 0:        // Image HDU
            os << *(GFitsImage*)hdu.m_data;
            break;
        case 1:        // ASCII Table HDU
            os << *(GFitsAsciiTable*)hdu.m_data;
            break;
        case 2:        // Binary Table HDU
            os << *(GFitsBinTable*)hdu.m_data;
            break;
        default:
            break;
        }
    }

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GFitsHDU                    =
 =                                                                         =
 ==========================================================================*/
