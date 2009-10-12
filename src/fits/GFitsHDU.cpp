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
/**
 * @file GFitsHDU.cpp
 * @brief GFitsHDU class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsHDU.hpp"
#include "GFitsImageFlt.hpp"
#include "GFitsImageDbl.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPEN      "GFitsHDU::open(int)"
#define G_SAVE      "GFitsHDU::save()"
#define G_HEADER    "GFitsHDU::header()"
#define G_DATA      "GFitsHDU::data()"
#define G_COLUMN    "GFitsHDU::column(const std::string&)"
#define G_NEW_IMAGE "GFitsHDU::new_image()"

/* __ Definitions ________________________________________________________ */
#define HT_IMAGE       0
#define HT_ASCII_TABLE 1
#define HT_BIN_TABLE   2

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                     GFitsHDU constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsHDU::GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Image HDU Constructor
 *
 * @param[in] image Image from which the HDU should be constructed
 *
 * Note that the HDU constructor copies the image, hence any change in the
 * original object after copying will not be reflected in the copied object.
 ***************************************************************************/
GFitsHDU::GFitsHDU(const GFitsImage& image)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy image into HDU data and set data type
    m_data = image.clone();
    m_type = HT_IMAGE;

    // Allocate new header
    m_header = new GFitsHeader();

    // Set image header keywords
    m_header->update(GFitsHeaderCard("XTENSION", "'IMAGE   '",
                                     "IMAGE extension"));
    m_header->update(GFitsHeaderCard("BITPIX", image.bitpix(),
                                     "number of bits per data pixel"));
    m_header->update(GFitsHeaderCard("NAXIS", image.naxis(),
                                     "number of data axes"));
    for (int i = 0; i < image.naxis(); ++i) {
        std::ostringstream s_key;
        std::ostringstream s_comment;
        s_key     << "NAXIS" << (i+1);
        s_comment << "length of data axis " << (i+1);
        m_header->update(GFitsHeaderCard(s_key.str(), image.naxes(i),
                                         s_comment.str()));
    }
    m_header->update(GFitsHeaderCard("PCOUNT", 0,
                                     "required keyword; must = 0"));
    m_header->update(GFitsHeaderCard("GCOUNT", 1,
                                     "required keyword; must = 1"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief ASCII Table HDU Constructor
 *
 * @param[in] table Table from which the HDU should be constructed
 ***************************************************************************/
GFitsHDU::GFitsHDU(const GFitsAsciiTable& table)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy table into HDU data and set data type
    m_data = table.clone();
    m_type = HT_ASCII_TABLE;

    // Allocate new header
    m_header = new GFitsHeader();

    // Set image header keywords
    /*
    m_header->update(GFitsHeaderCard("XTENSION", "'BINTABLE'",
                                     "binary table extension"));
    m_header->update(GFitsHeaderCard("BITPIX", 8,
                                     "8-bit bytes"));
    m_header->update(GFitsHeaderCard("NAXIS", 2,
                                     "2-dimensional binary table"));
    m_header->update(GFitsHeaderCard("NAXIS1", 0,
                                     "width of table in bytes"));
    m_header->update(GFitsHeaderCard("NAXIS2", table.rows(),
                                     "number of rows in table"));
    m_header->update(GFitsHeaderCard("PCOUNT", 0,
                                     "size of special data area"));
    m_header->update(GFitsHeaderCard("GCOUNT", 1,
                                     "one data group (required keyword)"));
    m_header->update(GFitsHeaderCard("TFIELDS", 1,
                                     "number of fields in each row"));
*/
    // Return
    return;
}


/***********************************************************************//**
 * @brief Binary Table HDU Constructor
 *
 * @param[in] table Table from which the HDU should be constructed
 ***************************************************************************/
GFitsHDU::GFitsHDU(const GFitsBinTable& table)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy table into HDU data and set data type
    m_data = table.clone();
    m_type = HT_BIN_TABLE;

    // Allocate new header
    m_header = new GFitsHeader();

    // Set image header keywords
    /*
    m_header->update(GFitsHeaderCard("XTENSION", "'BINTABLE'",
                                     "binary table extension"));
    m_header->update(GFitsHeaderCard("BITPIX", 8,
                                     "8-bit bytes"));
    m_header->update(GFitsHeaderCard("NAXIS", 2,
                                     "2-dimensional binary table"));
    m_header->update(GFitsHeaderCard("NAXIS1", 0,
                                     "width of table in bytes"));
    m_header->update(GFitsHeaderCard("NAXIS2", table.rows(),
                                     "number of rows in table"));
    m_header->update(GFitsHeaderCard("PCOUNT", 0,
                                     "size of special data area"));
    m_header->update(GFitsHeaderCard("GCOUNT", 1,
                                     "one data group (required keyword)"));
    m_header->update(GFitsHeaderCard("TFIELDS", 1,
                                     "number of fields in each row"));
*/
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] hdu HDU from which the instance should be constructed
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


/***********************************************************************//**
 * @brief Destructor
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

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] hdu HDU which should be assigned
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

/***********************************************************************//**
 * @brief Get HDU extension name (EXTNAME keyword)
 ***************************************************************************/
std::string GFitsHDU::extname(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Set HDU extension name (EXTNAME keyword)
 *
 * @param[in] extname Name of HDU
 *
 * This method sets the extension name of the HDU. The extension name will
 * be saved in the 'EXTNAME' header keyword. The header attached to the
 * HDU will be automatically updated by this method.
 ***************************************************************************/
void GFitsHDU::extname(const std::string& extname)
{
    // Set extension name
    m_name = extname;

    // If no header exists then create one now
    if (m_header == NULL) m_header = new GFitsHeader();

    // Update header
    m_header->update(GFitsHeaderCard("EXTNAME", extname,
                                     "name of this extension"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get HDU extension number
***************************************************************************/
int GFitsHDU::extno(void) const
{ 
    return m_hdunum; 
}


/***********************************************************************//**
 * @brief Get HDU extension type
***************************************************************************/
int GFitsHDU::exttype(void) const
{
    return m_type;
}


/***********************************************************************//**
 * @brief Get pointer to header
 *
 * @exception GException::fits_no_header
 *            No header was found in HDU
 ***************************************************************************/
GFitsHeader* GFitsHDU::header(void) const
{
    // If no header pointer is available then throw exception
    if (m_header == NULL)
        throw GException::fits_no_header(G_HEADER, "No header found");

    // Return header pointer
    return m_header;
}


/***********************************************************************//**
 * @brief Get pointer to data
 *
 * @exception GException::fits_no_data
 *            No data was found in HDU
 ***************************************************************************/
GFitsData* GFitsHDU::data(void) const
{
    // If no header pointer is available then throw exception
    if (m_data == NULL)
        throw GException::fits_no_data(G_DATA, "No data found");

    // Return data pointer
    return m_data;
}


/***********************************************************************//**
 * @brief Get pointer to header card
 *
 * @param[in] keyname Name of header card
 ***************************************************************************/
GFitsHeaderCard* GFitsHDU::card(const std::string& keyname) const
{
    return m_header->card(keyname);
}


/***********************************************************************//**
 * @brief Get pointer to header card
 *
 * @param[in] cardno Number of card in header
 ***************************************************************************/
GFitsHeaderCard* GFitsHDU::card(const int& cardno) const
{
    return m_header->card(cardno);
}


/***********************************************************************//**
 * @brief Return pointer to column of table
 *
 * @param[in] colname Name of FITS table column to be returned
 *
 * @exception GException::fits_HDU_not_a_table
 *            HDU is not a table
 * @exception GException::fits_unknown_HDU_type
 *            HDU type is not known
 *
 * If this method is called for an image a 'fits_HDU_not_a_table' exception
 * will be thrown.
 ***************************************************************************/
GFitsTableCol* GFitsHDU::column(const std::string& colname) const
{
    // Initialise pointer
    GFitsTableCol* ptr = NULL;

    // Get pointer to table column
    switch (m_type) {
    case HT_IMAGE:              // Image HDU
        throw GException::fits_HDU_not_a_table(G_COLUMN, m_type);
        break;
    case HT_ASCII_TABLE:        // ASCII Table HDU
        ptr = ((GFitsAsciiTable*)this->data())->column(colname);
        break;
    case HT_BIN_TABLE:          // Binary Table HDU
        ptr = ((GFitsBinTable*)this->data())->column(colname);
        break;
    default:
        throw GException::fits_unknown_HDU_type(G_COLUMN, m_type);
        break;
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Return pointer to column of table
 *
 * @param[in] colnum Number FITS table column to be returned (starting from 0)
 *
 * If this method is called for an image a 'fits_HDU_not_a_table' exception
 * will be thrown.
 ***************************************************************************/
GFitsTableCol* GFitsHDU::column(const int& colnum) const
{
    // Initialise pointer
    GFitsTableCol* ptr = NULL;

    // Get pointer to table column
    switch (m_type) {
    case HT_IMAGE:              // Image HDU
        throw GException::fits_HDU_not_a_table(G_COLUMN, m_type);
        break;
    case HT_ASCII_TABLE:        // ASCII Table HDU
        ptr = ((GFitsAsciiTable*)this->data())->column(colnum);
        break;
    case HT_BIN_TABLE:          // Binary Table HDU
        ptr = ((GFitsBinTable*)this->data())->column(colnum);
        break;
    default:
        throw GException::fits_unknown_HDU_type(G_COLUMN, m_type);
        break;
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Setup minimal primary HDU
 ***************************************************************************/
void GFitsHDU::primary(void)
{
    // Free any allocated header and data
    free_members();

    // Create primary image in memory
    int status = 0;
    __fitsfile* fptr;
    status = __ffinit(&fptr, "mem://", &status);
    status = __ffcrim(fptr, 8, 0, NULL, &status);

    // Open HDU
    this->open(fptr,1);

    // Detach FITS file pointer
    this->m_fitsfile.HDUposition = 0;
    this->m_fitsfile.Fptr        = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GFitsHDU private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 *
 * Sets all class members to well defined values.
 ***************************************************************************/
void GFitsHDU::init_members(void)
{
    // Initialise members
    m_fitsfile.HDUposition = 0;
    m_fitsfile.Fptr        = NULL;
    m_hdunum               = 0;
    m_name.clear();
    m_type                 = 0;
    m_header               = NULL;
    m_data                 = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] hdu HDU to be copied
 *
 * Assumes that all memory has been freed correctly before calling.
 ***************************************************************************/
void GFitsHDU::copy_members(const GFitsHDU& hdu)
{
    // Copy members
    m_fitsfile = hdu.m_fitsfile;
    m_hdunum   = hdu.m_hdunum;
    m_name     = hdu.m_name;
    m_type     = hdu.m_type;
    if (hdu.m_header != NULL) m_header = hdu.m_header->clone();
    if (hdu.m_data   != NULL) m_data   = hdu.m_data->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
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


/***********************************************************************//**
 * @brief Connect HDU to FITS file
 *
 * @param[in] fptr FITS file pointer to which the HDU should be connected
 *
 * Connects the HDU and the associated data area to the specified FITS file.
 * Note that the header area does not require a FITS pointer, hence no
 * connection is required.
 ***************************************************************************/
void GFitsHDU::connect(__fitsfile* fptr)
{
    // First connect HDU
    m_fitsfile = *fptr;
    m_hdunum   = m_fitsfile.HDUposition;

    // Then connect data
    if (m_data != NULL) m_data->connect(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate new FITS image and return memory pointer
 *
 * Depending on the number of bits per pixel, a FITS image is allocated
 * and the pointer is returned. The following FITS image classes are
 * handled:
 * GFitsBytImage (bitpix=8)
 * GFitsShtImage (bitpix=16)
 * GFitsLngImage (bitpix=32)
 * GFitsLlgImage (bitpix=64)
 * GFitsFltImage (bitpix=-32)
 * GFitsDblImage (bitpix=-64)
 * The information about the number of bits per pixels is extracted from
 * the actual HDU.
 ***************************************************************************/
GFitsImage* GFitsHDU::new_image(void)
{
    // Initialise return value
    GFitsImage* image = NULL;

    // Get number of bits per pixel
    int status =   0;
    int bitpix = -64;
    status     = __ffgipr(&m_fitsfile, 0, &bitpix, NULL, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_NEW_IMAGE, status);

    // Allocate bitpix dependent image
    switch (bitpix) {
    case 8:
        image = new GFitsImageDbl();  // TO BE REPLACED BY CORRECT CLASS
        break;
    case 16:
        image = new GFitsImageDbl();  // TO BE REPLACED BY CORRECT CLASS
        break;
    case 32:
        image = new GFitsImageDbl();  // TO BE REPLACED BY CORRECT CLASS
        break;
    case 64:
        image = new GFitsImageDbl();  // TO BE REPLACED BY CORRECT CLASS
        break;
    case -32:
        image = new GFitsImageFlt();
        break;
    case -64:
        image = new GFitsImageDbl();
        break;
    default:
        throw GException::fits_bad_bitpix(G_OPEN, m_type);
        break;
    }

    // Return image pointer
    return image;
}


/***********************************************************************//**
 * @brief Open HDU
 *
 * @param[in] fptr FITS file pointer
 * @param[in] hdunum Number of HDU (starting from 0)
 *
 * Opens an (existing) HDU in the FITS file. This method does NOT create any
 * HDU if it does not exist. Opening consists of fetching all header cards
 * (by opening an associated GFitsHeader instance) and of opening the data
 * area (which can be of type Image or Table)
 ***************************************************************************/
void GFitsHDU::open(__fitsfile* fptr, int hdunum)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(fptr, hdunum, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Store information. 
    // Note that the HDU number is stored in m_fitsfile->HDUposition !!!
    m_fitsfile = *fptr;
    m_hdunum   = hdunum;

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
    case HT_IMAGE:              // Image HDU
        m_data = new_image();
        break;
    case HT_ASCII_TABLE:        // ASCII Table HDU
        m_data = new GFitsAsciiTable();
        break;
    case HT_BIN_TABLE:          // Binary Table HDU
        m_data = new GFitsBinTable();
        break;
    default:
        throw GException::fits_unknown_HDU_type(G_OPEN, m_type);
        break;
    }
    m_data->open(&m_fitsfile);

    // Get HDU name from header
    try {
        m_name = strip_whitespace(m_header->string("EXTNAME"));
    }
    catch (GException::fits_key_not_found) {
        m_name.clear();
    }
    if (m_name.length() == 0) {
        if (hdunum == 1)
            m_name = "Primary";
        else
            m_name = "NoName";
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Saves HDU
 *
 * Save the HDU to the FITS file. In case that the HDU does not exist in the
 * file it will be created (by the corresponding opening routines of the
 * GFitsData instance).
 * In the special case that no first HDU exists an empty primary image is
 * created.
 ***************************************************************************/
void GFitsHDU::save(void)
{
//cout << "GFitsHDU::save entry" << endl;

    // Save data and header if they exist
    if (m_data != NULL) { 
        m_data->connect(&m_fitsfile);
        m_data->save();
        if (m_header != NULL)
            m_header->save(&m_fitsfile);
    }

    // ... otherwise make sure that the FITS file has at least a primary
    // image.
    else {

        // Special case: The first HDU does not exist. A FITS file needs at
        // least one primary HDU. Thus an empty image is created.
        if (m_fitsfile.HDUposition == 0) {
            int status = 0;
            status     = __ffcrim(&m_fitsfile, 8, 0, NULL, &status);
            if (status != 0)
                throw GException::fits_error(G_SAVE, status);
            if (m_header != NULL)
                m_header->save(&m_fitsfile);
        }

    } // endelse: no data were available

//cout << "GFitsHDU::save exit" << endl;
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GFitsHDU friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream into which the HDU will be dumped
 * @param[in] hdu HDU to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GFitsHDU& hdu)
{
    // Put header in stream
    os << "=== GFitsHDU ===" << std::endl;
    os << " HDU number ................: " << hdu.m_hdunum << std::endl;
    os << " HDU name ..................: " << hdu.m_name << std::endl;
    os << " HDU type ..................: ";
    switch (hdu.m_type) {
    case 0:        // Image HDU
        os << "Image" << std::endl;
        break;
    case 1:        // ASCII Table HDU
        os << "ASCII table" << std::endl;
        break;
    case 2:        // Binary Table HDU
        os << "Binary table" << std::endl;
        break;
    default:
        os << "Unknown" << std::endl;
        break;
    }

    // Put FITS header in stream
    if (hdu.m_header != NULL)
        os << *hdu.m_header;

    // Put FITS data in stream
    if (hdu.m_data != NULL) {
        switch (hdu.m_type) {
        case HT_IMAGE:              // Image HDU
            os << *(GFitsImage*)hdu.m_data;
            break;
        case HT_ASCII_TABLE:        // ASCII Table HDU
            os << *(GFitsAsciiTable*)hdu.m_data;
            break;
        case HT_BIN_TABLE:          // Binary Table HDU
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
