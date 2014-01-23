/***************************************************************************
 *                   GFitsHDU.cpp  - FITS HDU handling class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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
 * @file GFitsHDU.cpp
 * @brief GFitsHDU class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsCfitsio.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsHeaderCard.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPEN                                          "GFitsHDU::open(int)"
#define G_SAVE                                             "GFitsHDU::save()"
#define G_HEADER                                         "GFitsHDU::header()"
#define G_DATA                                             "GFitsHDU::data()"
#define G_COLUMN                             "GFitsHDU::column(std::string&)"
#define G_CONNECT                                  "GFitsHDU::connect(void*)"
#define G_MOVE_TO_HDU                               "GFitsHDU::move_to_hdu()"
#define G_GET_HDU_TYPE                             "GFitsHDU::get_hdu_type()"
#define G_NEW_IMAGE                                   "GFitsHDU::new_image()"

/* __ Definitions ________________________________________________________ */
#define DEBUG 0

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
GFitsHDU::GFitsHDU(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] hdu HDU from which the instance should be constructed.
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
GFitsHDU::~GFitsHDU(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] hdu HDU which should be assigned.
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set HDU extension name (EXTNAME keyword)
 *
 * @param[in] extname Name of HDU.
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
    //if (m_header == NULL) m_header = new GFitsHeader();

    // Update header
    m_header.append(GFitsHeaderCard("EXTNAME", extname,
                                    "name of this extension"));

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Connect HDU to FITS file
 *
 * @param[in] vptr FITS file pointer.
 *
 * Connects the HDU to the file specified by the FITS file pointer. Sets
 * also the HDU number (or extension number, starting from 0). This method
 * does nothing if the file pointer in not valid.
 ***************************************************************************/
void GFitsHDU::connect(void* vptr)
{
    // Continue only if file pointer is valid
    if (vptr != NULL) {

        // Connect HDU by copying the file pointer
        FPTR_COPY(m_fitsfile, vptr);

        // Extract HDU number from file pointer
        m_hdunum = FPTR(m_fitsfile)->HDUposition;

        // Connect data
        data_connect(vptr);

    } // endif: file pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Move FITS file pointer to HDU
 *
 * @exception GException::fits_file_not_open
 *            No FITS file has been opened.
 *
 * Moves to FITS file pointer to the actual HDU. This operation should
 * preceed any FITS file manipulation.
 ***************************************************************************/
void GFitsHDU::move_to_hdu(void)
{
    // Throw an exception if FITS file is not open
    if (FPTR(m_fitsfile)->Fptr == NULL) {
        throw GException::fits_file_not_open(G_MOVE_TO_HDU, 
              "Open file before moving file pointer to HDU.");
    }

    // Move to HDU
    gammalib::fits_move_to_hdu(G_MOVE_TO_HDU, m_fitsfile, m_hdunum+1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get HDU type from FITS file
 *
 * @exception GException::fits_file_not_open
 *            No FITS file has been opened.
 * @exception GException::fits_hdu_not_found
 *            Requested HDU not found.
 ***************************************************************************/
GFitsHDU::HDUType GFitsHDU::get_hdu_type(void) const
{
    // Initialise HDU type
    int type = 0;

    // Throw an exception if FITS file is not open
    if (FPTR(m_fitsfile)->Fptr == NULL) {
        throw GException::fits_file_not_open(G_GET_HDU_TYPE, 
              "Open file before requesting HDU type.");
    }

    // Get HDU type
    int status = 0;
    status     = __ffghdt(FPTR(m_fitsfile), &type, &status);
    if (status != 0) {
        throw GException::fits_error(G_GET_HDU_TYPE, status);
    }

    // Return HDU type
    return (HDUType)type;
}


/***********************************************************************//**
 * @brief Open HDU
 *
 * @param[in] vptr FITS file pointer.
 * @param[in] hdunum Number of HDU (starting from 0).
 *
 * @exception GException::fits_file_not_open
 *            FITS file pointer does not point to an open FITS file
 * @exception GException::fits_error
 *            Unable to open FITS HDU.
 *
 * Opens an (existing) HDU in the FITS file. This method does NOT create any
 * HDU if it does not exist. Opening consists of fetching all header cards
 * (by opening an associated GFitsHeader instance) and of opening the data
 * area (which can be of type Image or Table)
 ***************************************************************************/
void GFitsHDU::open(void* vptr, int hdunum)
{
    // Verify that FITS file pointer is valid
    if (vptr == NULL) {
        throw GException::fits_file_not_open(G_OPEN,
              "FITS file pointer does not point to an open FITS file.");
    }

    // Move to HDU
    gammalib::fits_move_to_hdu(G_OPEN, vptr, hdunum+1);

    // Save the FITS file pointer and the HDU number
    FPTR_COPY(m_fitsfile, vptr);
    m_hdunum = hdunum;

    // Load HDU header
    m_header.load(FPTR(m_fitsfile));

    // Open HDU data area
    data_open(FPTR(m_fitsfile));

    // Get HDU name from header. If no name was found and this is the primary
    // HDU then set the name to "Primary", otherwise to "NoName".
    if (m_header.contains("EXTNAME")) {
        m_name = gammalib::strip_whitespace(m_header.string("EXTNAME"));
    }
    else {
        m_name.clear();
    }
    if (m_name.length() == 0) {
        if (hdunum == 0) {
            m_name = "Primary";
        }
        else {
            m_name = "NoName";
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Saves HDU
 *
 * Save the HDU to the FITS file. In case that the HDU does not exist in the
 * file it will be created (by the save_data() method). In the special case
 * that no first HDU exists, an empty primary image is created.
 *
 * @todo Put the m_header.save(FPTR(m_fitsfile)) call in the data_save()
 * methods
 ***************************************************************************/
void GFitsHDU::save(void)
{
    // Debug header
    #if DEBUG
    std::cout << "GFitsHDU::save() -->" << std::endl;
    #endif

    // Save data. This method allows saving past the file and will
    // physically append a new HDU to a FITS file. It has to be called
    // before the saving of the header cards, since header card
    // saving requires the HDU to exist in the FITS file.
    data_save();

    // Save header cards (assumes implicitely that the HDU exists in
    // FITS file)
    m_header.save(FPTR(m_fitsfile));

    // Debug trailer
    #if DEBUG
    std::cout << "<-- GFitsHDU::save" << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print basic HDU information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing basic HDU information.
 ***************************************************************************/
std::string GFitsHDU::print_hdu(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append HDU information
        result.append(gammalib::parformat("HDU number"));
        result.append(gammalib::str(m_hdunum)+"\n");
        result.append(gammalib::parformat("HDU name")+m_name);

    } // endif: chatter was not silent

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return typecode as string
 *
 * @param[in] type Type code.
 ***************************************************************************/
std::string GFitsHDU::typecode(int type) const
{
    // Allocate string
    std::string result;

    // Set typecode
    switch (type) {
    case __TBIT:
        result = "bit";
        break;
    case __TBYTE:
        result = "unsigned byte";
        break;
    case __TSBYTE:
        result = "signed byte";
        break;
    case __TLOGICAL:
        result = "boolean";
        break;
    case __TSTRING:
        result = "string";
        break;
    case __TUSHORT:
        result = "unsigned short integer";
        break;
    case __TSHORT:
        result = "short integer";
        break;
    case __TULONG:
        result = "unsigned long integer";
        break;
    case __TLONG:
        result = "long integer";
        break;
    case __TFLOAT:
        result = "single precision floating point";
        break;
    case __TLONGLONG:
        result = "long long integer";
        break;
    case __TDOUBLE:
        result = "double precision floating point";
        break;
    default:
        result = "unsupported format";
        break;
    }

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
 *
 * Sets all class members to well defined values.
 ***************************************************************************/
void GFitsHDU::init_members(void)
{
    // Allocate FITS file pointer. As the pointer is of type void the C-style
    // memory allocation function malloc is used as this function does not
    // deal with types but always returns a void pointer
    m_fitsfile = malloc(sizeof(__fitsfile));

    // Initialise FITS file pointer
    FPTR(m_fitsfile)->HDUposition = 0;
    FPTR(m_fitsfile)->Fptr        = NULL;

    // Initialise members
    m_hdunum = 0;
    m_name.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] hdu HDU to be copied.
 *
 * Assumes that all memory has been freed correctly before calling.
 ***************************************************************************/
void GFitsHDU::copy_members(const GFitsHDU& hdu)
{
    // Copy members
    FPTR_COPY(m_fitsfile, hdu.m_fitsfile);
    m_hdunum = hdu.m_hdunum;
    m_name   = hdu.m_name;
    m_header = hdu.m_header;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsHDU::free_members(void)
{
    // Free memory. Note that the malloc function was used for allocation,
    // hence the free function needs to be called for freeing.
    if (m_fitsfile != NULL) free(m_fitsfile);

    // Signal free pointers
    m_fitsfile = NULL;

    // Return
    return;
}
