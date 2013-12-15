/***************************************************************************
 *                    GFits.cpp - FITS file access class                   *
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
 * @file GFits.cpp
 * @brief FITS file access class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>
#include "GException.hpp"
#include "GFitsCfitsio.hpp"
#include "GFits.hpp"
#include "GFitsImageByte.hpp"
#include "GFitsImageSByte.hpp"
#include "GFitsImageUShort.hpp"
#include "GFitsImageShort.hpp"
#include "GFitsImageULong.hpp"
#include "GFitsImageLong.hpp"
#include "GFitsImageLongLong.hpp"
#include "GFitsImageFloat.hpp"
#include "GFitsImageDouble.hpp"
#include "GFitsAsciiTable.hpp"
#include "GFitsBinTable.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT1                                               "GFits::at(int&)"
#define G_AT2                                       "GFits::at(std::string&)"
#define G_IMAGE1                                         "GFits::image(int&)"
#define G_IMAGE2                                 "GFits::image(std::string&)"
#define G_TABLE1                                         "GFits::table(int&)"
#define G_TABLE2                                 "GFits::table(std::string&)"
#define G_SET1                                  "GFits::set(int&, GFitsHDU&)"
#define G_SET2                          "GFits::set(std::string&, GFitsHDU&)"
#define G_INSERT1            "GFits::insert(const int& extno, GFitsHDU& hdu)"
#define G_INSERT2                "GFits::insert(std::string&, GFitsHDU& hdu)"
#define G_REMOVE1                                       "GFits::remove(int&)"
#define G_REMOVE2                               "GFits::remove(std::string&)"
#define G_OPEN                                    "GFits::open(std::string&)"
#define G_SAVE                                           "GFits::save(bool&)"
#define G_SAVETO                         "GFits::saveto(std::string&, bool&)"
#define G_FREE_MEM                                    "GFits::free_members()"
#define G_NEW_IMAGE                                      "GFits::new_image()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_DELETE_EMPTY_FITS_FILES          //!< Do not write empty FITS files

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GFits::GFits(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS file constructor
 *
 * @param[in] filename FITS file name.
 * @param[in] create Create file if it does not exist (default=false)?
 *
 * Construct an object by opening a FITS file. If the file does not exist it
 * will be created if @p create is set to true.
 ***************************************************************************/
GFits::GFits(const std::string& filename, const bool& create)
{
    // Initialise class members
    init_members();

    // Open specified FITS file
    open(filename, create);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
GFits::GFits(const GFits& fits)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFits::~GFits(void)
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
 * @brief  Assignment operator
 *
 * @param[in] fits FITS file.
 * @return FITS file.
 ***************************************************************************/
GFits& GFits::operator=(const GFits& fits)
{
    // Execute only if object is not identical
    if (this != &fits) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(fits);

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
 * @brief Clear FITS file
 ***************************************************************************/
void GFits::clear(void)
{
    // Close file and free all members 
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone FITS file
 *
 * @return Pointer to deep copy of FITS file.
 ***************************************************************************/
GFits* GFits::clone(void) const 
{
    return new GFits(*this);
}


/***********************************************************************//**
 * @brief Get pointer to HDU
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return Pointer to HDU.
 *
 * @exception GException::out_of_range
 *            Extension number is out of range.
 * @exception GException::fits_hdu_not_found
 *            No HDU found for specified extension number.
 *
 * Returns a pointer to the HDU with the specified extension number @p extno.
 * An exception is thrown if the HDU is not valid (i.e. NULL).
 ***************************************************************************/
GFitsHDU* GFits::at(const int& extno)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (extno < 0 || extno >= size()) {
        throw GException::out_of_range(G_AT1, "Extension number", extno, size());
    }
    #endif

    // Get HDU pointer
    GFitsHDU* ptr = m_hdu[extno];

    // Throw an error if HDU has not been found
    if (ptr == NULL) {
        std::string extname = "extno="+gammalib::str(extno);
        throw GException::fits_hdu_not_found(G_AT1, extname);
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to HDU (const version)
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return Pointer to HDU.
 *
 * @exception GException::out_of_range
 *            Extension number is out of range.
 * @exception GException::fits_hdu_not_found
 *            No HDU found for specified extension number.
 *
 * Returns a pointer to the HDU with the specified extension number @p extno.
 * An exception is thrown if the HDU is not valid (i.e. NULL).
 ***************************************************************************/
const GFitsHDU* GFits::at(const int& extno) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (extno < 0 || extno >= size()) {
        throw GException::out_of_range(G_AT1, "Extension number", extno, size());
    }
    #endif

    // Get HDU pointer
    GFitsHDU* ptr = m_hdu[extno];

    // Throw an error if HDU has not been found
    if (ptr == NULL) {
        std::string extname = "extno="+gammalib::str(extno);
        throw GException::fits_hdu_not_found(G_AT1, extname);
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to HDU
 *
 * @param[in] extname Name of HDU extension.
 * @return Pointer to HDU.
 *
 * @exception GException::fits_hdu_not_found
 *            No HDU with specified name has been found.
 *
 * Returns a pointer to the HDU with the specified @p extname. An exception
 * is thrown if the HDU is not valid (i.e. NULL).
 ***************************************************************************/
GFitsHDU* GFits::at(const std::string& extname)
{
    // Get extenion number
    int extno = this->extno(extname);
    if (extno == -1) {
        throw GException::fits_hdu_not_found(G_AT2, extname);
    }
    
    // Get HDU pointer
    GFitsHDU* ptr = m_hdu[extno];

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to HDU (const version)
 *
 * @param[in] extname Name of HDU extension.
 * @return Pointer to HDU.
 *
 * @exception GException::fits_hdu_not_found
 *            No HDU with specified name has been found.
 *
 * Returns a pointer to the HDU with the specified @p extname. An exception
 * is thrown if the HDU is not valid (i.e. NULL).
 ***************************************************************************/
const GFitsHDU* GFits::at(const std::string& extname) const
{
    // Get extenion number
    int extno = this->extno(extname);
    if (extno == -1) {
        throw GException::fits_hdu_not_found(G_AT2, extname);
    }
    
    // Get HDU pointer
    GFitsHDU* ptr = m_hdu[extno];

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to image HDU
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return Pointer to FITS image.
 *
 * @exception GException::fits_hdu_not_image
 *            Requested HDU is not an image.
 *
 * Returns a pointer to the image HDU with extension number extno.
 ***************************************************************************/
GFitsImage* GFits::image(const int& extno)
{
    // Get HDU pointer
    GFitsImage* ptr = dynamic_cast<GFitsImage*>(at(extno));

    // Throw an error if HDU is not an image
    if (ptr == NULL) {
        throw GException::fits_hdu_not_image(G_IMAGE1,
                                             "(extno="+gammalib::str(extno)+")",
                                             at(extno)->exttype());
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to image HDU (const version)
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return Pointer to FITS image.
 *
 * @exception GException::fits_hdu_not_image
 *            Requested HDU is not an image.
 *
 * Returns a pointer to the image HDU with extension number extno.
 ***************************************************************************/
const GFitsImage* GFits::image(const int& extno) const
{
    // Get HDU pointer
    const GFitsImage* ptr = dynamic_cast<const GFitsImage*>(at(extno));

    // Throw an error if HDU is not an image
    if (ptr == NULL) {
        throw GException::fits_hdu_not_image(G_IMAGE1,
                                             "(extno="+gammalib::str(extno)+")",
                                             at(extno)->exttype());
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to image HDU
 *
 * @param[in] extname Name of HDU extension.
 * @return Pointer to FITS image.
 *
 * @exception GException::fits_hdu_not_image
 *            Requested HDU is not an image.
 *
 * Returns a pointer to the image HDU with extension name extname.
 ***************************************************************************/
GFitsImage* GFits::image(const std::string& extname)
{
    // Get HDU pointer
    GFitsImage* ptr = dynamic_cast<GFitsImage*>(at(extname));

    // Throw an error if HDU is not an image
    if (ptr == NULL) {
        throw GException::fits_hdu_not_image(G_IMAGE2,
                                             extname,
                                             at(extname)->exttype());
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to image HDU (const version)
 *
 * @param[in] extname Name of HDU extension.
 * @return Pointer to FITS image.
 *
 * @exception GException::fits_hdu_not_image
 *            Requested HDU is not an image.
 *
 * Returns a pointer to the image HDU with extension name extname.
 ***************************************************************************/
const GFitsImage* GFits::image(const std::string& extname) const
{
    // Get HDU pointer
    const GFitsImage* ptr = dynamic_cast<const GFitsImage*>(at(extname));

    // Throw an error if HDU is not an image
    if (ptr == NULL) {
        throw GException::fits_hdu_not_image(G_IMAGE2,
                                             extname,
                                             at(extname)->exttype());
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to table HDU
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return Pointer to FITS table.
 *
 * @exception GException::fits_hdu_not_table
 *            Requested HDU is not a table.
 *
 * Returns a pointer to the table HDU with extension number extno.
 ***************************************************************************/
GFitsTable* GFits::table(const int& extno)
{
    // Get HDU pointer
    GFitsTable* ptr = dynamic_cast<GFitsTable*>(at(extno));

    // Throw an error if HDU is not an image
    if (ptr == NULL) {
        throw GException::fits_hdu_not_table(G_TABLE1,
                                             "(extno="+gammalib::str(extno)+")",
                                             at(extno)->exttype());
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to table HDU (const version)
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return Pointer to FITS table.
 *
 * @exception GException::fits_hdu_not_table
 *            Requested HDU is not a table.
 *
 * Returns a pointer to the table HDU with extension number extno.
 ***************************************************************************/
const GFitsTable* GFits::table(const int& extno) const
{
    // Get HDU pointer
    const GFitsTable* ptr = dynamic_cast<const GFitsTable*>(at(extno));

    // Throw an error if HDU is not an image
    if (ptr == NULL) {
        throw GException::fits_hdu_not_table(G_TABLE1,
                                             "(extno="+gammalib::str(extno)+")",
                                             at(extno)->exttype());
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to table HDU
 *
 * @param[in] extname Name of HDU extension.
 * @return Pointer to FITS table.
 *
 * @exception GException::fits_hdu_not_table
 *            Requested HDU is not a table.
 *
 * Returns a pointer to the table HDU with extension name extname.
 ***************************************************************************/
GFitsTable* GFits::table(const std::string& extname)
{
    // Get HDU pointer
    GFitsTable* ptr = dynamic_cast<GFitsTable*>(at(extname));

    // Throw an error if HDU is not an image
    if (ptr == NULL) {
        throw GException::fits_hdu_not_table(G_TABLE2,
                                             extname,
                                             at(extname)->exttype());
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to table HDU
 *
 * @param[in] extname Name of HDU extension.
 * @return Pointer to FITS table.
 *
 * @exception GException::fits_hdu_not_table
 *            Requested HDU is not a table.
 *
 * Returns a pointer to the table HDU with extension name extname.
 ***************************************************************************/
const GFitsTable* GFits::table(const std::string& extname) const
{
    // Get HDU pointer
    const GFitsTable* ptr = dynamic_cast<const GFitsTable*>(at(extname));

    // Throw an error if HDU is not an image
    if (ptr == NULL) {
        throw GException::fits_hdu_not_table(G_TABLE2,
                                             extname,
                                             at(extname)->exttype());
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Set HDU for the specified extension number
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @param[in] hdu HDU.
 * @return Pointer to cloned HDU.
 *
 * @exception GException::out_of_range
 *            Extension number out of range.
 * @exception GException::invalid_argument
 *            Attempt to insert non-image HDU in first slot.
 ***************************************************************************/
GFitsHDU* GFits::set(const int& extno, const GFitsHDU& hdu)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (extno < 0 || extno >= size()) {
        throw GException::out_of_range(G_SET1, "Extension number", extno, size());
    }
    #endif

    // Throw an exception if a non-image extension should be set in the
    // first slot
    if (extno == 0 && hdu.exttype() != GFitsHDU::HT_IMAGE) {
        std::string msg = "Attempt to set a table extension as the first"
                          " extension of a FITS file.\n"
                          "The first extension of a FITS file must be an"
                          " image, hence use the next slot to set a"
                          " table.";
        throw GException::invalid_argument(G_SET1, msg);
    }

    // Delete any existing HDU
    if (m_hdu[extno] != NULL) delete m_hdu[extno];

    // Assign new HDU by cloning
    m_hdu[extno] = hdu.clone();

    // If FITS file exists then connect cloned HDU to FITS file
    if (m_fitsfile != NULL) {
        __fitsfile fptr  = *(FPTR(m_fitsfile));
        fptr.HDUposition = extno;
        m_hdu[extno]->connect(&fptr);
    }
    else {
        m_hdu[extno]->extno(extno);
    }

    // Return pointer to HDU
    return (m_hdu[extno]);
}


/***********************************************************************//**
 * @brief Set HDU for the specified extension name
 *
 * @param[in] extname Name of HDU extension.
 * @param[in] hdu HDU.
 * @return Pointer to cloned HDU.
 *
 * @exception GException::invalid_argument
 *            Extension name not found.
 ***************************************************************************/
GFitsHDU* GFits::set(const std::string& extname, const GFitsHDU& hdu)
{
    // Get extenion number
    int extno = this->extno(extname);
    if (extno == -1) {
        std::string msg = "No extension with name \""+extname+"\" found in"
                          " FITS file.\n"
                          "Please specify a valid extension name.";
        throw GException::invalid_argument(G_SET2, msg);
    }

    // Set HDU and return pointer
    return (set(extno, hdu));
}


/***********************************************************************//**
 * @brief Append HDU to FITS file
 *
 * @param[in] hdu HDU.
 * @return Pointer to appended HDU.
 *
 * Append HDU to the next free position in a FITS file. In case that no HDU
 * exists so far in the FITS file and if the HDU to append is NOT an image,
 * an empty primary image will be inserted as first HDU in the FITS file.
 * This guarantees the compatibility with the FITS standard.
 ***************************************************************************/
GFitsHDU* GFits::append(const GFitsHDU& hdu)
{
    // Debug header
    #if defined(G_DEBUG)
    std::cout << "GFits::append(";
    switch (hdu.exttype()) {
    case GFitsHDU::HT_IMAGE:
        std::cout << "GFitsImage";
        break;
    case GFitsHDU::HT_ASCII_TABLE:
        std::cout << "GFitsAsciiTable";
        break;
    case GFitsHDU::HT_BIN_TABLE:
        std::cout << "GFitsBinTable";
        break;
    default:
        std::cout << "<unknown header type>";
        break;
    }
    std::cout << ") (size=" << size() << ") -->" << std::endl;
    #endif

    // Determine next free HDU number
    int n_hdu = size();

    // Add primary image if required
    if (n_hdu == 0 && hdu.exttype() != GFitsHDU::HT_IMAGE) {

        // Allocate primary image
        GFitsHDU* primary = new_primary();

        // If FITS file exists then connect primary image to FITS file
        if (m_fitsfile != NULL) {
            __fitsfile fptr  = *(FPTR(m_fitsfile));
            fptr.HDUposition = n_hdu;
            primary->connect(&fptr);
        }

        // Debug information
        #if defined(G_DEBUG)
        std::cout << "Append primary image to HDU." << std::endl;
        #endif

        // Push back primary image
        m_hdu.push_back(primary);

        // Increment HDU number
        n_hdu++;
    }

    // Clone FITS HDU
    GFitsHDU* ptr = hdu.clone();

    // Append HDU if it's valid
    if (ptr != NULL) {

        // If FITS file exists then connect cloned HDU to FITS file
        if (m_fitsfile != NULL) {
            __fitsfile fptr  = *(FPTR(m_fitsfile));
            fptr.HDUposition = n_hdu;
            ptr->connect(&fptr);
        }
        else {
            ptr->extno(n_hdu);
        }

        // Debug information
        #if defined(G_DEBUG)
        std::cout << "Append HDU (extno=" << n_hdu << ")." << std::endl;
        #endif

        // Push back HDU
        m_hdu.push_back(ptr);

        // Debug trailer
        #if defined(G_DEBUG)
        std::cout << "<-- GFits::append" << std::endl;
        #endif
        
    } // endif: HDU was valid

    // Return
    return ptr;
}


/***********************************************************************//**
 * @brief Set HDU for the specified extension number
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @param[in] hdu HDU.
 * @return Pointer to cloned HDU.
 *
 * @exception GException::out_of_range
 *            Extension number out of range.
 * @exception GException::invalid_argument
 *            Attempt to insert non-image HDU in first slot.
 ***************************************************************************/
GFitsHDU* GFits::insert(const int& extno, const GFitsHDU& hdu)
{
    // Compile option: raise exception if extension number is out of range
    #if defined(G_RANGE_CHECK)
    if (extno < 0 || extno >= size()) {
        throw GException::out_of_range(G_INSERT1, "Extension number", extno, size());
    }
    #endif

    // Throw an exception if a non-image extension should be inserted in the
    // first slot
    if (extno == 0 && hdu.exttype() != GFitsHDU::HT_IMAGE) {
        std::string msg = "Attempt to insert a table extension as the first"
                          " extension of a FITS file.\n"
                          "The first extension of a FITS file must be an"
                          " image, hence use the next slot to insert a"
                          " table.";
        throw GException::invalid_argument(G_INSERT1, msg);
    }

    // Create deep copy of HDU
    GFitsHDU* ptr = hdu.clone();

    // Inserts deep copy of HDU
    m_hdu.insert(m_hdu.begin() + extno, ptr);

    // If FITS file exists then connect cloned HDU to FITS file
    if (m_fitsfile != NULL) {
        __fitsfile fptr  = *(FPTR(m_fitsfile));
        fptr.HDUposition = extno;
        ptr->connect(&fptr);
    }
    else {
        ptr->extno(extno);
    }

    // Update extno for all subsequent HDUs
    for (int i = extno + 1; i < size(); ++i) {
        m_hdu[i]->extno(i);
    }

    // Return pointer to HDU
    return ptr;
}


/***********************************************************************//**
 * @brief Insert HDU into FITS file
 *
 * @param[in] extname Extension name.
 * @param[in] hdu Extension.
 * @return Pointer to deep copy of extension @p hdu.
 *
 * @exception GException::invalid_argument
 *            Extension name not found.
 *
 * Inserts an extension @p hdu into the FITS file before the extension with
 * the specified @p extname.
 ***************************************************************************/
GFitsHDU* GFits::insert(const std::string& extname, const GFitsHDU& hdu)
{
    // Get extenion number
    int extno = this->extno(extname);
    if (extno == -1) {
        std::string msg = "No extension with name \""+extname+"\" found in"
                          " FITS file.\n"
                          "Please specify a valid extension name.";
        throw GException::invalid_argument(G_INSERT2, msg);
    }

    // Insert HDU and return pointer
    return (insert(extno, hdu));
}


/***********************************************************************//**
 * @brief Remove HDU from FITS file
 *
 * @param[in] extno Extension number [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Specified @p extno is out of range.
 *
 * @todo Handle HDU update in FITS file.
 ***************************************************************************/
void GFits::remove(const int& extno)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (extno < 0 || extno >= size()) {
        throw GException::out_of_range(G_REMOVE1, "Extension number", extno, size());
    }
    #endif

    // Throw an exception if removal would lead to a non-image extension
    // in the first slot
    if (extno == 0 && size() > 1 && m_hdu[1]->exttype() != GFitsHDU::HT_IMAGE) {
        std::string msg = "Attempt to remove primary image from a FITS file"
                          " with a table extension as second extension.\n"
                          "The removal of the primary image would result in"
                          " having a table as the first extension of the FITS"
                          " file, which is not a valid FITS file.";
        throw GException::invalid_argument(G_REMOVE1, msg);
    }

    // Erase HDU from FITS file
    m_hdu.erase(m_hdu.begin() + extno);

    // Update extno for all subsequent HDUs
    for (int i = extno; i < size(); ++i) {
        m_hdu[i]->extno(i);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove HDU from FITS file
 *
 * @param[in] extname Name of HDU extension.
 *
 * @exception GException::invalid_argument
 *            Specified @p extname not found in FITS file.
 ***************************************************************************/
void GFits::remove(const std::string& extname)
{
    // Get extenion number
    int extno = this->extno(extname);
    if (extno == -1) {
        std::string msg = "No extension with name \""+extname+"\" found in"
                          " FITS file.\n"
                          "Please specify a valid extension name.";
        throw GException::invalid_argument(G_REMOVE2, msg);
    }

    // Remove extension
    remove(extno);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append FITS file
 *
 * @param[in] fits FITS file.
 *
 * Append all extension of @p fits file to FITS file.
 *
 * @todo Handle HDU update in FITS file.
 ***************************************************************************/
void GFits::extend(const GFits& fits)
{
    // Do nothing if FITS file is empty
    if (!fits.isempty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = fits.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all HDUs and append pointers to deep copies 
        for (int i = 0; i < num; ++i) {

            // Clone HDU
            GFitsHDU* ptr = m_hdu[i]->clone();

            // Push back HDU on stack
            m_hdu.push_back(ptr);

            // Retrieve extno
            int extno = m_hdu.size()-1;

            // If FITS file exists then connect cloned HDU to FITS file
            if (m_fitsfile != NULL) {
                __fitsfile fptr  = *(FPTR(m_fitsfile));
                fptr.HDUposition = extno;
                ptr->connect(&fptr);
            }
            else {
                ptr->extno(extno);
            }

        } // endfor: looped over all HDUs

    } // endif: FITS file was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return extension number for specified extension name
 *
 * @param[in] extname Name of HDU extension.
 * @return Extension number (-1 of not found).
 *
 * Returns the extension number for a specified extension name @p extname. If
 * the extension name if "PRIMARY" an extension number of 0 is returned.
 * If the extension is not found, -1 is returned.
 ***************************************************************************/
int GFits::extno(const std::string& extname) const
{
    // Initialise extension number
    int extno = -1;

    // Return primary HDU if requested ...
    if (gammalib::toupper(extname) == "PRIMARY") {
        if (size() > 0) {
            extno = 0;
        }
    }

    // ... otherwise search for specified extension
    else {
        for (int i = 0; i < size(); ++i) {
            if (m_hdu[i]->extname() == extname) {
                extno = i;
                break;
            }
        }
    }

    // Return extno
    return extno;
}


/***********************************************************************//**
 * @brief Open or (optionally) create FITS file
 *
 * @param[in] filename Name of FITS file to be opened.
 * @param[in] create Create FITS file if it does not exist (default=false)
 *
 * @exception GException::fits_already_opened
 *            Class instance contains already an opened FITS file.
 *            Close file before opening a new one using GFits::close().
 * @exception GException::fits_open_error 
 *            Unable to open the specified file.
 * @exception GException::fits_hdu_not_found
 *            Requested HDU not found.
 * @exception GException::fits_error
 *            Unable to determine number of HDUs in the FITS file.
 *
 * This method opens all HDUs that are found in the specified FITS file.
 * If the file does not exist, and if create=true, a new FITS file is created.
 * For each HDU, a GFitsHDU object is associated to the GFits object.
 * The HDUs can then be accessed using the hdu(const std::string&) or
 * hdu(int extno) method.
 * Any environment variable present in the filename will be expanded.
 ***************************************************************************/
void GFits::open(const std::string& filename, const bool& create)
{
    // Remove any HDUs
    m_hdu.clear();

    // Don't allow opening if another file is already open
    if (m_fitsfile != NULL) {
        throw GException::fits_already_opened(G_OPEN, m_filename);
    }

    // Expand environment variables
    std::string fname = gammalib::expand_env(filename);

    // Initialise FITS file as readwrite and non created
    m_readwrite = true;
    m_created   = false;

    // Try opening FITS file with readwrite access
    int status = 0;
    status     = __ffopen(FHANDLE(m_fitsfile), fname.c_str(), 1, &status);

    // If failed then try opening as readonly
    if (status == 104 || status == 112) {
        status      = 0;
        status      = __ffopen(FHANDLE(m_fitsfile), fname.c_str(), 0, &status);
        m_readwrite = false;
    }

    // If failed and if we are allowed to create a new FITS file then create
    // FITS file now
    if (create && status == 104) {
        status      = 0;
        status      = __ffinit(FHANDLE(m_fitsfile), fname.c_str(), &status);
        m_readwrite = true;
        m_created   = true;
    }

	// Throw special exception if status=202 (keyword not found). This error
	// may occur if the file is opened with an expression
	if (status == 202) {
        throw GException::fits_open_error(G_OPEN, fname, status,
		                  "Keyword not found when opening file.");
    }

    // Throw any other error
    else if (status != 0) {
        throw GException::fits_open_error(G_OPEN, fname, status);
    }

    // Store FITS file attributes
    m_filename = fname;

    // Determine number of HDUs
    int num_hdu = 0;
    status = __ffthdu(FPTR(m_fitsfile), &num_hdu, &status);
    if (status != 0) {
        throw GException::fits_error(G_OPEN, status);
    }

    // Open and append all HDUs
    for (int i = 0; i < num_hdu; ++i) {

        // Move to HDU and get HDU type
        int type = gammalib::fits_move_to_hdu(G_OPEN, m_fitsfile, i+1);

        // Perform type dependent HDU allocation
        GFitsHDU* hdu = NULL;
        switch (type) {
        case GFitsHDU::HT_IMAGE:
            hdu = new_image();
            break;
        case GFitsHDU::HT_ASCII_TABLE:
            hdu = new GFitsAsciiTable;
            break;
        case GFitsHDU::HT_BIN_TABLE:
            hdu = new GFitsBinTable;
            break;
        default:
            std::string msg = "Unknown HDU type \""+gammalib::str(type)+"\"";
            throw GException::fits_invalid_type(G_OPEN, msg);
            break;
        }

        // Open HDU
        hdu->open(FPTR(m_fitsfile), i);

        // Append HDU
        m_hdu.push_back(hdu);

    } // endfor: looped over all HDUs

    // Return
    return;
}


/***********************************************************************//**
 * @brief Saves FITS file
 *
 * @param[in] clobber Overwrite existing FITS file (true/false).
 *
 * @exception GException::fits_file_exist
 *            Attemting to overwrite an existing file without having specified
 *            clobber=true.
 * @exception GException::fits_file_not_open
 *            FITS file needs to be opened before saving.
 *
 * Saves all HDUs of an open FITS file to disk. After saving, the FITS file
 * remains open. Invoke the close() method if explicit closing is needed.
 * Note that de-allocation of the GFits object also closes the FITS file.
 *
 * In the special case that no first HDU exists an empty primary image is
 * created.
 ***************************************************************************/
void GFits::save(const bool& clobber)
{
    // Debug header
    #if defined(G_DEBUG)
    std::cout << "GFits::save (size=" << size() << ") -->" << std::endl;
    #endif

    // Initialise cfitsio status
    int status = 0;

    // If we attempt to save an existing file without overwriting permission
    // then throw an error
    if (!m_created && !clobber) {
        throw GException::fits_file_exist(G_SAVE, m_filename);
    }

    // If no FITS file has been opened then throw an error
    if (m_fitsfile == NULL) {
        throw GException::fits_file_not_open(G_SAVE, 
              "Either open FITS file before saving or use saveto() method.");
    }

    // If no HDUs exist then save an empty primary image.
    if (size() == 0) {
        status     = __ffmahd(FPTR(m_fitsfile), 1, NULL, &status);
        status     = __ffcrim(FPTR(m_fitsfile), 8, 0, NULL, &status);
        if (status != 0) {
            throw GException::fits_error(G_SAVE, status);
        }
    }

    // ... otherwise save all HDUs
    else {
        for (int i = 0; i < size(); ++i) {
            m_hdu[i]->extno(i);
            m_hdu[i]->save();
        }
    }

    // Determine number of HDUs in FITS file
    int num_hdu = 0;
    status = __ffthdu(FPTR(m_fitsfile), &num_hdu, &status);
    if (status != 0) {
        throw GException::fits_error(G_SAVE, status);
    }

    // Delete all excedent HDUs (these may be their following removal
    // of HDUs from the GFits object)
    for (int i = num_hdu-1; i >= size(); --i) {
        status = __ffdhdu(FPTR(m_fitsfile), NULL, &status);
        if (status != 0) {
            throw GException::fits_error(G_SAVE, status);
        }
    }

    // Signal that file needs not to be created anymore
    m_created = false;
    
    // Debug trailer
    #if defined(G_DEBUG)
    std::cout << "<-- GFits::save" << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Saves to specified FITS file
 *
 * @param[in] filename Filename.
 * @param[in] clobber Overwrite existing FITS file (true/false).
 *
 * @exception GException::fits_file_exist
 *            Specified file exists already. Overwriting requires
 *            clobber=true.
 *
 * Saves object into a specific FITS file.
 * Any environment variable present in the filename will be expanded.
 ***************************************************************************/
void GFits::saveto(const std::string& filename, const bool& clobber)
{
    // Expand environment variables
    std::string fname = gammalib::expand_env(filename);

    // Debug header
    #if defined(G_DEBUG)
    std::cout << "GFits::saveto(\"" << fname << "\", " << clobber << ")"
              << " (size=" << size() << ") -->" << std::endl;
    #endif

    // If overwriting has been specified then remove any existing file ...
    if (clobber) {
        std::remove(fname.c_str());
    }

    // ... otherwise, if file exists then throw an exception
    else if (gammalib::file_exists(fname)) {
        throw GException::fits_file_exist(G_SAVETO, fname);
    }

    // Create or open FITS file
    GFits new_fits;
    new_fits.open(fname, true);

    // Append all headers
    for (int i = 0; i < size(); ++i) {
        new_fits.append(*m_hdu[i]);
    }

    // Save new FITS file
    new_fits.save();

    // Close new FITS file
    new_fits.close();

    // Debug trailer
    #if defined(G_DEBUG)
    std::cout << "<-- GFits::saveto" << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close FITS file
 *
 * Closing detaches a FITS file from the GFits object and returns a clean
 * empty object.
 ***************************************************************************/
void GFits::close(void)
{
    // Close file and free all members 
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print FITS information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing FITS information.
 ***************************************************************************/
std::string GFits::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GFits ===");

        // Append file information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("History"));
        if (m_created) {
            result.append("new file");
        }
        else {
            result.append("existing file");
        }
        result.append("\n"+gammalib::parformat("Mode"));
        if (m_readwrite) {
            result.append("read/write");
        }
        else {
            result.append("read only");
        }
        result.append("\n"+gammalib::parformat("Number of HDUs"));
        result.append(gammalib::str(size()));

        // NORMAL: Append HDUs
        if (chatter >= NORMAL) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n");
                result.append(m_hdu[i]->print(gammalib::reduce(chatter)));
            }
        }

    } // endif: chatter was not silent

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
 ***************************************************************************/
void GFits::init_members(void)
{
    // Initialise GFits members
    m_hdu.clear();
    m_filename.clear();
    m_fitsfile  = NULL;
    m_readwrite = true;
    m_created   = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief  Copy class members
 *
 * @param fits Object to be copied
 *
 * The method does not copy the FITS filename and FITS file pointer.
 * This prevents that several copies of the FITS file pointer exist in
 * different instances of GFits, which would lead to confusion since one
 * instance could close the file while for another it still would be
 * opened. The rule ONE INSTANCE - ONE FILE applies.
  ***************************************************************************/
void GFits::copy_members(const GFits& fits)
{
    // Reset FITS file attributes
    m_filename.clear();
    m_fitsfile  = NULL;
    m_readwrite = true;
    m_created   = true;

    // Clone HDUs
    m_hdu.clear();
    for (int i = 0; i < fits.m_hdu.size(); ++i) {
        m_hdu.push_back((fits.m_hdu[i]->clone()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * Closes the FITS file if it had been opened and deallocate all HDUs.
 *
 * If the G_DELETE_EMPTY_FITS_FILES option is defined, files without HDUs
 * or corrupted files will be deleted. This prevents leaving corrupted files
 * on disk (yet, corrupted files may be generated by another application,
 * thus this is not 100% safe; better make the code solid against reading
 * corrupted FITS files).
 ***************************************************************************/
void GFits::free_members(void)
{
    // If FITS file has been opened then close it now
    if (m_fitsfile != NULL) {

        // If file has been created but not yet save then delete the file
        // now. We do not worry about the status in this case.
        if (m_created) {
            int status = 0;
            __ffdelt(FPTR(m_fitsfile), &status);
        }

        // Compile option: If there are no HDUs then delete the file (don't
        // worry about error)
        #if defined(G_DELETE_EMPTY_FITS_FILES)
        else if (size() == 0) {
            int status = 0;
            __ffdelt(FPTR(m_fitsfile), &status);
        }
        #endif

        // ... otherwise close the file
        else {
            int status = 0;
            status     = __ffclos(FPTR(m_fitsfile), &status);
            if (status == 252) {
                int new_status = 0;
                __ffdelt(FPTR(m_fitsfile), &new_status);
                if (new_status != 0) {
                    throw GException::fits_error(G_FREE_MEM, new_status);
                }
            }
            else if (status != 0) {
                throw GException::fits_error(G_FREE_MEM, status);
            }

        } // endelse: there was an open FITS file

    } // endif: FITS file had been opened

    // Free HDUs
    for (int i = 0; i < m_hdu.size(); ++i) {
        if (m_hdu[i] != NULL) {
            delete m_hdu[i];
            m_hdu[i] = NULL;
        }
    }

    // Clear HDUs
    m_hdu.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate new FITS image and return memory pointer
 *
 * Depending on the number of bits per pixel, a FITS image is allocated
 * and the pointer is returned. The following FITS image classes are
 * handled:
 * GFitsImageByte     (bitpix=8)
 * GFitsImageShort    (bitpix=16)
 * GFitsImageLong     (bitpix=32)
 * GFitsImageLongLong (bitpix=64)
 * GFitsImageFloat    (bitpix=-32)
 * GFitsImageDouble   (bitpix=-64)
 * The information about the number of bits per pixels is extracted from
 * the actual HDU.
 *
 * @todo Additional code is needed to detect unsigned integer images. This
 * code may be insprired by the code used for table columns as the unsigned
 * information is stored in the BZERO keyword.
 ***************************************************************************/
GFitsImage* GFits::new_image(void)
{
    // Initialise return value
    GFitsImage* image = NULL;

    // Get number of bits per pixel
    int status =   0;
    int bitpix = -64;
    //status     = __ffgipr(FPTR(m_fitsfile), 0, &bitpix, NULL, NULL, &status);
    //status     = __ffgidt(FPTR(m_fitsfile), &bitpix, &status);
    status     = __ffgiet(FPTR(m_fitsfile), &bitpix, &status);
    if (status != 0) {
        throw GException::fits_error(G_NEW_IMAGE, status);
    }

    // Check for unsigned image
    char      keyname[10];
    long long bzero;
    long long bscale;
    std::sprintf(keyname, "BZERO");
    if (__ffgky(FPTR(m_fitsfile), __TLONGLONG, keyname, &bzero, NULL, &status) != 0) {
        bzero = 0;
    }
    std::sprintf(keyname, "BSCALE");
    if (__ffgky(FPTR(m_fitsfile), __TLONGLONG, keyname, &bscale, NULL, &status) != 0) {
        bscale = 0;
    }
    if (bitpix == 8 && bzero == -128 && bscale == 1) {
        bitpix = 10;
    }
    else if (bitpix == 16 && bzero == 32768u && bscale == 1) {
        bitpix = 20;
    }
    else if (bitpix == 32 && bzero == 2147483648u && bscale == 1) {
        bitpix = 40;
    }

    // Allocate bitpix dependent image
    switch (bitpix) {
    case 8:
        image = new GFitsImageByte;
        break;
    case 10:
        image = new GFitsImageSByte;
        break;
    case 16:
        image = new GFitsImageShort;
        break;
    case 20:
        image = new GFitsImageUShort;
        break;
    case 32:
        image = new GFitsImageLong;
        break;
    case 40:
        image = new GFitsImageULong;
        break;
    case 64:
        image = new GFitsImageLongLong;
        break;
    case -32:
        image = new GFitsImageFloat;
        break;
    case -64:
        image = new GFitsImageDouble;
        break;
    default:
        throw GException::fits_bad_bitpix(G_NEW_IMAGE, bitpix);
        break;
    }

    // Return image pointer
    return image;
}


/***********************************************************************//**
 * @brief Return minimal primary HDU
 *
 * Creates a primary HDU in memory and open it using the GFitsHDU::open()
 * method.
 ***************************************************************************/
GFitsImage* GFits::new_primary(void)
{
    // Allocate an empty image
    GFitsImage* image = new GFitsImageByte;

    // Create primary image in memory
    int status = 0;
    __fitsfile* fptr;
    status = __ffinit(&fptr, "mem://", &status);
    status = __ffcrim(fptr, 8, 0, NULL, &status);

    // Open HDU
    image->open(fptr,0);

    // Close FITS file in memory
    status = __ffclos(fptr, &status);
    if (status == 252) {
        status = 0;
        status = __ffdelt(fptr, &status);
    }

    // Initialise FITS file pointer
    FPTR(image->m_fitsfile)->HDUposition = 0;
    FPTR(image->m_fitsfile)->Fptr        = NULL;

    // Return image
    return image;
}


/*==========================================================================
 =                                                                         =
 =                         FITS utility functions                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Move to FITS extension
 *
 * @param[in] caller Name of caller.
 * @param[in] vptr FITS file void pointer.
 * @param[in] hdunum HDU number (optional)
 *
 * @exception GException::fits_error
 *            cfitsio error occured.
 *
 * If @p hdunum is >0, moves the FITS file void pointer to the HDU
 * specified by @p hdunum. Otherwise, the FITS file void pointer is moved
 * to the HDU specified by the @p HDUposition attribute of the void pointer.
 ***************************************************************************/
int gammalib::fits_move_to_hdu(const std::string& caller, void* vptr,
                               const int& hdunum)
{
    // Initialise status and HDU type
    int status = 0;
    int type   = 0;

    // Set HDU position
    int position = (hdunum > 0) ? hdunum : (FPTR(vptr)->HDUposition)+1;

    // Move to HDU
    status = __ffmahd(FPTR(vptr), position, &type, &status);

    // Throw exception in case of an error
    if (status != 0) {
        std::string msg = "Unable to move FITS file pointer to extension"
                          " number "+
                          gammalib::str(position-1)+".";
        throw GException::fits_error(caller, status, msg);
    }

    // Return HDU type
    return type;
}
