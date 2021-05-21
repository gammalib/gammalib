/***************************************************************************
 *               GFitsImage.cpp - Abstract FITS image base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2021 by Juergen Knoedlseder                         *
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
 * @file GFitsImage.cpp
 * @brief Abstract FITS image base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GFitsCfitsio.hpp"
#include "GFits.hpp"
#include "GFitsImage.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAXES                                      "GFitsImage::naxes(int)"
#define G_OPEN_IMAGE                                "GFitsImage::open(void*)"
#define G_LOAD_IMAGE           "GFitsImage::load_image(int,void*,void*,int*)"
#define G_SAVE_IMAGE                      "GFitsImage::save_image(int,void*)"
#define G_OFFSET_1D                                "GFitsImage::offset(int&)"
#define G_OFFSET_2D                           "GFitsImage::offset(int&,int&)"
#define G_OFFSET_3D                      "GFitsImage::offset(int&,int&,int&)"
#define G_OFFSET_4D                 "GFitsImage::offset(int&,int&,int&,int&)"

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
 *
 * Construct instance of an empty image. No header cards are present in an
 * empty image.
 ***************************************************************************/
GFitsImage::GFitsImage(void) : GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief 1D image constructor
 *
 * @param[in] bitpix Number of Bits per pixel (negative is floating point).
 * @param[in] nx Number of pixels.
 *
 * Construct 1D instance of GFitsImage by specifying the number of pixels.
 * This method also adds the relevant header cards.
 ***************************************************************************/
GFitsImage::GFitsImage(const int& bitpix, const int& nx) : GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Store number of Bits per pixel
    m_bitpix = bitpix;

    // Set number of axes
    m_naxis = 1;

    // Set image dimensions
    m_naxes      = new long[m_naxis];
    m_naxes[0]   = nx;
    m_num_pixels = nx;

    // Initialise header
    init_image_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief 2D image constructor
 *
 * @param[in] bitpix Number of Bits per pixel (negative is floating point).
 * @param[in] nx Number of pixels in first dimension.
 * @param[in] ny Number of pixels in second dimension.
 *
 * Construct 2D instance of GFitsImage by specifying the number of pixels
 * in each dimension. This method also adds the relevant header cards.
 ***************************************************************************/
GFitsImage::GFitsImage(const int& bitpix, const int& nx, const int& ny) :
            GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Store number of Bits per pixel
    m_bitpix = bitpix;

    // Set number of axes
    m_naxis = 2;

    // Set image dimensions
    m_naxes      = new long[m_naxis];
    m_naxes[0]   = nx;
    m_naxes[1]   = ny;
    m_num_pixels = nx * ny;

    // Initialise header
    init_image_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief 3D image constructor
 *
 * @param[in] bitpix Number of Bits per pixel (negative is floating point).
 * @param[in] nx Number of pixels in first dimension.
 * @param[in] ny Number of pixels in second dimension.
 * @param[in] nz Number of pixels in third dimension.
 *
 * Construct 3D instance of GFitsImage by specifying the number of pixels
 * in each dimension. This method also adds the relevant header cards.
 ***************************************************************************/
GFitsImage::GFitsImage(const int& bitpix, const int& nx, const int& ny,
                       const int& nz) : GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Store number of Bits per pixel
    m_bitpix = bitpix;

    // Set number of axes
    m_naxis = 3;

    // Set image dimensions
    m_naxes      = new long[m_naxis];
    m_naxes[0]   = nx;
    m_naxes[1]   = ny;
    m_naxes[2]   = nz;
    m_num_pixels = nx * ny * nz;

    // Initialise header
    init_image_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief 4D image constructor
 *
 * @param[in] bitpix Number of Bits per pixel (negative is floating point).
 * @param[in] nx Number of pixels in first dimension.
 * @param[in] ny Number of pixels in second dimension.
 * @param[in] nz Number of pixels in third dimension.
 * @param[in] nt Number of pixels in forth dimension.
 *
 * Construct 4D instance of GFitsImage by specifying the number of pixels
 * in each dimension. This method also adds the relevant header cards.
 ***************************************************************************/
GFitsImage::GFitsImage(const int& bitpix, const int& nx, const int& ny,
                       const int& nz, const int& nt) : GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Store number of Bits per pixel
    m_bitpix = bitpix;

    // Set number of axes
    m_naxis = 4;

    // Set image dimensions
    m_naxes      = new long[m_naxis];
    m_naxes[0]   = nx;
    m_naxes[1]   = ny;
    m_naxes[2]   = nz;
    m_naxes[3]   = nt;
    m_num_pixels = nx * ny * nz * nt;

    // Initialise header
    init_image_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] bitpix Number of Bits per pixel (negative is floating point).
 * @param[in] naxes Vector of number of pixels in each dimension.
 *
 * Construct instance of GFitsImage by specifying the image dimension and
 * the number of pixels in each dimension. This method also adds the relevant
 * header cards.
 ***************************************************************************/
GFitsImage::GFitsImage(const int& bitpix, const std::vector<int>& naxes) : GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Store number of Bits per pixel
    m_bitpix = bitpix;

    // Store number of axes
    m_naxis = naxes.size();

    // Copy number of pixels in each dimension and calculate the total
    // number of pixels
    if (m_naxis > 0) {
        m_naxes      = new long[m_naxis];
        m_num_pixels = 1;
        for (int i = 0; i < m_naxis; ++i) {
            m_naxes[i]    = naxes[i];
            m_num_pixels *= naxes[i];
        }
    }

    // Initialise header
    init_image_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] image FITS image.
 ***************************************************************************/
GFitsImage::GFitsImage(const GFitsImage& image) : GFitsHDU(image)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(image);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFitsImage::~GFitsImage(void)
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
 * @param[in] image FITS image.
 * @return FITS image.
 ***************************************************************************/
GFitsImage& GFitsImage::operator=(const GFitsImage& image)
{
    // Execute only if object is not identical
    if (this != &image) {

        // Copy base class members
        this->GFitsHDU::operator=(image);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(image);

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
 * @brief Return dimension of an image axis
 *
 * @param[in] axis Image axis [0,...,naxis()-1].
 *
 * @exception GException::out_of_range
 *            Image axis not valid.
 ***************************************************************************/
int GFitsImage::naxes(const int& axis) const
{
    // Check if axis is within the range
    #if defined(G_RANGE_CHECK)
    if (axis < 0 || axis >= naxis()) {
        throw GException::out_of_range(G_NAXES, "Image axis", axis, naxis());
    }
    #endif

    // Get axis dimension
    int dim = m_naxes[axis];

    // Return axis dimension
    return dim;
}


/***********************************************************************//**
 * @brief Set nul value
 *
 * @param[in] value Nul value.
 *
 * @todo To correctly reflect the nul value in the data, the image should
 * be reloaded. However, the image may have been changed, so in principle
 * saving is needed. However, we may not want to store the image, hence saving
 * is also not desired. We thus have to develop a method to update the
 * image information for a new nul value in place ...
 ***************************************************************************/
void GFitsImage::nulval(const void* value)
{
    // Allocate nul value
    alloc_nulval(value);

    // Update image
    //TODO

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print column information
 *
 * @param[in] chatter Chattiness.
 * @return String containing column information.
 *
 * @todo Format and cfitsio information is mainly for debugging. This could
 * be vanish in a more stable version of the code, or it could be compiled
 * in conditionally using a debug option.
 ***************************************************************************/
std::string GFitsImage::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GFitsImage ===");

        // Append HDU information
        result.append("\n"+print_hdu(chatter));

        // Append image dimensions
        result.append("\n"+gammalib::parformat("Image type"));
        result.append(typecode(type()));
        result.append("\n"+gammalib::parformat("Number of dimensions"));
        result.append(gammalib::str(naxis()));
        result.append("\n"+gammalib::parformat("Number of image pixels"));
        result.append(gammalib::str(npix()));
        for (int i = 0; i < naxis(); ++i) {
            result.append("\n"+gammalib::parformat("Number of bins in "+gammalib::str(i)));
            result.append(gammalib::str(naxes(i)));
        }

        // NORMAL: Append header information
        if (chatter >= NORMAL) {
            result.append(+"\n"+m_header.print(gammalib::reduce(chatter)));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsImage::init_members(void)
{
    // Initialise members
    m_bitpix     = 8;
    m_naxis      = 0;
    m_naxes      = NULL;
    m_num_pixels = 0;
    m_anynul     = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] image FITS image to copy
 ***************************************************************************/
void GFitsImage::copy_members(const GFitsImage& image)
{
    // Copy attributes
    m_bitpix     = image.m_bitpix;
    m_naxis      = image.m_naxis;
    m_num_pixels = image.m_num_pixels;
    m_anynul     = image.m_anynul;

    // Copy axes
    m_naxes = NULL;
    if (image.m_naxes != NULL && m_naxis > 0) {
        m_naxes = new long[m_naxis];
        for (int i = 0; i < m_naxis; ++i) {
            m_naxes[i] = image.m_naxes[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsImage::free_members(void)
{
    // Free memory
    if (m_naxes != NULL) delete [] m_naxes;

    // Mark memory as free
    m_naxes = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Open FITS image
 *
 * @param[in] vptr FITS file pointer
 *
 * Open FITS image in FITS file. Opening means connecting the FITS file
 * pointer to the image and reading the image and axes dimensions.
 ***************************************************************************/
void GFitsImage::data_open(void* vptr)
{
    // Open image
    open_image(vptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save FITS image
 *
 * Saves the image into the FITS file.
 ***************************************************************************/
void GFitsImage::data_save(void)
{
    // Save image
    save_image(type(), pixels());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close FITS image
 *
 * Closing a FITS image resets the object into its initial state. Closing
 * does NOT save the image into the FITS file. Use the save method for this
 * purpose.
 *
 * @todo Not sure that this is efficient at this level since the pixel array
 * will not be deallocated!!!
 ***************************************************************************/
void GFitsImage::data_close(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Connect FITS image
 *
 * @param[in] vptr FITS file pointer
 *
 * Connects a FITS file pointer to an image. This method actually does
 * nothing since any GFitsImage can directly access the FITS file pointer
 * that is stored in the GFitsHDU base class.
 ***************************************************************************/
void GFitsImage::data_connect(void* vptr)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise image header
 *
 * Initialises the image header by setting the default header cards. This
 * method requires the members m_bitpix, m_naxis, and m_naxes to be set
 * previously.
 ***************************************************************************/
void GFitsImage::init_image_header(void)
{
    // Set image header keywords
    m_header.append(GFitsHeaderCard("XTENSION", "IMAGE   ",
                                    "IMAGE extension"));
    m_header.append(GFitsHeaderCard("BITPIX", bitpix(),
                                    "number of bits per data pixel"));
    m_header.append(GFitsHeaderCard("NAXIS", naxis(),
                                    "number of data axes"));
    for (int i = 0; i < naxis(); ++i) {
        std::ostringstream s_key;
        std::ostringstream s_comment;
        s_key     << "NAXIS" << (i+1);
        s_comment << "length of data axis " << (i+1);
        m_header.append(GFitsHeaderCard(s_key.str(), naxes(i),
                                        s_comment.str()));
    }
    m_header.append(GFitsHeaderCard("PCOUNT", 0,
                                    "required keyword; must = 0"));
    m_header.append(GFitsHeaderCard("GCOUNT", 1,
                                    "required keyword; must = 1"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Open Image
 *
 * @param[in] vptr FITS file void pointer.
 *
 * @exception GException::fits_error
 *            FITS error.
 *
 * Open FITS image in FITS file. Opening means connecting the FITS file
 * pointer to the image and reading the image and axes dimensions.
 ***************************************************************************/
void GFitsImage::open_image(void* vptr)
{
    // Move to HDU
    gammalib::fits_move_to_hdu(G_OPEN_IMAGE, vptr);

    // Save the FITS file pointer and the HDU number
    FPTR_COPY(m_fitsfile, vptr);
    m_hdunum = FPTR(vptr)->HDUposition;

    // Get the image dimensions
    int status = 0;
    status = __ffgidm(FPTR(m_fitsfile), &m_naxis, &status);
    if (status != 0) {
        throw GException::fits_error(G_OPEN_IMAGE, status);
    }

    // Reset number of image pixels
    m_num_pixels = 0;

    // Get the axes dimensions
    if (m_naxis > 0) {

        // Allocate memory for axes dimensions
        if (m_naxes != NULL) delete [] m_naxes;
        m_naxes      = new long[m_naxis];

        // Get the axes dimensions
        status = __ffgisz(FPTR(m_fitsfile), m_naxis, m_naxes, &status);
        if (status != 0) {
            throw GException::fits_error(G_OPEN_IMAGE, status);
        }

        // Calculate number of image pixels
        m_num_pixels = 1;
        for (int i = 0; i < m_naxis; ++i) {
            m_num_pixels *= m_naxes[i];
        }

    } // endif: there is an image

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load FITS image
 *
 * @param[in] datatype Datatype of pixels to be saved.
 * @param[in] pixels Pixel array to be saved.
 * @param[in] nulval Pointer to pixel nul value.
 * @param[out] anynul Number of nul values encountered during loading.
 *
 * @exception GException::fits_error
 *            FITS error.
 *
 * Load image pixels from FITS file.
 ***************************************************************************/
void GFitsImage::load_image(int datatype, const void* pixels,
                            const void* nulval, int* anynul)
{
    // Move to HDU
    move_to_hdu();

    // Load the image pixels (if there are some ...)
    if (m_naxis > 0) {
        long* fpixel = new long[m_naxis];
        long* lpixel = new long[m_naxis];
        long* inc    = new long[m_naxis];
        for (int i = 0; i < m_naxis; ++i) {
            fpixel[i] = 1;
            lpixel[i] = m_naxes[i];
            inc[i]    = 1;
        }
        int status = 0;
        status     = __ffgsv(FPTR(m_fitsfile), datatype, fpixel, lpixel, inc,
                             (void*)nulval, (void*)pixels, anynul, &status);
        delete [] fpixel;
        delete [] lpixel;
        delete [] inc;
        if (status != 0) {
            throw GException::fits_error(G_LOAD_IMAGE, status);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save FITS image
 *
 * @param[in] datatype Datatype of pixels to be saved
 * @param[in] pixels Pixel array to be saved
 *
 * @exception GException::runtime_error
 *            No FITS file has been opened.
 * @exception GException::fits_error
 *            FITS error.
 *
 * Save image pixels into FITS file. In case that the HDU does not exist it
 * is created. In case that the pixel array is empty no data are saved; all
 * image pixels will be empty in this case.
 ***************************************************************************/
void GFitsImage::save_image(int datatype, const void* pixels)
{
    // Throw an exception if FITS file is not open
    if (FPTR(m_fitsfile)->Fptr == NULL) {
        std::string msg = "FITS file not open. Please open the FITS file "
                          "before saving the image.";
        throw GException::runtime_error(G_SAVE_IMAGE, msg);
    }

    // Move to HDU. We use here an explicit cfitsio moveto function since we
    // want to recover the error code ...
    int status = 0;
    int type   = 0;
    status     = __ffmahd(FPTR(m_fitsfile), m_hdunum+1, &type, &status);

    // If move was successful but HDU type in file differs from HDU type
    // of object then replace the HDU in the file
    if (status == 0 && type != exttype()) {
        status = __ffdhdu(FPTR(m_fitsfile), NULL, &status);
        if (status != 0) {
            throw GException::fits_error(G_SAVE_IMAGE, status);
        }
        status = __ffiimg(FPTR(m_fitsfile), m_bitpix, m_naxis, m_naxes, &status);
        //status = __ffiimgll(FPTR(m_fitsfile), m_bitpix, m_naxis, m_naxes, &status);
        if (status != 0) {
            throw GException::fits_error(G_SAVE_IMAGE, status);
        }
    }

    // If HDU does not yet exist in file then create it now
    if (status == 107) {
        status = 0;
        status = __ffcrim(FPTR(m_fitsfile), m_bitpix, m_naxis, m_naxes, &status);
        if (status != 0) {
            throw GException::fits_error(G_SAVE_IMAGE, status);
        }
    }
    else if (status != 0) {
        throw GException::fits_error(G_SAVE_IMAGE, status);
    }

    // If HDU seems to be empty then create it now. This is only needed for the
    // primary HDU, since __ffmahd gives no error if the primary HDU is empty.
    // By checking the number of keywords in the HDU we detect an empty HDU ...
    int num = 0;
    status  = __ffghsp(FPTR(m_fitsfile), &num, NULL, &status);
    if (status != 0) {
        throw GException::fits_error(G_SAVE_IMAGE, status);
    }
    if (num == 0) {
        status = __ffcrim(FPTR(m_fitsfile), m_bitpix, m_naxis, m_naxes, &status);
        if (status != 0) {
            throw GException::fits_error(G_SAVE_IMAGE, status);
        }
    }

    // Make sure that the image on disk has the right size by resizing it
    status = __ffrsim(FPTR(m_fitsfile), m_bitpix, m_naxis, m_naxes, &status);
    if (status != 0) {
        throw GException::fits_error(G_SAVE_IMAGE, status);
    }

    // Save the image pixels (if there are some ...)
    if (m_naxis > 0 && pixels != NULL) {
        long* fpixel = new long[m_naxis];
        long* lpixel = new long[m_naxis];
        for (int i = 0; i < m_naxis; ++i) {
            fpixel[i] = 1;
            lpixel[i] = m_naxes[i];
        }
        status = __ffpss(FPTR(m_fitsfile), datatype, fpixel, lpixel,
                         (void*)pixels, &status);
        delete [] fpixel;
        delete [] lpixel;
        if (status != 0) {
            throw GException::fits_error(G_SAVE_IMAGE, status);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fetch image pixels
 *
 * Fetch the image pixels. This function is in general called if pixel
 * values should be read or written yet no pixel array is allocated. In case
 * that pixels existed already before they will be deleted before fetching
 * new ones.
 * There are two possibilities to fetch the pixels:
 * (1) In case that a FITS file is attached to the image, the pixel array
 * will be loaded from the FITS file using the load_image() method.
 * (2) In case that no FITS file is attached, a new pixel array will be
 * allocated that is initalised to zero.
 ***************************************************************************/
void GFitsImage::fetch_data(void)
{
    // Fetch only if there are pixels in image
    if (m_num_pixels > 0) {

        // Allocate and initialise fresh memory
        alloc_data();
        init_data();

        // If a FITS file is attached then load pixels from FITS file.
        if (FPTR(m_fitsfile)->Fptr != NULL) {
            load_image(type(), ptr_data(), ptr_nulval(), &m_anynul);
        }

    } // endif: there were pixels available

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pixel offset
 *
 * @param[in] ix Pixel index [0,...,m_num_pixels-1].
 *
 * @exception GException::out_of_range
 *            Pixel index is outside valid range.
 ***************************************************************************/
int GFitsImage::offset(const int& ix) const
{
    // Check if axis is within the range
    #if defined(G_RANGE_CHECK)
    if (ix < 0 || ix >= m_num_pixels) {
        throw GException::out_of_range(G_OFFSET_1D, "Pixel index", ix, m_num_pixels);
    }
    #endif

    // Return index
    return ix;
}


/***********************************************************************//**
 * @brief Return 2D pixel offset
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 *
 * @exception GException::invalid_argument
 *            Pixel array has less than 2 dimensions.
 * @exception GException::out_of_range
 *            Image axis not valid.
 *
 * Computes of offset in the pixel array of a 2D coordinate. This method is
 * only applicable to pixels arrays that have at least 2 dimenions.
 ***************************************************************************/
int GFitsImage::offset(const int& ix, const int& iy) const
{
    // Operator is only valid for 2D images
    if (m_naxis < 2) {
        std::string msg = "2D pixel access operator used for image with "+
                          gammalib::str(m_naxis)+" dimensions. Please use "
                          "the correct pixel access operator.";
        throw GException::invalid_argument(G_OFFSET_2D, msg);
    }

    // Check if axis is within the range
    #if defined(G_RANGE_CHECK)
    if (ix < 0 || ix >= m_naxes[0]) {
        throw GException::out_of_range(G_OFFSET_2D,
                                       "Pixel index for first dimension",
                                       ix, m_naxes[0]);
    }
    if (iy < 0 || iy >= m_naxes[1]) {
        throw GException::out_of_range(G_OFFSET_2D,
                                       "Pixel index for second dimension",
                                       iy, m_naxes[1]);
    }
    #endif

    // Return offset
    return (ix + iy * m_naxes[0]);
}


/***********************************************************************//**
 * @brief Return 3D pixel offset
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 *
 * @exception GException::invalid_argument
 *            Pixel array has less than 3 dimensions.
 * @exception GException::out_of_range
 *            Image axis not valid.
 *
 * Computes of offset in the pixel array of a 3D coordinate. This method is
 * only applicable to pixels arrays that have at least 3 dimenions.
 ***************************************************************************/
int GFitsImage::offset(const int& ix, const int& iy, const int& iz) const
{
    // Operator is only valid for 3D images
    if (m_naxis < 3) {
        std::string msg = "3D pixel access operator used for image with "+
                          gammalib::str(m_naxis)+" dimensions. Please use "
                          "the correct pixel access operator.";
        throw GException::invalid_argument(G_OFFSET_3D, msg);
    }

    // Check if axis is within the range
    #if defined(G_RANGE_CHECK)
    if (ix < 0 || ix >= m_naxes[0]) {
        throw GException::out_of_range(G_OFFSET_3D,
                                       "Pixel index for first dimension",
                                       ix, m_naxes[0]);
    }
    if (iy < 0 || iy >= m_naxes[1]) {
        throw GException::out_of_range(G_OFFSET_3D,
                                       "Pixel index for second dimension",
                                       iy, m_naxes[1]);
    }
    if (iz < 0 || iz >= m_naxes[2]) {
        throw GException::out_of_range(G_OFFSET_3D,
                                       "Pixel index for third dimension",
                                       iz, m_naxes[2]);
    }
    #endif

    // Return offset
    return (ix + m_naxes[0] * (iy + iz * m_naxes[1]));
}


/***********************************************************************//**
 * @brief Return 4D pixel offset
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 * @param[in] it Pixel index in forth dimension (starting from 0).
 *
 * @exception GException::invalid_argument
 *            Pixel array has less than 4 dimensions.
 * @exception GException::out_of_range
 *            Image axis not valid.
 *
 * Computes of offset in the pixel array of a 4D coordinate. This method is
 * only applicable to pixels arrays that have at least 4 dimenions.
 ***************************************************************************/
int GFitsImage::offset(const int& ix, const int& iy, const int& iz,
                       const int& it) const
{
    // Operator is only valid for 4D images
    if (m_naxis < 4) {
        std::string msg = "4D pixel access operator used for image with "+
                          gammalib::str(m_naxis)+" dimensions. Please use "
                          "the correct pixel access operator.";
        throw GException::invalid_argument(G_OFFSET_4D, msg);
    }

    // Check if axis is within the range
    #if defined(G_RANGE_CHECK)
    if (ix < 0 || ix >= m_naxes[0]) {
        throw GException::out_of_range(G_OFFSET_4D,
                                       "Pixel index for first dimension",
                                       ix, m_naxes[0]);
    }
    if (iy < 0 || iy >= m_naxes[1]) {
        throw GException::out_of_range(G_OFFSET_4D,
                                       "Pixel index for second dimension",
                                       iy, m_naxes[1]);
    }
    if (iz < 0 || iz >= m_naxes[2]) {
        throw GException::out_of_range(G_OFFSET_4D,
                                       "Pixel index for third dimension",
                                       iz, m_naxes[2]);
    }
    if (it < 0 || it >= m_naxes[3]) {
        throw GException::out_of_range(G_OFFSET_4D,
                                       "Pixel index for third dimension",
                                       it, m_naxes[3]);
    }
    #endif

    // Return offset
    return (ix + m_naxes[0] * (iy + m_naxes[1] * (iz + it *  m_naxes[2])));
}
