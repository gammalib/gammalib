/***************************************************************************
 *             GSkymap.cpp  -  Class that implements a sky map             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
#include "GSkymap.hpp"
#include "GWcsHPX.hpp"
#include "GFits.hpp"
#include "GFitsTableDblCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT_HPX "GSkymap::GSkymap(std::string,std::string,int,std::string,int)"
#define G_OP_ACCESS                              "GSkymap::operator(int,int)"
#define G_READ                               "GSkymap::read(const GFitsHDU*)"
#define G_PIX2DIR                                     "GSkymap::pix2dir(int)"
#define G_DIR2PIX                                 "GSkymap::dir2pix(GSkyDir)"
#define G_OMEGA                                         "GSkymap::omega(int)"
#define G_READ_HEALPIX                     "GSkymap::read_healpix(GFitsHDU*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_READ_HEALPIX_DEBUG 0  // Debug read_healpix

/* __ Prototype __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                      GSkymap constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GSkymap::GSkymap(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS file constructor
 *
 * @param[in] filename FITS file from which sky map should be instantiated.
 ***************************************************************************/
GSkymap::GSkymap(const std::string& filename)
{
    // Initialise class members for clean destruction
    init_members();

    // Load skymap
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix constructor
 *
 * @param[in] wcs World Coordinate System (HPX).
 * @param[in] coords Coordinate System (CEL or GAL).
 * @param[in] nside Nside parameter.
 * @param[in] order Pixel ordering (RING or NEST).
 * @param[in] nmaps Number of maps in set (default=1).
 *
 * @exception GException::wcs_invalid 
 *            Invalid wcs parameter.
 * @exception GException::skymap_bad_nmaps 
 *            Invalid nmaps parameter.
 * @exception GException::wcs_hpx_bad_nside 
 *            Invalid nside parameter.
 * @exception GException::wcs_bad_coords 
 *            Invalid coordsys parameter.
 * @exception GException::wcs_hpx_bad_ordering 
 *            Invalid ordering parameter.
 ***************************************************************************/
GSkymap::GSkymap(const std::string& wcs, const std::string& coords,
                 const int& nside, const std::string& order,
                 const int nmaps)
{
    // Initialise class members for clean destruction
    init_members();

    // Check if wcs is HPX
    if (toupper(wcs) != "HPX")
        throw GException::wcs_invalid(G_CONSTRUCT_HPX, wcs,
                                      "WCS parameter must be 'HPX'.");

    // Check if nmaps parameter is >0
    if (nmaps < 1)
        throw GException::skymap_bad_nmaps(G_CONSTRUCT_HPX, nmaps,
                                           "nmaps parameter must be >0.");

    // Allocate WCS
    m_wcs = new GWcsHPX(nside, order, coords);

    // Set number of maps
    m_num_pixels = ((GWcsHPX*)m_wcs)->npix();
    m_num_maps   = nmaps;

    // Allocate pixels
    alloc_pixels();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Image constructor
 *
 * @param[in] wcs World Coordinate System (HPX).
 * @param[in] coords Coordinate System (CEL or GAL).
 * @param[in] dir Centre of skymap.
 * @param[in] nlon Number of pixels in longitude.
 * @param[in] nlat Number of pixels in latitude.
 * @param[in] dlon Size of pixels in longitude.
 * @param[in] dlat Size of pixels in latitude.
 * @param[in] nmaps Number of maps in set (default=1).
 ***************************************************************************/
GSkymap::GSkymap(const std::string& wcs, const std::string& coords,
                 GSkyDir& dir, const int& nlon, const int& nlat,
                 const double& dlon, const double& dlat, const int nmaps)
{
    // Initialise class members for clean destruction
    init_members();

    //TODO: Implement non HEALPix constructor

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] map Sky map from which class should be instantiated.
 ***************************************************************************/
GSkymap::GSkymap(const GSkymap& map)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(map);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkymap::~GSkymap(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GSkymap operators                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Pixel access operator
 *
 * @param[in] pixel pixel number (starting from 0).
 * @param[in] element vector element number (starting from 0).
 ***************************************************************************/
double& GSkymap::operator() (int pixel, int element)
{
    // Throw error if pixel is not in range
    #if defined(G_RANGE_CHECK)
    if (pixel   < 0 || pixel   >= m_num_pixels ||
        element < 0 || element >= m_num_maps)
        throw GException::out_of_range(G_OP_ACCESS, pixel, element,
                                       m_num_pixels, m_num_maps);
    #endif

    // Return reference to pixel value
    return m_pixels[pixel*m_num_maps+element];
}


/***********************************************************************//**
 * @brief Pixel access operator
 *
 * @param[in] pixel pixel number (starting from 0).
 * @param[in] element vector element number (starting from 0).
 ***************************************************************************/
const double& GSkymap::operator() (int pixel, int element) const
{
    // Throw error if pixel is not in range
    #if defined(G_RANGE_CHECK)
    if (pixel   < 0 || pixel   >= m_num_pixels ||
        element < 0 || element >= m_num_maps)
        throw GException::out_of_range(G_OP_ACCESS, pixel, element, 
                                       m_num_pixels, m_num_maps);
    #endif

    // Return reference to pixel value
    return m_pixels[pixel*m_num_maps+element];
}


/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] map GSkymap instance to be assigned
 ***************************************************************************/
GSkymap& GSkymap::operator= (const GSkymap& map)
{
    // Execute only if object is not identical
    if (this != &map) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(map);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GSkymap public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Load skymap from FITS file.
 *
 * @param[in] filename FITS file into which the skymap will be saved.
 *
 * Loads HEALPix and non HEALPix skymaps. First searches for HEALPix map in
 * FITS file by scanning all HDUs for PIXTYPE=HEALPIX. If no HEALPix map has
 * been found then search load first non-empty image.
 ***************************************************************************/
void GSkymap::load(const std::string& filename)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Open FITS file
    GFits fits(filename);

    // Get number of HDUs
    int num = fits.num_hdus();

    // Initialize pointer to HDU and load flag
    GFitsHDU* hdu    = NULL;
    int       loaded = 0;

    // First search for HEALPix extension. We can skip the first extension
    // since this is always an image and a HEALPix map is stored in a
    // binary table
    for (int extno = 1; extno < num; ++extno) {

        // Get pointer to HDU
        hdu = fits.hdu(extno);

        // Check if PIXTYPE keyword equals "HEALPIX"
        try {
            if (hdu->card("PIXTYPE")->string() == "HEALPIX") {
                read_healpix(hdu);
                loaded = 1;
                break;
            }
        }
        catch (GException::fits_key_not_found &e) {
        }

    } // endfor: looped over HDUs

    // If we have not found a HEALPIX map then search now for image.
    // Skip empty images
    if (loaded == 0) {
        for (int extno = 0; extno < num; ++extno) {
            //TODO: Implement non HEALPix loading
        }
    }

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save skymap into FITS file.
 *
 * @param[in] filename FITS file into which the skymap will be saved.
 * @param[in] clobber Overwrite existing file (0=false, 1=true).
 ***************************************************************************/
void GSkymap::save(const std::string& filename, int clobber)
{
    // Continue only if we have data to save
    if (m_wcs != NULL) {

        // Initialise HDU pointer
        GFitsHDU* hdu = NULL;

        // Case A: Skymap is Healpix
        if (m_wcs->type() == "HPX")
            hdu = create_healpix_hdu();

        // Case B: Skymap is not Healpix
        else {
            //TODO: Implement non HEALPix saving
        }

        // Create FITS file and save it to disk
        if (hdu != NULL) {
            GFits fits;
            fits.append_hdu(*hdu);
            fits.saveto(filename, clobber);
        }

    } // endif: we had data to save

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read skymap from FITS HDU
 *
 * @param[in] hdu FITS HDU from which skymap should be read.
 ***************************************************************************/
void GSkymap::read(const GFitsHDU* hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Continue only if HDU pointer is valid
    if (hdu != NULL) {

        // Initialize load flag
        int loaded = 0;

        // Try load as HEALPix map
        try {
            if (hdu->card("PIXTYPE")->string() == "HEALPIX") {
                read_healpix(hdu);
                loaded = 1;
            }
        }
        catch (GException::fits_key_not_found &e) {
        }

        // ... otherwise try loading as non HEALPix map
        if (loaded == 0) {
            //TODO: Implement non HEALPix loading
        }

    } // endif: HDU pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write skymap into FITS HDU
 *
 * @param[in] hdu FITS HDU into which skymap should be written.
 ***************************************************************************/
void GSkymap::write(GFitsHDU* hdu)
{
    // TODO: Skymap writing

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pix Pixel number (0,1,...,num_pixels).
 *
 * @exception GException::wcs
 *            No valid WCS found.
 ***************************************************************************/
GSkyDir GSkymap::pix2dir(const int& pix)
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL)
        throw GException::wcs(G_PIX2DIR, "No valid WCS found.");

    // Return sky direction
    return (m_wcs->pix2dir(pix));
}


/***********************************************************************//**
 * @brief Returns pixel index for a given sky direction
 *
 * @param[in] dir Sky direction.
 *
 * @exception GException::wcs
 *            No valid WCS found.
 ***************************************************************************/
int GSkymap::dir2pix(GSkyDir dir) const
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL)
        throw GException::wcs(G_DIR2PIX, "No valid WCS found.");

    // Return pixel index
    return (m_wcs->dir2pix(dir));
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pix Pixel number (0,1,...,num_pixels).
 *
 * @exception GException::wcs
 *            No valid WCS found.
 ***************************************************************************/
double GSkymap::omega(const int& pix) const
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL)
        throw GException::wcs(G_OMEGA, "No valid WCS found.");

    // Return solid angle
    return (m_wcs->omega(pix));
}


/***********************************************************************//**
 * @brief Returns number of pixels
 ***************************************************************************/
int GSkymap::npix(void) const
{
    // Return number of pixels
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Returns number of maps
 ***************************************************************************/
int GSkymap::nmaps(void) const
{
    // Return number of maps
    return m_num_maps;
}


/*==========================================================================
 =                                                                         =
 =                          GSkymap private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkymap::init_members(void)
{
    // Initialise members
    m_num_pixels = 0;
    m_num_maps   = 0;
    m_wcs        = NULL;
    m_pixels     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate skymap pixels
 ***************************************************************************/
void GSkymap::alloc_pixels(void)
{
    // Compute data size
    int size = m_num_pixels * m_num_maps;

    // Continue only if there are pixels
    if (size > 0) {

        // Allocate pixels and initialize them to 0
        m_pixels = new double[size];
        for (int i = 0; i < size; ++i)
            m_pixels[i] = 0.0;

    } // endif: there were pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Sky map from which members should be copied.
 ***************************************************************************/
void GSkymap::copy_members(const GSkymap& map)
{
    // Copy attributes
    m_num_pixels = map.m_num_pixels;
    m_num_maps   = map.m_num_maps;

    // Clone WCS if it is valid
    if (map.m_wcs != NULL) m_wcs = map.m_wcs->clone();

    // Compute data size
    int size = m_num_pixels * m_num_maps;

    // Copy pixels
    if (size > 0 && map.m_pixels != NULL) {
        alloc_pixels();
        for (int i = 0; i <  size; ++i)
            m_pixels[i] = map.m_pixels[i];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkymap::free_members(void)
{
    // Free memory
    if (m_wcs    != NULL) delete m_wcs;
    if (m_pixels != NULL) delete [] m_pixels;

    // Signal free pointers
    m_wcs        = NULL;
    m_pixels     = NULL;
    m_num_pixels = 0;
    m_num_maps   = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Healpix data from FITS table.
 *
 * @param[in] hdu FITS HDU containing the Healpix data.
 *
 * @exception GException::skymap
 *            HEALPix HDU is not a binary table.
 *
 * HEALPix data may be stored in various formats depending on the 
 * application that has writted the data. HEALPix IDL, for example, may
 * store the data in vectors of length 1024 if the number of pixels is
 * a multiple of 1024. On the other hand, vectors may also be used to store
 * several HEALPix maps into a single column. Alternatively, multiple maps
 * may be stored in multiple columns.
 ***************************************************************************/
void GSkymap::read_healpix(const GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Verify that HDU contains binary table
        if (hdu->exttype() != 2) //TODO: Implement unambigous method in GFitsHDU
            throw GException::skymap(G_READ_HEALPIX,
                                     "HEALPix HDU is not a binary table.");

        // Determine number of rows and columns in table
        int nrows = ((GFitsTable*)hdu->data())->nrows();
        int ncols = ((GFitsTable*)hdu->data())->ncols();
        #if G_READ_HEALPIX_DEBUG
        std::cout << "nrows=" << nrows << " ncols=" << ncols << std::endl;
        #endif

        // Allocate Healpix WCS
        m_wcs = new GWcsHPX;

        // Read WCS information from FITS header
        m_wcs->read(hdu);

        // Set number of pixels based on NSIDE parameter
        m_num_pixels = ((GWcsHPX*)m_wcs)->npix();
        #if G_READ_HEALPIX_DEBUG
        std::cout << "m_num_pixels=" << m_num_pixels << std::endl;
        #endif

        // Number of map pixels has to be a multiple of the number of
        // rows in column
        if (m_num_pixels % nrows != 0)
            throw GException::skymap_bad_size(G_READ_HEALPIX, nrows,
                                              m_num_pixels);

        // Determine vector length for HEALPix data storage
        int nentry = m_num_pixels / nrows;
        #if G_READ_HEALPIX_DEBUG
        std::cout << "nentry=" << nentry << std::endl;
        #endif

        // Determine number of maps from NBRBINS keyword. If keyword was
        // not found then determine the number of maps that fit into
        // all columns. Only count columns that can full hold the map.
        try {
            m_num_maps = hdu->card("NBRBINS")->integer();
        }
        catch (GException::fits_key_not_found &e) {
            m_num_maps = 0;
            for (int icol = 0; icol < ncols; ++icol) {
                GFitsTableCol* col = hdu->column(icol);
                if (col->number() % nentry == 0)
                    m_num_maps += col->number() / nentry;
            }
        }
        #if G_READ_HEALPIX_DEBUG
        std::cout << "m_num_maps=" << m_num_maps << std::endl;
        #endif

        // Allocate pixels to hold the map
        alloc_pixels();

        // Initialise map counter
        int imap = 0;

        // Loop over all columns
        for (int icol = 0; icol < ncols; ++icol) {

            // Get next column
            GFitsTableCol* col = hdu->column(icol);

            // Only consider columns that can fully hold maps
            if (col->number() % nentry == 0) {

                // Determine number of maps in column
                int num = col->number() / nentry;

                // Loop over all maps in column
                int inx_start = 0;
                int inx_end   = nentry;
                for (int i = 0; i < num; ++i) {

                    // Load map
                    double *ptr = m_pixels + imap;
                    for (int row = 0; row < col->length(); ++row) {
                        for (int inx = inx_start; inx < inx_end; ++inx, ptr+=m_num_maps)
                                *ptr = col->real(row,inx);
                    }
                    #if G_READ_HEALPIX_DEBUG
                    std::cout << "Load map=" << imap << " index="
                              << inx_start << "-" << inx_end << std::endl;
                    #endif

                    // Increment index range
                    inx_start  = inx_end;
                    inx_end   += nentry;

                    // Increment map counter
                    imap++;

                    // Break if we have loaded all maps
                    if (imap >= m_num_maps)
                        break;

                } // endfor: looped over all maps in column
            } // endif: column could fully hold maps

            // Break if we have loaded all maps
            if (imap >= m_num_maps)
                break;

        } // endfor: looped over all columns

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create FITS HDU containing Healpix data
 *
 * Returns pointer to HDU that contains the Healpix data. Deallocation of the
 * GFitsHDU object has to be done by the client.
 ***************************************************************************/
GFitsHDU* GSkymap::create_healpix_hdu(void)
{
    // Initialise result to NULL pointer
    GFitsHDU* hdu = NULL;

    // Compute size of Healpix data
    int size = m_num_pixels * m_num_maps;

    // Continue only if we have pixels
    if (size > 0) {

        // Set number of rows and columns
        int rows   = m_num_pixels;
        int number = m_num_maps;

        // Create column to hold Healpix data
        GFitsTableDblCol column = GFitsTableDblCol("DATA", rows, number);

        // Fill data into column
        double* ptr = m_pixels;
        for (int row = 0; row < rows; ++row) {
            for (int inx = 0; inx < number; ++inx, ++ptr)
                column(row,inx) = *ptr;
        }

        // Create HDU that contains Healpix map in a binary table
        GFitsBinTable table = GFitsBinTable(rows);
        table.append_column(column);
        hdu = new GFitsHDU(table);

    } // endif: there were pixels

    // ... otherwise create an empty header
    else
        hdu = new GFitsHDU;

    // Set extension name
    hdu->extname("HEALPIX");

    // If we have WCS information then write into FITS header
    if (m_wcs != NULL) m_wcs->write(hdu);

    // Set additional keywords
    hdu->card("NBRBINS", m_num_maps, "Number of HEALPix maps");

    // Return HDU
    return hdu;
}


/*==========================================================================
 =                                                                         =
 =                              GSkymap friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] map Sky map to put in output stream.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSkymap& map)
{
    // Put header in stream
    os << "=== GSkymap ===" << std::endl;

    // Put WCS information in stream
    if (map.m_wcs != NULL) {
        if (map.m_wcs->type() == "HPX")
            os << *((GWcsHPX*)map.m_wcs);
        //os << std::endl;
    }

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                      Other functions used by GSkymap                    =
 =                                                                         =
 ==========================================================================*/
