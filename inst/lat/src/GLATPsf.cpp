/***************************************************************************
 *             GLATPsf.cpp  -  Fermi LAT point spread function             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATPsf.cpp
 * @brief Fermi LAT point spread function class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATPsf.hpp"
#include "GLATResponse.hpp"
#include "GLATPointing.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                             "GLATPsf::read(const GFits* file)"
#define G_READ_PSF                           "GLATPsf::read_psf(GFitsTable*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
const double scale_factor_log = 39.810717;
const double scale_base_log   = 0.15848932;
const double scale_front_c0   = 3.77e-4;
const double scale_front_c1   = 5.8e-2;
const double scale_back_c0    = 1.3e-3;
const double scale_back_c1    = 9.6e-2;
const int    angle_num        = 5000;    // Number of angles for PSF
const double angle_min        = 0.0001;  // Minimum angular separation (rad)
const double angle_bin        = 0.01;    // Angular separation binning (rad)
const double max_sep          = pihalf;  // Maximum angular separation (rad)


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATPsf::GLATPsf(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename FITS file name.
 *
 * Construct instance by loading the point spread function from FITS file.
 ***************************************************************************/
GLATPsf::GLATPsf(const std::string filename)
{
    // Initialise class members
    init_members();

    // Load energy dispersion from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsf::GLATPsf(const GLATPsf& psf)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(psf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATPsf::~GLATPsf(void)
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
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsf& GLATPsf::operator= (const GLATPsf& psf)
{
    // Execute only if object is not identical
    if (this != &psf) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(psf);

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
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATPsf::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATPsf* GLATPsf::clone(void) const
{
    return new GLATPsf(*this);
}


/***********************************************************************//**
 * @brief Load point spread function from FITS file
 *
 * @param[in] filename FITS file.
 ***************************************************************************/
void GLATPsf::load(const std::string filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read point spread function from file
    read(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save point spread function into FITS file
 *
 * @param[in] filename FITS file.
 * @param[in] clobber Overwrite existing file?.
 ***************************************************************************/
void GLATPsf::save(const std::string filename, bool clobber)
{
    // Open FITS file
    GFits fits(filename);

    // Write point spread function into file
    write(fits);

    // Close FITS file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read point spread function from FITS file
 *
 * @param[in] fits FITS file.
 *
 * @exception GException::fits_hdu_not_found
 *            Effective area HDU not found in FITS file
 *
 * @todo Implement reading of scaling parameters
 ***************************************************************************/
void GLATPsf::read(const GFits& fits)
{
    // Clear instance
    clear();

    // Get pointer to effective area HDU
    GFitsTable* hdu_rpsf  = fits.table("RPSF");
    GFitsTable* hdu_scale = fits.table("PSF_SCALING_PARAMS");
    if (hdu_rpsf == NULL)
        throw GException::fits_hdu_not_found(G_READ, "RPSF");
    if (hdu_scale == NULL)
        throw GException::fits_hdu_not_found(G_READ, "PSF_SCALING_PARAMS");

    // Read point spread function
    read_psf(hdu_rpsf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write point spread function into FITS file
 *
 * @param[in] fits FITS file.
 ***************************************************************************/
void GLATPsf::write(GFits& fits) const
{
    // Write energy dispersion
    write_psf(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point spread function information
 ***************************************************************************/
std::string GLATPsf::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATPsf ===");
    result.append("\n"+parformat("Number of energy bins")+str(nenergies()));
    result.append("\n"+parformat("Number of cos theta bins")+str(ncostheta()));

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATPsf::init_members(void)
{
    // Initialise members
    m_rpsf_bins.clear();
    m_ncore.clear();
    m_sigma.clear();
    m_gcore.clear();
    m_gtail.clear();
    m_scale.clear();
    m_ltcube_logE = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GLATPsf::copy_members(const GLATPsf& psf)
{
    // Copy members
    m_rpsf_bins   = psf.m_rpsf_bins;
    m_ncore       = psf.m_ncore;
    m_sigma       = psf.m_sigma;
    m_gcore       = psf.m_gcore;
    m_gtail       = psf.m_gtail;
    m_scale       = psf.m_scale;
    m_ltcube_logE = psf.m_ltcube_logE;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATPsf::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read point spread function from FITS table
 *
 * @param[in] hdu FITS table pointer.
 *
 * @exception GException::fits_column_not_found
 *            Effective area column not found
 * @exception GLATException::inconsistent_response
 *            Inconsistent response table encountered
 ***************************************************************************/
void GLATPsf::read_psf(const GFitsTable* hdu)
{
    // Clear arrays
    m_ncore.clear();
    m_sigma.clear();
    m_gcore.clear();
    m_gtail.clear();

    // Get energy and cos theta binning
    m_rpsf_bins.read(hdu);

    // Continue only if there are bins
    int size = m_rpsf_bins.size();
    if (size > 0) {

        // Allocate arrays
        m_ncore.reserve(size);
        m_sigma.reserve(size);
        m_gcore.reserve(size);
        m_gtail.reserve(size);

        // Get pointer to columns
        GFitsTableCol* ncore = ((GFitsTable*)hdu)->column("NCORE");
        GFitsTableCol* sigma = ((GFitsTable*)hdu)->column("SIGMA");
        GFitsTableCol* gcore = ((GFitsTable*)hdu)->column("GCORE");
        GFitsTableCol* gtail = ((GFitsTable*)hdu)->column("GTAIL");
        if (ncore == NULL)
            throw GException::fits_column_not_found(G_READ_PSF, "NCORE");
        if (sigma == NULL)
            throw GException::fits_column_not_found(G_READ_PSF, "SIGMA");
        if (gcore == NULL)
            throw GException::fits_column_not_found(G_READ_PSF, "GCORE");
        if (gtail == NULL)
            throw GException::fits_column_not_found(G_READ_PSF, "GTAIL");

        // Check consistency of columns
        if (ncore->number() != size)
            throw GLATException::inconsistent_response(G_READ_PSF,
                                                       ncore->number(), size);
        if (sigma->number() != size)
            throw GLATException::inconsistent_response(G_READ_PSF,
                                                       sigma->number(), size);
        if (gcore->number() != size)
            throw GLATException::inconsistent_response(G_READ_PSF,
                                                       gcore->number(), size);
        if (gtail->number() != size)
            throw GLATException::inconsistent_response(G_READ_PSF,
                                                       gtail->number(), size);

        // Copy data
        for (int i = 0; i < size; ++i) {
            m_ncore.push_back(ncore->real(0,i));
            m_sigma.push_back(sigma->real(0,i));
            m_gcore.push_back(gcore->real(0,i));
            m_gtail.push_back(gtail->real(0,i));
        }

    } // endif: there were bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write point spread function into FITS file
 *
 * @param[in] file FITS file.
 *
 * This method does not write anything if the instance is empty.
 ***************************************************************************/
void GLATPsf::write_psf(GFits& file) const
{
    // Continue only if there are bins
    int size = m_rpsf_bins.size();
    if (size > 0) {

        // Create new binary table
        GFitsBinTable* hdu_rpsf = new GFitsBinTable;

        // Set table attributes
        hdu_rpsf->extname("RPSF");

        // Write boundaries into table
        m_rpsf_bins.write(hdu_rpsf);

        // Allocate floating point vector columns
        GFitsTableFloatCol col_ncore = GFitsTableFloatCol("NCORE",  1, size);
        GFitsTableFloatCol col_sigma = GFitsTableFloatCol("SIGMA",  1, size);
        GFitsTableFloatCol col_gcore = GFitsTableFloatCol("GCORE",  1, size);
        GFitsTableFloatCol col_gtail = GFitsTableFloatCol("GTAIL",  1, size);

        // Fill columns
        for (int i = 0; i < size; ++i) {
            col_ncore(0,i) = m_ncore[i];
            col_sigma(0,i) = m_sigma[i];
            col_gcore(0,i) = m_gcore[i];
            col_gtail(0,i) = m_gtail[i];
        }

        // Append columns to table
        hdu_rpsf->append_column(col_ncore);
        hdu_rpsf->append_column(col_sigma);
        hdu_rpsf->append_column(col_gcore);
        hdu_rpsf->append_column(col_gtail);
        
        // Append HDU to FITS file
        file.append(hdu_rpsf);

    } // endif: there were data to write

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] psf Point spread function.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATPsf& psf)
{
     // Write point spread function in output stream
    os << psf.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GLATPsf& psf)
{
    // Write point spread function into logger
    log << psf.print();

    // Return logger
    return log;
}
