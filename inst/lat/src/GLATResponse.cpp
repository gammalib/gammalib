/***************************************************************************
 *               GLATResponse.cpp  -  GLAST LAT Response class             *
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
 * @file GLATResponse.cpp
 * @brief GLATResponse class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GVector.hpp"
#include "GLATResponse.hpp"
#include "GLATPointing.hpp"
#include "GFits.hpp"
#include "GFitsHDU.hpp"
#include "GFitsImageDbl.hpp"
#include "GFitsTableFltCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD            "load(const std::string&, const std::string&)"
#define G_GET_FITS_VECTOR "GLATResponse::get_fits_vector(GFitsHDU*,std::string&,int)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GLATResponse::GLATResponse(void) : GResponse()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
GLATResponse::GLATResponse(const GLATResponse& rsp) : GResponse(rsp)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATResponse::~GLATResponse(void)
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
 * @brief Assignment operator
 *
 * @param[in] rsp Response to be assigned
 ***************************************************************************/
GLATResponse& GLATResponse::operator= (const GLATResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GResponse::operator=(rsp);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(rsp);

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
 * @brief Return value of instrument response function.
 *
 * @param[in] obsDir Observed photon direction
 * @param[in] obsEng Observed energy of photon
 * @param[in] srcDir True photon direction
 * @param[in] srcEng True energy of photon
 * @param[in] pnt Pointer to instrument pointing information
 * @param[in] time Photon arrival time
 *
 * @todo Implement method (just a dummy for the moment)
 ***************************************************************************/
double GLATResponse::irf(GSkyDir& obsDir, const GEnergy& obsEng,
                         GSkyDir& srcDir, const GEnergy& srcEng,
                         const GPointing* pnt, const GTime& time)
{
    // Get IRF components
    double irf;
    irf  = psf(obsDir, obsEng, srcDir, srcEng, pnt, time);
    irf *= aeff(obsDir, obsEng, srcDir, srcEng, pnt, time);
    irf *= edisp(obsDir, obsEng, srcDir, srcEng, pnt, time);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Set the path to the calibration database.
 *
 * @param[in] caldb Absolute path to calibration database
 *
 * NOTE: So far no check is done on whether the path exists!
 ***************************************************************************/
void GLATResponse::set_caldb(const std::string& caldb)
{
    // Simply copy path
    /// @todo Some check should be done on whether the path exists
    m_caldb = caldb;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load a specified LAT response function.
 *
 * @param[in] rspname Name of response
 * @param[in] rsptype Type of response ('front' or 'back')
 *
 * Loads the specified GLAST LAT response from the calibration database and
 * performs some pre-calculations for faster response determination.
 ***************************************************************************/
void GLATResponse::load(const std::string& rspname, const std::string& rsptype)
{
    // Store response name
    m_rspname = rspname;

    // Convert response type string to section
    if (rsptype == "front")
        m_section = 0;
    else if (rsptype == "back")
        m_section = 1;
    else
        throw GException::rsp_invalid_type(G_LOAD, rsptype);

    // Store response type
    m_rsptype = rsptype;

    // Initialise effective area
    aeff_init();

    // Initialise PSF
    psf_init();

    // Initialise energy dispersion
    edisp_init();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save response
 *
 * @param[in] rspname Filename into which the response will be saved
 ***************************************************************************/
void GLATResponse::save(const std::string& rspname) const
{
    // Allocate FITS file
    GFits file;
    
    // Open FITS file
    file.open(rspname);
    
    // Save effective area
    aeff_append(file);
    
    // Append PSF HDUs
    psf_append(file);
    
    // Save energy dispersions
    edisp_append(file);

    // Save FITS file
    file.save();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATResponse::init_members(void)
{
    // Initialise effective area members
    aeff_init_members();

    // Initialise PSF members
    psf_init_members();
    
    // Initialize energy dispersion members
    edisp_init_members();

    // By default use HANDOFF response database.
    char* handoff = getenv("HANDOFF_IRF_DIR");
    if (handoff != NULL)
        m_caldb.assign(handoff);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
void GLATResponse::copy_members(const GLATResponse& rsp)
{
    // Copy effective area members
    aeff_copy_members(rsp);
    
    // Copy PSF members
    psf_copy_members(rsp);

    // Copy energy dispersion members
    edisp_copy_members(rsp);
    
    // Copy remaining members
    m_caldb = rsp.m_caldb;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATResponse::free_members(void)
{
    // Free effective area members
    aeff_free_members();

    // Free PSF  members
    psf_free_members();
    
    // Free energy dispersion members
    edisp_free_members();
    
    // Free remaining members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone response class
***************************************************************************/
GLATResponse* GLATResponse::clone(void) const
{
    return new GLATResponse(*this);
}


/***********************************************************************//**
 * @brief Load floating point data from vector column into vector
 *
 * @param[in] hdu Pointer to HDU from which data should be loaded
 * @param[in] colname Column name from which data should be loaded
 * @param[in] row Table row from which data should be loaded
 ***************************************************************************/
GVector GLATResponse::get_fits_vector(const GFitsHDU* hdu, 
                                      const std::string& colname, 
                                      int row)
{
    // Get pointer to
    GFitsTableCol* ptr = hdu->column(colname);
    if (ptr == NULL)
        throw GException::fits_column_not_found(G_GET_FITS_VECTOR, colname);

    // Determine number of entries
    int num = ptr->number();

    // Load data into vector
    GVector data(num);
    for (int i = 0; i < num; ++i)
        data(i) = ptr->real(row,i);

    // Return vector
    return data;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
