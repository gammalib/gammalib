/***************************************************************************
 *               GLATResponse.cpp  -  GLAST LAT Response class             *
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
#include "GVector.hpp"
#include "GLATResponse.hpp"
#include "GFits.hpp"
#include "GFitsTableFltCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_INIT_AEFF       "GLATResponse::init_aeff()"
#define G_INIT_PSF        "GLATResponse::init_psf()"
#define G_INIT_EDISP      "GLATResponse::init_edisp()"
#define G_GET_FITS_VECTOR "GLATResponse::get_fits_vector(const GFitsHDU*, const std::string&, int)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                   GLATResponse constructors/destructors                 =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GLATResponse::GLATResponse() : GResponse()
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
GLATResponse::GLATResponse(const GLATResponse& rsp) : GResponse(rsp)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***************************************************************************
 *                               Destructor                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GLATResponse::~GLATResponse()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GLATResponse operators                         =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
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
 =                        GLATResponse public methods                      =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *               Return value of instrument response function              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GLATResponse::irf(const GSkyDir& obsDir, const double& obsEng,
                         const GSkyDir& srcDir, const double& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const double& time)
{
    // Return IRF value
    return 0.0;
}


/***************************************************************************
 *                           Return effective area                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GLATResponse::aeff(const GSkyDir& obsDir, const double& obsEng,
                          const GSkyDir& srcDir, const double& srcEng,
                          const GSkyDir& instPntDir, const double& instPosAng,
                          const double& time)
{
    // Return Aeff value
    return 0.0;
}


/***************************************************************************
 *                   Return value of point spread function                 *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GLATResponse::psf(const GSkyDir& obsDir, const double& obsEng,
                         const GSkyDir& srcDir, const double& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const double& time)
{
    // Return PSF value
    return 0.0;
}


/***************************************************************************
 *                     Return value of energy dispersion                   *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GLATResponse::edisp(const GSkyDir& obsDir, const double& obsEng,
                           const GSkyDir& srcDir, const double& srcEng,
                           const GSkyDir& instPntDir, const double& instPosAng,
                           const double& time)
{
    // Return Edisp value
    return 0.0;
}


/***************************************************************************
 *                         Set calibration database                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @brief Set the path to the calibration database.
 */
void GLATResponse::set_caldb(std::string caldb)
{
    // Simply copy path
    // NOTE: Some check should be done on whether the path exists !!!
    m_caldb = caldb;

    // Return
    return;
}


/***************************************************************************
 *                              Load response                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @brief Load a specified LAT response function.
 */
void GLATResponse::load(std::string rspname)
{
    // Store response name
    m_rspname = rspname;

    // Initialise effective area
    init_aeff();

    // Initialise PSF
    init_psf();

    // Initialise energy dispersion
    init_edisp();
}


/***************************************************************************
 *                              Save response                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @brief Save me my friend ...
 */
void GLATResponse::save(std::string rspname) const
{
}


/*==========================================================================
 =                                                                         =
 =                        GLATResponse private methods                     =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                        Initialise effective area                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::init_aeff(void)
{
    // Build filename
    std::string filename = "aeff_"  + m_rspname + ".fits";

    // Open FITS file
    GFits file;
    file.open(m_caldb + "/" + filename);

    // Get pointer to effective area HDU
    GFitsHDU* hdu = file.hdu("EFFECTIVE AREA");
    if (hdu == NULL)
        throw GException::fits_hdu_not_found(G_INIT_AEFF, "EFFECTIVE AREA");

    // Get the data
    GVector energ_lo  = get_fits_vector(hdu, "ENERG_LO");
    GVector energ_hi  = get_fits_vector(hdu, "ENERG_HI");
    GVector ctheta_lo = get_fits_vector(hdu, "CTHETA_LO");
    GVector ctheta_hi = get_fits_vector(hdu, "CTHETA_HI");
    GVector effarea   = get_fits_vector(hdu, "EFFAREA");

    // Return
    return;
}


/***************************************************************************
 *                             Initialise PSF                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::init_psf(void)
{
    // Build filename
    std::string filename = "psf_"  + m_rspname + ".fits";

    // Open FITS file
    GFits file;
    file.open(m_caldb + "/" + filename);

    // Get pointer to PSF HDU
    GFitsHDU* hdu = file.hdu("RPSF");
    if (hdu == NULL)
        throw GException::fits_hdu_not_found(G_INIT_PSF, "RPSF");

    // Get the data
    GVector energ_lo  = get_fits_vector(hdu, "ENERG_LO");
    GVector energ_hi  = get_fits_vector(hdu, "ENERG_HI");
    GVector ctheta_lo = get_fits_vector(hdu, "CTHETA_LO");
    GVector ctheta_hi = get_fits_vector(hdu, "CTHETA_HI");
    GVector ncore     = get_fits_vector(hdu, "NCORE");
    GVector sigma     = get_fits_vector(hdu, "SIGMA");
    GVector gcore     = get_fits_vector(hdu, "GCORE");
    GVector gtail     = get_fits_vector(hdu, "GTAIL");

    // Return
    return;
}


/***************************************************************************
 *                      Initialise energy dispersion                       *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::init_edisp(void)
{
    // Build filename
    std::string filename = "edisp_"  + m_rspname + ".fits";

    // Open FITS file
    GFits file;
    file.open(m_caldb + "/" + filename);

    // Get pointer towards energy dispersion HDU
    GFitsHDU* hdu = file.hdu("ENERGY DISPERSION");
    if (hdu == NULL)
        throw GException::fits_hdu_not_found(G_INIT_EDISP, "ENERGY DISPERSION");

    // Get the data
    GVector energ_lo  = get_fits_vector(hdu, "ENERG_LO");
    GVector energ_hi  = get_fits_vector(hdu, "ENERG_HI");
    GVector ctheta_lo = get_fits_vector(hdu, "CTHETA_LO");
    GVector ctheta_hi = get_fits_vector(hdu, "CTHETA_HI");
    GVector dnorm     = get_fits_vector(hdu, "DNORM");
    GVector ltail     = get_fits_vector(hdu, "LTAIL");
    GVector rwidth    = get_fits_vector(hdu, "RWIDTH");
    GVector nr2       = get_fits_vector(hdu, "NR2");
    GVector lt2       = get_fits_vector(hdu, "LT2");
    GVector rt2       = get_fits_vector(hdu, "RT2");

    // ...

    // Return
    return;
}


/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::init_members(void)
{
    // Initialise members
    //m_aeff_name.clear();
    //m_psf_name.clear();
    //m_edisp_name.clear();

    // By default use HANDOFF response database.
    char* handoff = getenv("HANDOFF_IRF_DIR");
    if (handoff != NULL)
        m_caldb.assign(handoff);


    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::copy_members(const GLATResponse& rsp)
{
    // Copy attributes
    //m_aeff_name  = rsp.m_aeff_name;
    //m_psf_name   = rsp.m_psf_name;
    //m_edisp_name = rsp.m_edisp_name;

    // Copy other membres

    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GLATResponse::free_members(void)
{
    // Free memory

    // Signal free pointers

    // Return
    return;
}


/***************************************************************************
 *          Load floating point data from vector column into vector        *
 * ----------------------------------------------------------------------- *
 * All information is stored in the first table row using vector columns.  *
 ***************************************************************************/
GVector GLATResponse::get_fits_vector(const GFitsHDU* hdu, const std::string& colname, int row)
{
    // Get pointer to
    GFitsTableCol* ptr = hdu->column(colname);
    if (ptr == NULL)
        throw GException::fits_column_not_found(G_GET_FITS_VECTOR, colname);
        
    // Determine number of entries
    int num = ptr->repeat();

    // Load data into vector
    GVector data(num);
    for (int i = 0; i < num; ++i)
        data(i) = ptr->real(row,i);

    // Return vector
    return data;
}


/*==========================================================================
 =                                                                         =
 =                           GLATResponse friends                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GLATResponse                  =
 =                                                                         =
 ==========================================================================*/
