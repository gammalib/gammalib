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
//#include <vector>
#include <iostream>                           // cout, cerr
#include "GException.hpp"
#include "GVector.hpp"
#include "GLATResponse.hpp"
#include "GFits.hpp"
#include "GFitsTableFltCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

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

    // Get pointer towards effective area HDU
    GFitsHDU* hdu = file.hdu("EFFECTIVE AREA");
    if (hdu == NULL) {
        cout << endl << "TEST ERROR: Unable to find HDU <EFFECTIVE AREA>." << endl;
        throw;
    }

    // Get pointers to data columns
    GFitsTableFltCol* energ_lo  = (GFitsTableFltCol*)hdu->column("ENERG_LO");
    GFitsTableFltCol* energ_hi  = (GFitsTableFltCol*)hdu->column("ENERG_HI");
    GFitsTableFltCol* ctheta_lo = (GFitsTableFltCol*)hdu->column("CTHETA_LO");
    GFitsTableFltCol* ctheta_hi = (GFitsTableFltCol*)hdu->column("CTHETA_HI");
    GFitsTableFltCol* effarea   = (GFitsTableFltCol*)hdu->column("EFFAREA");
    if (energ_lo  == NULL || energ_hi  == NULL ||
        ctheta_lo == NULL || ctheta_hi == NULL || effarea == NULL) {
        cout << endl << "TEST ERROR: Unable to get <EFFECTIVE AREA> columns." << endl;
        throw;
    }

    // ...

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

    // Get pointer towards PSF HDU
    GFitsHDU* hdu = file.hdu("RPSF");
    if (hdu == NULL) {
        cout << endl << "TEST ERROR: Unable to find HDU <RPSF>." << endl;
        throw;
    }

    // Get data
    GVector energ_lo  = get_fits_column(hdu, "ENERG_LO");
/*
    GVector energ_hi  = get_fits_column(hdu, "ENERG_HI");
    GVector ctheta_lo = get_fits_column(hdu, "CTHETA_LO");
    GVector ctheta_hi = get_fits_column(hdu, "CTHETA_HI");
    GVector ncore     = get_fits_column(hdu, "NCORE");
    GVector sigma     = get_fits_column(hdu, "SIGMA");
    GVector gcore     = get_fits_column(hdu, "GCORE");
    GVector gtail     = get_fits_column(hdu, "GTAIL");
*/
    // ...

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
    if (hdu == NULL) {
        cout << endl << "TEST ERROR: Unable to find HDU <ENERGY DISPERSION>." << endl;
        throw;
    }

    // Get pointers to data columns
    GFitsTableFltCol* energ_lo  = (GFitsTableFltCol*)hdu->column("ENERG_LO");
    GFitsTableFltCol* energ_hi  = (GFitsTableFltCol*)hdu->column("ENERG_HI");
    GFitsTableFltCol* ctheta_lo = (GFitsTableFltCol*)hdu->column("CTHETA_LO");
    GFitsTableFltCol* ctheta_hi = (GFitsTableFltCol*)hdu->column("CTHETA_HI");
    GFitsTableFltCol* dnorm     = (GFitsTableFltCol*)hdu->column("DNORM");
    GFitsTableFltCol* ltail     = (GFitsTableFltCol*)hdu->column("LTAIL");
    GFitsTableFltCol* rwidth    = (GFitsTableFltCol*)hdu->column("RWIDTH");
    GFitsTableFltCol* nr2       = (GFitsTableFltCol*)hdu->column("NR2");
    GFitsTableFltCol* lt2       = (GFitsTableFltCol*)hdu->column("LT2");
    GFitsTableFltCol* rt2       = (GFitsTableFltCol*)hdu->column("RT2");
    if (energ_lo  == NULL || energ_hi  == NULL ||
        ctheta_lo == NULL || ctheta_hi == NULL ||
        dnorm     == NULL || ltail     == NULL || rwidth == NULL ||
        nr2       == NULL || lt2       == NULL || rt2    == NULL) {
        cout << endl << "TEST ERROR: Unable to get <ENERGY DISPERSION> columns." << endl;
        throw;
    }

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
 *              Load floating point data from column into vector           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GVector GLATResponse::get_fits_column(const GFitsHDU* hdu, const std::string& colname)
{
    // Get pointer on column
    GFitsTableFltCol* ptr  = (GFitsTableFltCol*)hdu->column(colname);
    if (ptr == NULL)
        throw;

    // Load data into vector
    int num = ptr->repeat();
    GVector data(num);
    for (int i = 0; i < num; ++i)
        data(i) = ptr->real(0,i);
    cout << data << endl;

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
