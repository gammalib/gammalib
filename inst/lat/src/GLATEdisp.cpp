/***************************************************************************
 *         GLATEdisp.cpp  -  GLAST LAT Response class Edisp methods        *
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
 * @file GLATEdisp.cpp
 * @brief GLATResponse class energy dispersion implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GTools.hpp"
#include "GException.hpp"
#include "GLATResponse.hpp"
#include "GFitsBinTable.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INIT_EDISP                             "GLATResponse::init_edisp()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return energy dispersion (units: MeV^-1).
 *
 * @param[in] obsEng Observed energy of photon.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 *
 * @todo Needs to be implemented.
 ***************************************************************************/
double GLATResponse::edisp(const GEnergy& obsEng,
                           const GSkyDir& srcDir, const GEnergy& srcEng,
                           const GTime& srcTime, const GPointing& pnt) const
{
    // Return Edisp value
    return 1.0;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise energy dispersion
 ***************************************************************************/
void GLATResponse::edisp_init(void)
{
    // Build filename
    std::string filename = "edisp_"  + m_rspname + "_" + m_rsptype + ".fits";

    // Open FITS file
    GFits file;
    file.open(m_caldb + "/" + filename);

    // Get pointer towards energy dispersion HDU
    GFitsTable* hdu = file.table("ENERGY DISPERSION");
    if (hdu == NULL)
        throw GException::fits_hdu_not_found(G_INIT_EDISP, "ENERGY DISPERSION");

    // Get the energy and cos(theta) bins
    m_edisp_bins.read(hdu);

    // Get the data
    GVector dnorm  = get_fits_vector(hdu, "DNORM");
    GVector ltail  = get_fits_vector(hdu, "LTAIL");
    GVector rwidth = get_fits_vector(hdu, "RWIDTH");
    GVector nr2    = get_fits_vector(hdu, "NR2");
    GVector lt2    = get_fits_vector(hdu, "LT2");
    GVector rt2    = get_fits_vector(hdu, "RT2");

    // ...

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Edisp HDUs tu FOTS file
 *
 * @param[in] file FITS file into which the Edisp HDUs will be appended
 *
 * Append ? HDUs to the FITS file:
 * EBOUNDS (Edisp energy and zenith angle boundaries)
 * ...
 ***************************************************************************/
void GLATResponse::edisp_append(GFits& file) const
{
    // Get Edisp boundary table
    GFitsBinTable* hdu_bounds = new GFitsBinTable;
    hdu_bounds->extname("EBOUNDS");
    m_edisp_bins.write(hdu_bounds);

    // Append HDUs to FITS file
    file.append(hdu_bounds);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATResponse::edisp_init_members(void)
{
    // Initialise Edisp members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 ***************************************************************************/
void GLATResponse::edisp_copy_members(const GLATResponse& rsp)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
***************************************************************************/
void GLATResponse::edisp_free_members(void)
{
    // Return
    return;
}
