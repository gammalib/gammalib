/***************************************************************************
 *          GLATAeff.cpp  -  GLAST LAT Response class Aeff methods         *
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
 * @file GLATAeff.cpp
 * @brief GLATResponse class effective area implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GTools.hpp"
#include "GException.hpp"
#include "GLATResponse.hpp"
#include "GFitsDblImage.hpp"
#include "GFitsTableFltCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_INIT_AEFF "GLATResponse::init_aeff()"

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
 * @brief Return effective area
 *
 * @param[in] obsDir Observed photon direction
 * @param[in] obsEng Observed energy of photon
 * @param[in] srcDir True photon direction
 * @param[in] srcEng True energy of photon
 * @param[in] instPntDir Instrument pointing direction (e.g. z-axis)
 * @param[in] instPosAng Instrument position angle
 * @param[in] time Photon arrival time
 ***************************************************************************/
double GLATResponse::aeff(const GSkyDir& obsDir, const double& obsEng,
                          const GSkyDir& srcDir, const double& srcEng,
                          const GSkyDir& instPntDir, const double& instPosAng,
                          const double& time)
{
    // Return Aeff value
    return 0.0; //!< @todo Not yet implemented
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise effective area
 ***************************************************************************/
void GLATResponse::aeff_init(void)
{
    // Build filename
    std::string filename = "aeff_"  + m_rspname + "_" + m_rsptype + ".fits";

    // Open FITS file
    GFits file;
    file.open(m_caldb + "/" + filename);

    // Get pointer to effective area HDU
    GFitsHDU* hdu = file.hdu("EFFECTIVE AREA");
    if (hdu == NULL)
        throw GException::fits_hdu_not_found(G_INIT_AEFF, "EFFECTIVE AREA");

    // Get the energy and cos(theta) bins
    m_aeff_bins.load(hdu);

    // Get the data
    GVector effarea = get_fits_vector(hdu, "EFFAREA");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Aeff HDUs tu FOTS file
 *
 * @param[in] file FITS file into which the Aeff HDUs will be appended
 *
 * Append 2 HDUs to the FITS file:
 * ABOUNDS (Aeff energy and zenith angle boundaries)
 * AEFF (Aeff array)
 ***************************************************************************/
void GLATResponse::aeff_append(GFits& file) const
{
    // Get Aeff boundary table
    GFitsHDU hdu_bounds;
    m_aeff_bins.save(&hdu_bounds);
    hdu_bounds.extname("ABOUNDS");

    // Append HDUs to FITS file
    file.append_hdu(hdu_bounds);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATResponse::aeff_init_members(void)
{
    // Initialise Aeff members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 ***************************************************************************/
void GLATResponse::aeff_copy_members(const GLATResponse& rsp)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
***************************************************************************/
void GLATResponse::aeff_free_members(void)
{
    // Return
    return;
}
