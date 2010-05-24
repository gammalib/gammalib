/***************************************************************************
 *             GCTAAeff.cpp  -  CTA Response class Aeff methods            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAAeff.cpp
 * @brief GCTAResponse class effective area implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <iostream>
//#include "GTools.hpp"
//#include "GException.hpp"
#include "GCTAResponse.hpp"
//#include "GFitsImageDbl.hpp"
//#include "GFitsTableFltCol.hpp"

/* __ Method name definitions ____________________________________________ */

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
 * @brief Return effective area un units of cm2
 *
 * @param[in] obsDir Observed photon direction
 * @param[in] obsEng Observed energy of photon
 * @param[in] srcDir True photon direction
 * @param[in] srcEng True energy of photon
 * @param[in] instPntDir Instrument pointing direction (e.g. z-axis)
 * @param[in] instPosAng Instrument position angle
 * @param[in] time Photon arrival time
 *
 * @todo Implement method (just a dummy for the moment)
 ***************************************************************************/
double GCTAResponse::aeff(const GSkyDir& obsDir, const GEnergy& obsEng,
                          const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GSkyDir& instPntDir, const double& instPosAng,
                          const GTime& time)
{
    // Return Aeff value
    return 1.0;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/
