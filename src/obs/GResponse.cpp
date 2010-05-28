/***************************************************************************
 *               GResponse.cpp  -  Response abstract base class            *
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

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <iostream>
#include <string>
#include <unistd.h>           // access() function
#include "GException.hpp"
#include "GResponse.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CALDB                              "GResponse::caldb(std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GResponse::GResponse(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp Response from which the instance should be built.
 ***************************************************************************/
GResponse::GResponse(const GResponse& rsp)
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
GResponse::~GResponse(void)
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
 * @param[in] rsp Response which should be assigned.
 ***************************************************************************/
GResponse& GResponse::operator= (const GResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

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
 * @param[in] obsDir Pointer to observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 *
 * This method implements the default and complete instrument response
 * function (IRF). It may be overwritted by a specific method in the derived
 * class that drops response terms that are not used.
 ***************************************************************************/
double GResponse::irf(const GInstDir& obsDir, const GEnergy& obsEng,
                      const GTime& obsTime,
                      const GSkyDir&  srcDir, const GEnergy& srcEng,
                      const GTime& srcTime, const GPointing& pnt) const
{
    // Get IRF components
    double irf;
    irf  =   psf(obsDir, obsEng, obsTime, srcDir, srcEng, srcTime, pnt);
    irf *=  aeff(obsDir, obsEng, obsTime, srcDir, srcEng, srcTime, pnt);
    irf *= edisp(obsDir, obsEng, obsTime, srcDir, srcEng, srcTime, pnt);
    irf *= tdisp(obsDir, obsEng, obsTime, srcDir, srcEng, srcTime, pnt);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Set the path to the calibration database.
 *
 * @param[in] caldb Absolute path to calibration database
 *
 * This default method simply checks if the calibration database directory
 * exists. If the directory exists, the path will be stored. No checking is
 * implemented that checks for the consistency of the calibration database.
 ***************************************************************************/
void GResponse::caldb(const std::string& caldb)
{
    // Check if calibration database directory is accessible
    if (access(caldb.c_str(), R_OK) != 0)
        throw GException::caldb_not_found(G_CALDB, caldb);
    
    // Store the path to the calibration database
    m_caldb = caldb;

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
void GResponse::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_rspname.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response members which should be copied.
 ***************************************************************************/
void GResponse::copy_members(const GResponse& rsp)
{
    // Copy attributes
    m_caldb   = rsp.m_caldb;
    m_rspname = rsp.m_rspname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GResponse::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
