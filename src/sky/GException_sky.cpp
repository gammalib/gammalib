/***************************************************************************
 *               GException_sky.cpp  -  sky exception handlers             *
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

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief General GSkymap error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Error message.
 ***************************************************************************/
GException::skymap::skymap(std::string origin,
                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Bad sky map parameter
 *
 * @param[in] origin Method that throws the error.
 * @param[in] par Map parameter.
 * @param[in] message Error message.
 ***************************************************************************/
GException::skymap_bad_par::skymap_bad_par(std::string origin,
                                           int         par,
                                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid map parameter ("+str(par)+"). " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Incompatible data size
 *
 * @param[in] origin Method that throws the error.
 * @param[in] size Encountered data size.
 * @param[in] expected Expected data size.
 * @param[in] message Error message.
 ***************************************************************************/
GException::skymap_bad_size::skymap_bad_size(std::string origin,
                                             int         size,
                                             int         expected,
                                             std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Data size ("+str(size)+") is incompatible with "
                "expected size ("+str(expected)+"). " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief General GWcs error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Error message.
 ***************************************************************************/
GException::wcs::wcs(std::string origin,
                     std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid WCS
 *
 * @param[in] origin Method that throws the error.
 * @param[in] wcs WCS.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::wcs_invalid::wcs_invalid(std::string origin,
                                     std::string wcs,
                                     std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid WCS ("+wcs+"). " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Coordinate system invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] coordsys Specified coordinate system.
 ***************************************************************************/
GException::wcs_bad_coords::wcs_bad_coords(std::string origin,
                                           std::string coordsys)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid WCS coordinate system ("+coordsys+"). "
                "Should be one of 'EQU', 'CEL' or 'GAL'.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief No projection function declared
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::wcs_no_proj_fct::wcs_no_proj_fct(std::string origin,
                                             std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "No projection function declared. " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix resolution invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] nside Specified nside parameter.
 ***************************************************************************/
GException::wcs_hpx_bad_nside::wcs_hpx_bad_nside(std::string origin,
                                                 int         nside)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid Healpix nside value ("+str(nside)+"). "
                "Should be one of 1,2,4,8,16,32,64,128,256,512,1024,2048,"
                "4096,8192.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix pixel ordering invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] ordering Specified pixel ordering.
 ***************************************************************************/
GException::wcs_hpx_bad_ordering::wcs_hpx_bad_ordering(std::string origin,
                                                       std::string ordering)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Invalid pixel ordering ("+ordering+"). "
                "Should be either 'RING' or 'NESTED'.";

    // Return
    return;
}
