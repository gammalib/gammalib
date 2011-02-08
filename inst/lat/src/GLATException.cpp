/***************************************************************************
 *                 GLATException.cpp  - LAT exception handler              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATException.cpp
 * @brief LAT exception handler interface implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GLATException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Exception message.
 ***************************************************************************/
const char* GLATExceptionHandler::what() const throw()
{
    std::string message = "*** ERROR in " + m_origin + ": " + m_message;
    return message.c_str();
}


/***********************************************************************//**
 * @brief Error while opening file.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Name of file that could not be opened.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::file_open_error::file_open_error(std::string origin,
                                                std::string filename,
                                                std::string message)
{
    m_origin  = origin;
    m_message = "Unable to open file '"+filename+"'. "+message;
    return;
}


/***********************************************************************//**
 * @brief No sky pixels found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::no_sky::no_sky(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No sky pixels have been found. "+message;
    return;
}


/***********************************************************************//**
 * @brief No energy boundary information found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::no_ebds::no_ebds(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No energy boundaries have been found. "+message;
    return;
}


/***********************************************************************//**
 * @brief No Good Time Intervals found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::no_gti::no_gti(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No Good Time Intervals (GTIs) have been found. "+message;
    return;
}


/***********************************************************************//**
 * @brief No Livetime Cube found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::no_ltcube::no_ltcube(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No Livetime Cube have been found. "+message;
    return;
}


/***********************************************************************//**
 * @brief No energies set
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::no_energies::no_energies(std::string origin, std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Energy vector has not been setup."
                " Cannot access event information.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief No sky directions set
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::no_dirs::no_dirs(std::string origin, std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Sky direction vector has not been setup."
                " Cannot access event information.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief ROI is not a LAT ROI
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::bad_roi_type::bad_roi_type(std::string origin, std::string message)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "Specified ROI is not of type GLATRoi.";
    if (message.length() > 0)
        m_message += " "+message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid response found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::invalid_response::invalid_response(std::string origin,
                                                  std::string message)
{
    m_origin  = origin;
    m_message = "Invalid response encountered. " + message;
    return;
}


/***********************************************************************//**
 * @brief Incompatible source map
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] name Map name.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::wcs_incompatible::wcs_incompatible(std::string origin,
                                                  std::string name,
                                                  std::string message)
{
    m_origin  = origin;
    m_message = "Source map \""+name+"\" incompatible with counts map. " +
                message;
    return;
}


/***********************************************************************//**
 * @brief Diffuse model not found
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] name Diffuse model.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::diffuse_not_found::diffuse_not_found(std::string origin,
                                                    std::string name,
                                                    std::string message)
{
    m_origin  = origin;
    m_message = "Diffuse model \""+name+"\" not found. " +
                message;
    return;
}


/***********************************************************************//**
 * @brief Inconsistent response table
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] size Size.
 * @param[in] expect Expected size.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::inconsistent_response::inconsistent_response(std::string origin,
                                                            int         size,
                                                            int         expect,
                                                            std::string message)
{
    m_origin  = origin;
    m_message = "Inconsistent response table found. Expected "+str(expect)+
                " elements, found "+str(size)+". " + message;
    return;
}
