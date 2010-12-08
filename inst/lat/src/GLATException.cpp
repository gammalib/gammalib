/***************************************************************************
 *                 GLATException.cpp  - LAT exception handler              *
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
    static std::string message = "*** ERROR in " + m_origin + ": " + m_message;
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
 * @brief Instrument response not set.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GLATException::no_response::no_response(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "Attempt to model LAT response, but no response function"
                " has been found. Use GLATObservation::response() method"
                " to set response function. "+message;
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
