/***************************************************************************
 *                 GCTAException.cpp  - CTA exception handler              *
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
 * @file GCTAException.cpp
 * @brief CTA exception handler interface implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GCTAException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Exception message.
 ***************************************************************************/
const char* GCTAExceptionHandler::what() const throw()
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
GCTAException::file_open_error::file_open_error(std::string origin,
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
GCTAException::no_response::no_response(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "Attempt to model CTA response, but no response function"
                " has been found. Use GCTAObservation::response() method"
                " to set response function. "+message;
    return;
}


/***********************************************************************//**
 * @brief No sky pixels found.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GCTAException::no_sky::no_sky(std::string origin, std::string message)
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
GCTAException::no_ebds::no_ebds(std::string origin, std::string message)
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
GCTAException::no_gti::no_gti(std::string origin, std::string message)
{
    m_origin  = origin;
    m_message = "No Good Time Intervals (GTIs) have been found. "+message;
    return;
}
