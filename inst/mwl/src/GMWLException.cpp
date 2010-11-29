/***************************************************************************
 *          GMWLException.cpp  - Multi-wavelength exception handler        *
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
 * @file GMWLException.cpp
 * @brief MWL exception handler interface implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GMWLException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Exception message
 ***************************************************************************/
const char* GMWLExceptionHandler::what() const throw()
{
    std::string message = "*** ERROR in " + m_origin + ": " + m_message;
    return message.c_str();
}


/***********************************************************************//**
 * @brief Error while opening file
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] filename Name of file that could not be opened.
 * @param[in] message Optional error message.
 ***************************************************************************/
GMWLException::file_open_error::file_open_error(std::string origin,
                                                std::string filename,
                                                std::string message)
{
    m_origin  = origin;
    m_message = "Unable to open file \""+filename+"\". "+message;
    return;
}


/***********************************************************************//**
 * @brief File has invalid format
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] message Optional error message.
 ***************************************************************************/
GMWLException::bad_file_format::bad_file_format(std::string origin,
                                                std::string message)
{
    m_origin  = origin;
    m_message = "File has invalid format. "+message;
    return;
}


/***********************************************************************//**
 * @brief Invalid or unsupported unit encountered
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] unit Unit string.
 * @param[in] message Optional error message.
 ***************************************************************************/
GMWLException::invalid_unit::invalid_unit(std::string origin,
                                          std::string unit,
                                          std::string message)
{
    m_origin  = origin;
    m_message = "Invalid or unsupported unit \""+unit+"\" encountered. " +
                message;
    return;
}
