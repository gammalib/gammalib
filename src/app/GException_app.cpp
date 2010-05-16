/***************************************************************************
 *           GException_app.cpp  -  Application exception handlers         *
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
 * @brief Parameter file not found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename Filename that was not found.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::par_file_not_found::par_file_not_found(std::string origin,
                                                   std::string filename,
                                                   std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Parameter file '"+filename+"' not found. "
                "Make sure that PFILES environment variable is set "
                "correctly. " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unable to open parameter file
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename Filename that could not be opened.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::par_file_open_error::par_file_open_error(std::string origin,
                                                     std::string filename,
                                                     std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Unable to open parameter file '"+filename+"'. " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unable to determine users home directory
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::home_not_found::home_not_found(std::string origin,
                                           std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Unable to determine users home directory. " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unable to create pfiles directory in users home directory
 *
 * @param[in] origin Method that throws the error.
 * @param[in] home Users home directory.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::could_not_create_pfiles::could_not_create_pfiles(std::string origin,
                                                             std::string home,
                                                             std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Unable to create '"+home+"'. " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unable to change access rights for pfiles directory
 *
 * @param[in] origin Method that throws the error.
 * @param[in] home Users home directory.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::pfiles_not_accessible::pfiles_not_accessible(std::string origin,
                                                         std::string home,
                                                         std::string message)
{
    // Set origin and message
    m_origin  = origin;
    m_message = "Could not make '"+home+"' write accessible for "
                "writing of the applications parameter file. " + message;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Syntax error encountered in parameter file line
 *
 * @param[in] origin Method that throws the error.
 * @param[in] line Parameter file line in which syntax error occured.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::par_file_syntax_error::par_file_syntax_error(std::string origin,
                                                         std::string line,
                                                         std::string message)
{
    // Set origin and message
    m_origin  = origin;
    if (message.length() > 0)
        m_message = "Syntax error occured in the following line of the "
                    "parameter file ("+message+"): "+line;
    else
        m_message = "Syntax error occured in the following line of the "
                    "parameter file: "+line;

    // Return
    return;
}
