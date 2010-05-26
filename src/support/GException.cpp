/***************************************************************************
 *                   GException.cpp  -  exception handlers                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GException.cpp
 * @brief Exception handler interface implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Exception message.
 ***************************************************************************/
const char* GExceptionHandler::what() const throw()
{
    std::string message = "*** ERROR in " + m_origin + ": " + m_message;
    return message.c_str();
}


/***********************************************************************//**
 * @brief Memory allocation exception.
 *
 * @param[in] origin Name of method that has thrown the exception.
 * @param[in] num Requested amount of memory.
 ***************************************************************************/
GException::mem_alloc::mem_alloc(std::string origin, unsigned num)
{
    m_origin  = origin;
    m_message = "Memory allocation error (" + str((int)num) + " elements)";
}
