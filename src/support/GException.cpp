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

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"


/***************************************************************************
 *                       Private conversion functions                      *
 ***************************************************************************/
std::string str(int value)
{
    std::ostringstream s_value;
    s_value << value;
    return  s_value.str();
}
std::string str(double value)
{
    std::ostringstream s_value;
    s_value << std::scientific << value;
    return  s_value.str();
}


/***************************************************************************
 *                  Exception handler base class definition                *
 ***************************************************************************/
const char* GExceptionHandler::what() const throw()
{
    std::string message = "*** ERROR in " + m_origin + ": " + m_message;
    return message.c_str();
}


/***************************************************************************
 *                        Memory allocation exception                      *
 ***************************************************************************/
GException::mem_alloc::mem_alloc(std::string origin, unsigned num)
{
    m_origin  = origin;
    m_message = "Memory allocation error (" + str((int)num) + " elements)";
}


