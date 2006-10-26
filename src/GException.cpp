/***************************************************************************
 *                   GException.cpp  -  exception handler                  *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2006 by Jurgen Knodlseder                   *
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
#include <string>                             // string
#include <sstream>                            // ostringstream
using namespace std;


/***************************************************************************
 *                  Exception handler base class definition                *
 ***************************************************************************/
const char* GException::what() const throw()
{
  string message = "*** ERROR in " + m_origin + ": " + m_message;
  return message.c_str();
}


/***************************************************************************
 *                    Specific exception class definitions                 *
 ***************************************************************************/
mem_alloc::mem_alloc(string origin, unsigned num)
{
  ostringstream s_num;
  s_num << num;
  m_origin  = origin;
  m_message = "Memory allocation error (" + s_num.str() + " elements)";
}
empty::empty(string origin)
{
  m_origin  = origin;
  m_message = "Zero-size allocation";
}
