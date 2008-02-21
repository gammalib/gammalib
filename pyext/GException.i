/***************************************************************************
 *                     GException.i  -  exception handler                  *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
%module GException
%{
/* Put headers and other declarations here */
#include <string>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception
%}
#include "GException.hpp"
