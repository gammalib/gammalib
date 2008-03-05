/***************************************************************************
 *                          GTools.hpp  -  GammaLib tools                  *
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

#ifndef GTOOLS_HPP
#define GTOOLS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
//#include <sstream>                            // ostringstream
//#include <stdexcept>                          // exception

/* __ Namespaces _________________________________________________________ */

/* __ Prototypes ________________________________________________________ */
std::string strip_whitespace(const std::string& arg);
std::string str(const int& value);

#endif /* GTOOLS_HPP */
