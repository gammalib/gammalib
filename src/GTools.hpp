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

/* __ Constants __________________________________________________________ */
const double pi      = 3.1415926535897931159979635;
const double pihalf  = pi/2.0;
const double twopi   = 2.0*pi;
const double fourpi  = 4.0*pi;
const double deg2rad = 0.0174532925199432954743717;
const double rad2deg = 57.295779513082322864647722;
const double ln10    = 2.3025850929940459010936138;

/* __ Prototypes ________________________________________________________ */
std::string strip_whitespace(const std::string& arg);
std::string str(const int& value);

#endif /* GTOOLS_HPP */
