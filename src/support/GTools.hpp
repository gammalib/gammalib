/***************************************************************************
 *                       GTools.hpp  -  GammaLib tools                     *
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

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */
const double pi          =  3.141592653589793238462643383279502884197;
const double twopi       =  6.283185307179586476925286766559005768394;
const double fourpi      = 12.56637061435917295385057353311801153679;
const double pihalf      =  1.570796326794896619231321691639751442099;
const double inv_pihalf  =  0.6366197723675813430755350534900574;
const double inv_sqrt4pi =  0.2820947917738781434740397257803862929220;
const double deg2rad     =  0.0174532925199432954743717;
const double rad2deg     = 57.295779513082322864647722;
const double ln2         =  0.6931471805599453094172321214581766;
const double ln10        =  2.3025850929940456840179914546843642;
const double inv_ln2     =  1.4426950408889634073599246810018921;
const double onethird    =  1.0/3.0;
const double twothird    =  2.0/3.0;
const double fourthird   =  4.0/3.0;

/* __ Prototypes ________________________________________________________ */
std::string strip_whitespace(const std::string& arg);
std::string str(const int& value);
std::string toupper(const std::string& s);
std::string tolower(const std::string& s);
double      modulo(double v1, double v2);

#endif /* GTOOLS_HPP */
