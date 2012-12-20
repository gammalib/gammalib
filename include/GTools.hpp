/***************************************************************************
 *                        GTools.hpp - GammaLib tools                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GTools.hpp
 * @brief Gammalib tools definition
 * @author Juergen Knoedlseder
 */

#ifndef GTOOLS_HPP
#define GTOOLS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <cmath>

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */
const double pi           =  3.141592653589793238462643383279502884197;
const double twopi        =  6.283185307179586476925286766559005768394;
const double fourpi       = 12.56637061435917295385057353311801153679;
const double pihalf       =  1.570796326794896619231321691639751442099;
const double inv_pihalf   =  0.6366197723675813430755350534900574;
const double inv_sqrt4pi  =  0.2820947917738781434740397257803862929220;
const double pi2          = pi*pi;
const double deg2rad      =  0.0174532925199432954743717;
const double rad2deg      = 57.295779513082322864647722;
const double ln2          =  0.6931471805599453094172321214581766;
const double ln10         =  2.3025850929940456840179914546843642;
const double inv_ln2      =  1.4426950408889634073599246810018921;
const double onethird     =  1.0/3.0;
const double twothird     =  2.0/3.0;
const double fourthird    =  4.0/3.0;
const double sqrt_onehalf = sqrt(1.0/2.0);
const double sqrt_pihalf  = sqrt(pihalf);
const double sqrt_two     = sqrt(2.0);
const double MeV2erg      =  1.6021765e-6;
const double erg2MeV      =  624150.96;
const double pc2cm        =  3.08568025e18;
const double sec_in_day   = 86400.0;

/* __ Prototypes ________________________________________________________ */
std::string              strip_whitespace(const std::string& arg);
std::string              strip_chars(const std::string& arg,
                                     const std::string& chars);
std::string              expand_env(const std::string& arg);
std::string              str(const unsigned short int& value);
std::string              str(const unsigned int& value);
std::string              str(const unsigned long int& value);
std::string              str(const unsigned long long int& value);
std::string              str(const short int& value);
std::string              str(const int& value);
std::string              str(const long int& value);
std::string              str(const long long int& value);
std::string              str(const float& value);
std::string              str(const double& value);
char*                    tochar(const std::string& arg);
short                    toshort(const std::string& arg);
unsigned short           toushort(const std::string& arg);
int                      toint(const std::string& arg);
unsigned int             touint(const std::string& arg);
long                     tolong(const std::string& arg);
unsigned long            toulong(const std::string& arg);
long long                tolonglong(const std::string& arg);
unsigned long long       toulonglong(const std::string& arg);
float                    tofloat(const std::string& arg);
double                   todouble(const std::string& arg);
std::string              toupper(const std::string& s);
std::string              tolower(const std::string& s);
std::vector<std::string> split(const std::string& s, const std::string& sep);
std::string              fill(const std::string& s, int n);
std::string              left(const std::string& s, int n, char c = ' ');
std::string              right(const std::string& s, int n, char c = ' ');
std::string              center(const std::string& s, int n, char c = ' ');
std::string              parformat(const std::string& s, const int& indent = 0);
double                   modulo(double v1, double v2);
double                   arccos(const double& arg);
double                   plaw_photon_flux(const double& emin,
                                          const double& emax,
                                          const double& epivot,
                                          const double& gamma);
double                   plaw_energy_flux(const double& emin,
                                          const double& emax,
                                          const double& epivot,
                                          const double& gamma);
bool                     file_exists(const std::string& filename);
bool                     isinfinite(double x);
bool                     isnotanumber(double x);

#endif /* GTOOLS_HPP */
