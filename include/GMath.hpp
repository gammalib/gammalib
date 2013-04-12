/***************************************************************************
 *                     GMath.hpp - Mathematical functions                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GMath.hpp
 * @brief Mathematical function definitions
 * @author Juergen Knoedlseder
 */

#ifndef GMATH_HPP
#define GMATH_HPP

/* __ Includes ___________________________________________________________ */

/* __ Constants __________________________________________________________ */
namespace gammalib {
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
    const double sqrt_onehalf = std::sqrt(1.0/2.0);
    const double sqrt_pihalf  = std::sqrt(pihalf);
    const double sqrt_two     = std::sqrt(2.0);
}

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    double acos(const double& arg);
    double gammln(const double& arg);
    double modulo(const double& v1, const double& v2);
}

#endif /* GMATH_HPP */
