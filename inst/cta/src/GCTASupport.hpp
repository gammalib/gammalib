/***************************************************************************
 *                  GCTASupport.hpp - CTA support functions                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GCTASupport.hpp
 * @brief Definition of support function used by CTA classes
 * @author Juergen Knoedlseder
 */

#ifndef GCTASUPPORT_HPP
#define GCTASUPPORT_HPP

/* __ Includes ___________________________________________________________ */

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Typedefs ___________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    double cta_roi_arclength(const double& rad,     const double& dist,
                             const double& cosdist, const double& sindist,
                             const double& roi,     const double& cosroi);
}

#endif /* GCTASUPPORT_HPP */
