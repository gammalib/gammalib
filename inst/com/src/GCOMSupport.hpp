/***************************************************************************
 *                GCOMSupport.hpp - COMPTEL support functions              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMSupport.hpp
 * @brief Definition of support function used by COMPTEL classes
 * @author Juergen Knoedlseder
 */

#ifndef GCOMSUPPORT_HPP
#define GCOMSUPPORT_HPP

/* __ Includes ___________________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GSkyMap;

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    void          com_wcs_mer2car(GSkyMap& map);
    double        com_energy1(const double& energy, const double& phigeo);
    double        com_energy2(const double& energy, const double& phigeo);
    const double& com_exd2x(const int& id2);
    const double& com_exd2y(const int& id2);
    const double& com_exd2r(const int& id2);
}

#endif /* GCOMSUPPORT_HPP */
