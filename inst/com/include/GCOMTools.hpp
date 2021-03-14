/***************************************************************************
 *                     GCOMTools.hpp - COMPTEL tools                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GCOMTools.hpp
 * @brief Definition of COMPTEL tools
 * @author Juergen Knoedlseder
 */

#ifndef GCOMTOOLS_HPP
#define GCOMTOOLS_HPP

/* __ Includes ___________________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GTime;

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    GTime com_time(const int& tjd, const int& tics);
    int   com_tjd(const GTime& time);
    int   com_tics(const GTime& time);
}

#endif /* GCOMTOOLS_HPP */
