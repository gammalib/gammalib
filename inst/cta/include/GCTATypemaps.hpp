/***************************************************************************
 *                 GCTATypemaps.hpp - GammaLib CTA typemaps                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2020 by Juergen Knoedlseder                         *
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
 * @file GCTATypemaps.hpp
 * @brief Definition of GammaLib CTA typemaps
 * @author Juergen Knoedlseder
 */

#ifndef GCTATYPEMAPS_HPP
#define GCTATYPEMAPS_HPP

/* __ Includes ___________________________________________________________ */

/* __ Class code enumerations (used primarily to avoid dynamic casting) __ */
typedef enum {
    GCTA_CUBE_SOURCE_POINT,
    GCTA_CUBE_SOURCE_DIFFUSE
} GCTAClassCode;

/* __ Typemaps ___________________________________________________________ */

/* __ Prototypes _________________________________________________________ */

/* __ CTA specific definitions ___________________________________________ */
#define G_CTA_MJDREF 51544.5                 //!< Reference of CTA time frame

#endif /* GCTATYPEMAPS_HPP */
