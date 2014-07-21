/***************************************************************************
 *                    GTypemaps.hpp - GammaLib typemaps                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2014 by Juergen Knoedlseder                         *
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
 * @file GTypemaps.hpp
 * @brief Definition of GammaLib typemaps
 * @author Juergen Knoedlseder
 */

#ifndef GTYPEMAPS_HPP
#define GTYPEMAPS_HPP

/* __ Includes ___________________________________________________________ */

/* __ Chatter level enumerations _________________________________________ */
typedef enum {
    SILENT = 0,
    TERSE = 1,
    NORMAL = 2,
    EXPLICIT = 3,
    VERBOSE = 4
} GChatter;

/* __ Class code enumerations (used primarily to avoid dynamic casting) __ */
typedef enum {
    GMODEL_SPATIAL_POINT_SOURCE,
    GMODEL_SPATIAL_RADIAL,
    GMODEL_SPATIAL_ELLIPTICAL,
    GMODEL_SPATIAL_DIFFUSE
} GClassCode;

/* __ Typemaps ___________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    GChatter reduce(const GChatter& chatter);
}


/***********************************************************************//**
 * @brief Reduce chattiness by one level
 *
 * @param[in] chatter Chattiness.
 * @return Reduced chattiness.
 ***************************************************************************/
inline
GChatter gammalib::reduce(const GChatter& chatter)
{
    // Allocate reduced chattiness
    GChatter reduced;

    // Reduce chattiness
    switch (chatter) {
        case SILENT:
            reduced = SILENT;
            break;
        case TERSE:
            reduced = SILENT;
            break;
        case NORMAL:
            reduced = TERSE;
            break;
        case EXPLICIT:
            reduced = NORMAL;
            break;
        case VERBOSE:
            reduced = EXPLICIT;
            break;
        default:
            reduced = chatter;
            break;
    }
    
    // Return reduced chattiness
    return reduced;
}

#endif /* GTYPEMAPS_HPP */
