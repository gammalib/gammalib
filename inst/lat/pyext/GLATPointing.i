/***************************************************************************
 *                GLATPointing.i - Fermi/LAT pointing class                *
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
 * @file GLATPointing.i
 * @brief Fermi-LAT pointing class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATPointing.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GLATPointing
 *
 * @brief Interface for the Fermi-LAT pointing
 ***************************************************************************/
class GLATPointing {
public:
    // Constructors and destructors
    GLATPointing(void);
    GLATPointing(const GLATPointing& pnt);
    virtual ~GLATPointing(void);

    // Methods
    void           clear(void);
    GLATPointing*  clone(void) const;
    const GSkyDir& dir(void) const;
    void           dir(const GSkyDir& dir);
};


/***********************************************************************//**
 * @brief GLATPointing class extension
 ***************************************************************************/
%extend GLATPointing {
    GLATPointing copy() {
        return (*self);
    }
};
