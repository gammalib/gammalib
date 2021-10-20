/***************************************************************************
 *               GLATEventAtom.i - Fermi/LAT event atom class              *
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
 * @file GLATEventAtom.i
 * @brief Fermi/LAT event atom class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATEventAtom.hpp"
%}


/***********************************************************************//**
 * @class GLATEventAtom
 *
 * @brief Fermi/LAT event atom class
 ***************************************************************************/
class GLATEventAtom : public GEventAtom {
public:
    // Constructors and destructors
    GLATEventAtom(void);
    GLATEventAtom(const GLATEventAtom& atom);
    virtual ~GLATEventAtom(void);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GLATEventAtom*       clone(void) const;
    virtual std::string          classname(void) const;
    virtual const GLATInstDir&   dir(void) const;
    virtual const GEnergy&       energy(void) const;
    virtual const GTime&         time(void) const;
    virtual const GPolarization& polarization(void) const;
};


/***********************************************************************//**
 * @brief GLATEventAtom class extension
 ***************************************************************************/
%extend GLATEventAtom {
    GLATEventAtom copy() {
        return (*self);
    }
};
