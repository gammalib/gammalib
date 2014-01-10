/***************************************************************************
 *               GLATEventBin.i - Fermi/LAT event bin class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GLATEventBin.i
 * @brief Fermi/LAT event bin class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATEventBin.hpp"
%}


/***********************************************************************//**
 * @class GLATEventBin
 *
 * @brief Fermi/LAT event bin class
 ***************************************************************************/
class GLATEventBin : public GEventBin {

    // Friend classes
    friend class GLATEventCube;

public:
    // Constructors and destructors
    GLATEventBin(void);
    GLATEventBin(const GLATEventBin& bin);
    virtual ~GLATEventBin(void);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GLATEventBin*      clone(void) const;
    virtual double             size(void) const;
    virtual const GLATInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);

    // Other methods
    const double&        solidangle(void) const;
    const GEnergy&       ewidth(void) const;
    const double&        ontime(void) const;
    const int&           index(void) const;
    const int&           ipix(void) const;
    const int&           ieng(void) const;
    const GLATEventCube* cube(void) const;
};


/***********************************************************************//**
 * @brief GLATEventBin class extension
 ***************************************************************************/
%extend GLATEventBin {
    GLATEventBin copy() {
        return (*self);
    }
};
