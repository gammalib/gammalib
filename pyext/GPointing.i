/***************************************************************************
 *                  GPointing.i - Abstract pointing class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GPointing.i
 * @brief Abstract pointing class Python interface definition.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPointing.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GPointing
 *
 * @brief Abstract Pointing class
 ***************************************************************************/
class GPointing : public GBase {

public:
    // Constructors and destructors
    GPointing(void);
    GPointing(const GPointing& pnt);
    virtual ~GPointing(void);

    // Pure virtual methods
    virtual void           clear(void) = 0;
    virtual GPointing*     clone(void) const = 0;
    virtual const GSkyDir& dir(void) const = 0;
};


/***********************************************************************//**
 * @brief GPointing class extension
 ***************************************************************************/
%extend GPointing {
};
