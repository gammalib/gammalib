/***************************************************************************
 *      GFunctions.i - Single parameter functions abstract base class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GFunctions.i
 * @brief Single parameter functions abstract base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFunctions.hpp"
#include "GVector.hpp"
%}


/***********************************************************************//**
 * @class GFunctions
 *
 * @brief Single parameter vector function abstract base class
 ***************************************************************************/
class GFunctions {
public:
    // Constructors and destructors
    GFunctions(void);
    GFunctions(const GFunctions& functions);
    virtual ~GFunction(void);

    // Methods
    virtual int     size(void) const = 0;
    virtual GVector eval(const double& x) = 0;
};
