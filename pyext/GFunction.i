/***************************************************************************
 *       GFunction.i - Single parameter function abstract base class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @file GFunction.i
 * @brief Single parameter function abstract base class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFunction.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GFunction
 *
 * @brief Single parameter function abstract base class
 ***************************************************************************/
class GFunction {
public:
    // Constructors and destructors
    GFunction(void);
    GFunction(const GFunction& func);
    virtual ~GFunction(void);

    // Methods
    virtual double eval(double x) = 0;
};
