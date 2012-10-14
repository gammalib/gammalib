/***************************************************************************
 *         GSkyPixel.i  -  2D sky pixel index class SWIG definition        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GSkyPixel.i
 * @brief Sky pixel class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyPixel.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GSkyPixel
 *
 * @brief GSkyPixel class interface defintion
 ***************************************************************************/
class GSkyPixel {
public:
    // Constructors and destructors
    GSkyPixel(void);
    explicit GSkyPixel(double x, double y);
    GSkyPixel(const GSkyPixel& pixel);
    virtual ~GSkyPixel(void);

    // Methods
    void       clear(void);
    GSkyPixel* clone(void) const;
    void       x(const double& x);
    void       y(const double& y);
    double     x(void) const;
    double     y(void) const;
};


/***********************************************************************//**
 * @brief GSkyPixel class extension
 ***************************************************************************/
%extend GSkyPixel {
    char *__str__() {
        return tochar(self->print());
    }
    GSkyPixel copy() {
        return (*self);
    }
};
