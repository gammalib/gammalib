/***************************************************************************
 *                 GBilinear.i - Bilinear interpolator class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Juergen Knoedlseder                              *
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
 * @file GBilinear.i
 * @brief Bilinear interpolator class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GBilinear.hpp"
%}


/***********************************************************************//**
 * @class GBilinear
 *
 * @brief Bilinear interpolator class
 ***************************************************************************/
class GBilinear : public GBase {

public:
    // Constructors and destructors
    GBilinear(void);
    GBilinear(const GBilinear& interpolator);
    virtual ~GBilinear(void);

    // Operators
    double operator()(const double* array);

    // Methods
    void        clear(void);
    GBilinear*  clone(void) const;
    std::string classname(void) const;
    int&        index1(void);
    int&        index2(void);
    int&        index3(void);
    int&        index4(void);
    double&     weight1(void);
    double&     weight2(void);
    double&     weight3(void);
    double&     weight4(void);
};


/***********************************************************************//**
 * @brief GBilinear class extension
 ***************************************************************************/
%extend GBilinear {
    GBilinear copy() {
        return (*self);
    }
};
