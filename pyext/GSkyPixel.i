/***************************************************************************
 *                      GSkyPixel.i - Sky map pixel class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @brief Sky map pixel class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyPixel.hpp"
%}


/***********************************************************************//**
 * @class GSkyPixel
 *
 * @brief GSkyPixel class interface definition
 ***************************************************************************/
class GSkyPixel : public GBase {

public:
    // Constructors and destructors
    GSkyPixel(void);
    GSkyPixel(const int& index);
    GSkyPixel(const double& index);
    GSkyPixel(const int& x, const int& y);
    GSkyPixel(const double& x, const double& y);
    GSkyPixel(const GSkyPixel& pixel);
    virtual ~GSkyPixel(void);

    // Methods
    void          clear(void);
    GSkyPixel*    clone(void) const;
    std::string   classname(void) const;
    int           size(void) const;
    bool          is_1D(void) const;
    bool          is_2D(void) const;
    void          index(const double& index);
    void          x(const double& x);
    void          y(const double& y);
    void          xy(const double& x, const double& y);
    const double& index(void) const;
    const double& x(void) const;
    const double& y(void) const;
};


/***********************************************************************//**
 * @brief GSkyPixel class extension
 ***************************************************************************/
%extend GSkyPixel {
    bool __eq__(const GSkyPixel& pixel) const {
        return ((*self) == pixel);
    }
    bool __ne__(const GSkyPixel& pixel) const {
        return ((*self) != pixel);
    }
    GSkyPixel copy() {
        return (*self);
    }
};
