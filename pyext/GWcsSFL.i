/***************************************************************************
 *           GWcsSFL.i - Sanson-Flamsteed (SFL) projection class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GWcsSFL.i
 * @brief Sanson-Flamsteed (SFL) projection class Python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcsSFL.hpp"
%}


/***********************************************************************//**
 * @class GWcsSFL
 *
 * @brief Sanson-Flamsteed (SFL) projection class definition
 ***************************************************************************/
class GWcsSFL : public GWcs {
public:
    // Constructors and destructors
    GWcsSFL(void);
    GWcsSFL(const std::string& coords,
            const double& crval1, const double& crval2,
            const double& crpix1, const double& crpix2,
            const double& cdelt1, const double& cdelt2);
    GWcsSFL(const GWcsSFL& wcs);
    virtual ~GWcsSFL(void);

    // Implemented pure virtual methods
    virtual void        clear(void);
    virtual GWcsSFL*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
};


/***********************************************************************//**
 * @brief GWcsSFL class extension
 ***************************************************************************/
%extend GWcsSFL {
    GWcsSFL copy() {
        return (*self);
    }
};
