/***************************************************************************
 *                  GWcsMOL.i - Mollweide's projection class               *
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
 * @file GWcsMOL.i
 * @brief Mollweide's projection class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcsMOL.hpp"
%}


/***********************************************************************//**
 * @class GWcsMOL
 *
 * @brief Mollweide's projection class definition
 ***************************************************************************/
class GWcsMOL : public GWcs {

public:
    // Constructors and destructors
    GWcsMOL(void);
    GWcsMOL(const std::string& coords,
            const double& crval1, const double& crval2,
            const double& crpix1, const double& crpix2,
            const double& cdelt1, const double& cdelt2);
    GWcsMOL(const GWcsMOL& wcs);
    virtual ~GWcsMOL(void);

    // Implemented pure virtual methods
    virtual void        clear(void);
    virtual GWcsMOL*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
};


/***********************************************************************//**
 * @brief GWcsMOL class extension
 ***************************************************************************/
%extend GWcsMOL {
    GWcsMOL copy() {
        return (*self);
    }
};
