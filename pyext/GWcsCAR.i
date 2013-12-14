/***************************************************************************
 *            GWcsCAR.i  -  Plate carree (CAR) projection class            *
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
 * @file GWcsCAR.i
 * @brief Plate carree (CAR) projection class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcsCAR.hpp"
%}


/***********************************************************************//**
 * @class GWcsCAR
 *
 * @brief Plate carree (CAR) projection class definition
 ***************************************************************************/
class GWcsCAR : public GWcs {
public:
    // Constructors and destructors
    GWcsCAR(void);
    explicit GWcsCAR(const std::string& coords,
                     const double& crval1, const double& crval2,
                     const double& crpix1, const double& crpix2,
                     const double& cdelt1, const double& cdelt2);
    GWcsCAR(const GWcsCAR& wcs);
    virtual ~GWcsCAR(void);

    // Implemented pure virtual methods
    virtual void        clear(void);
    virtual GWcsCAR*    clone(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
};


/***********************************************************************//**
 * @brief GWcsCAR class extension
 ***************************************************************************/
%extend GWcsCAR {
    GWcsCAR copy() {
        return (*self);
    }
};
