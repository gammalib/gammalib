/***************************************************************************
 *               GWcsMER.i - Mercator's (MER) projection class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GWcsMER.i
 * @brief Mercator's (MER) projection class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcsMER.hpp"
%}


/***********************************************************************//**
 * @class GWcsMER
 *
 * @brief Mercator's (MER) projection class definition
 ***************************************************************************/
class GWcsMER : public GWcs {
public:
    // Constructors and destructors
    GWcsMER(void);
    explicit GWcsMER(const std::string& coords,
                     const double& crval1, const double& crval2,
                     const double& crpix1, const double& crpix2,
                     const double& cdelt1, const double& cdelt2);
    GWcsMER(const GWcsMER& wcs);
    virtual ~GWcsMER(void);

    // Implemented pure virtual methods
    virtual void        clear(void);
    virtual GWcsMER*    clone(void) const;
    virtual std::string code(void) const;
    virtual std::string name(void) const;
};


/***********************************************************************//**
 * @brief GWcsMER class extension
 ***************************************************************************/
%extend GWcsMER {
    GWcsMER copy() {
        return (*self);
    }
};
