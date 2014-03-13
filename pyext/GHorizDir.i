/***************************************************************************
 *                  GHorizDir.i - Horizontal direction class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Karl Kosack                                      *
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
 * @file GHorizDir.i
 * @brief Horizontal direction class interface definition
 * @author Karl Kosack 
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GHorizDir.hpp"
%}


/***********************************************************************//**
 * @class GHorizDir
 *
 * @brief Horizontal (Alt/Az) direction class
 ***************************************************************************/
class GHorizDir : public GBase {
public:
    // Constructors and destructors
    GHorizDir(void);
    GHorizDir(const GHorizDir& dir);
    virtual ~GHorizDir(void);

    // Methods
    void          clear(void);
    GHorizDir*    clone(void) const;
    void          altaz(const double& alt, const double& az);
    void          altaz_deg(const double& alt, const double& az);
    void          celvector(const GVector& vector);
    void          rotate_deg(const double& phi, const double& theta);
    const double& alt(void) const;
    const double& az(void) const;
    double        zenith(void) const;
    double        zenith_deg(void) const;
    double        alt_deg(void) const;
    double        az_deg(void) const;
    GVector       celvector(void) const;
    double        dist(const GHorizDir& dir) const;
    double        dist_deg(const GHorizDir& dir) const;
};


/***********************************************************************//**
 * @brief GHorizDir class extension
 ***************************************************************************/
%extend GHorizDir {
    bool __eq__(const GHorizDir& dir) const {
        return ((*self) == dir);
    }
    bool __ne__(const GHorizDir& dir) const {
        return ((*self) != dir);
    }
    GHorizDir copy() {
        return (*self);
    }
};
