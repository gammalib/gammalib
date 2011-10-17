/***************************************************************************
 *                    GCTADir.i  -  CTA direction class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @file GCTADir.i
 * @brief CTA camera direction class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTADir.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTADir
 *
 * @brief CTA camera direction class.
 ***************************************************************************/
class GCTADir {
public:
    // Constructors and destructors
    GCTADir(void);
    GCTADir(const GCTADir& dir);
    explicit GCTADir(const GSkyDir& dir, const GCTAPointing& pnt);
    explicit GCTADir(const GCTAInstDir& dir, const GCTAPointing& pnt);
    virtual ~GCTADir(void);

    // Methods
    void          clear(void);
    GCTADir*      clone(void) const;
    void          dir(const GSkyDir& dir, const GCTAPointing& pnt);
    void          dir(const GCTAInstDir& dir, const GCTAPointing& pnt);
    const double& theta(void) const;
    const double& phi(void) const;
    double        theta_deg(void) const;
    double        phi_deg(void) const;
    const double& costheta(void) const;
    const double& sintheta(void) const;
};


/***********************************************************************//**
 * @brief GCTADir class extension
 ***************************************************************************/
%extend GCTADir {
    char *__str__() {
        return tochar(self->print());
    }
    GCTADir copy() {
        return (*self);
    }
};
