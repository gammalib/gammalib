/***************************************************************************
 *                  GCTAAeff.i - CTA effective area base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCTAAeff.i
 * @brief CTA effective area base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAAeff.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAAeff
 *
 * @brief Abstract base class for the CTA effective area
 ***************************************************************************/
class GCTAAeff : public GBase {

public:
    // Constructors and destructors
    GCTAAeff(void);
    GCTAAeff(const GCTAAeff& cta);
    virtual ~GCTAAeff(void);

    // Pure virtual operators
    virtual double operator()(const double& logE, 
                              const double& theta = 0.0, 
                              const double& phi = 0.0,
                              const double& zenith = 0.0,
                              const double& azimuth = 0.0,
                              const bool&   etrue = true) const = 0;

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GCTAAeff*   clone(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual std::string filename(void) const = 0;
};


/***********************************************************************//**
 * @brief GCTAAeff class extension
 ***************************************************************************/
%extend GCTAAeff {
    char *__str__() {
        return tochar(self->print());
    }
};
