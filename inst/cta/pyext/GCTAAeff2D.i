/***************************************************************************
 *                  GCTAAeff2D.i - CTA 2D effective area class             *
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
 * @file GCTAAeff2D.hpp
 * @brief CTA 2D effective area class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAAeff2D.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAAeff2D
 *
 * @brief CTA 2D effective area class
 ***************************************************************************/
class GCTAAeff2D : public GCTAAeff {

public:
    // Constructors and destructors
    GCTAAeff2D(void);
    explicit GCTAAeff2D(const std::string& filename);
    GCTAAeff2D(const GCTAAeff2D& cta);
    virtual ~GCTAAeff2D(void);

    // Operators
    double operator()(const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void        clear(void);
    GCTAAeff2D* clone(void) const;
    void        load(const std::string& filename);
    std::string filename(void) const;

    // Methods
    void read(const GFits& file);
};


/***********************************************************************//**
 * @brief GCTAAeff2D class extension
 ***************************************************************************/
%extend GCTAAeff2D {
    GCTAAeff2D copy() {
        return (*self);
    }
};
