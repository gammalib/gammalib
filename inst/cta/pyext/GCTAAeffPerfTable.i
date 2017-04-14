/***************************************************************************
 *     GCTAAeffPerfTable.i - CTA performance table effective area class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
 * @file GCTAAeffPerfTable.hpp
 * @brief CTA performance table effective area class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAAeffPerfTable.hpp"
%}


/***********************************************************************//**
 * @class GCTAAeffPerfTable
 *
 * @brief CTA performance table effective area class
 ***************************************************************************/
class GCTAAeffPerfTable : public GCTAAeff {

public:
    // Constructors and destructors
    GCTAAeffPerfTable(void);
    explicit GCTAAeffPerfTable(const GFilename& filename);
    GCTAAeffPerfTable(const GCTAAeffPerfTable& cta);
    virtual ~GCTAAeffPerfTable(void);

    // Operators
    double operator()(const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void               clear(void);
    GCTAAeffPerfTable* clone(void) const;
    std::string        classname(void) const;
    void               load(const GFilename& filename);
    GFilename          filename(void) const;
    double             max(const double& logE,
                           const double& zenith,
                           const double& azimuth,
                           const bool&   etrue = true) const;
    GEbounds           ebounds(void) const;

    // Methods
    int            size(void) const;
    void           sigma(const double& sigma);
    const double&  sigma(void) const;
};


/***********************************************************************//**
 * @brief GCTAAeffPerfTable class extension
 ***************************************************************************/
%extend GCTAAeffPerfTable {
    GCTAAeffPerfTable copy() {
        return (*self);
    }
};
