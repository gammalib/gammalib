/***************************************************************************
 *   GCTAEdispPerfTable.i - CTA performance table energy dispersion class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Christoph Deil & Ellis Owen                      *
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
 * @file GCTAEdispPerfTable.i
 * @brief CTA performance table energy dispersion class definition
 * @author Christoph Deil & Ellis Owen
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEdispPerfTable.hpp"
%}


/***********************************************************************//**
 * @class GCTAEdispPerfTable
 *
 * @brief CTA performance energy dispersion function class
 ***************************************************************************/
class GCTAEdispPerfTable : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdispPerfTable(void);
    explicit GCTAEdispPerfTable(const std::string& filename);
    GCTAEdispPerfTable(const GCTAEdispPerfTable& psf);
    virtual ~GCTAEdispPerfTable(void);

    // Operators
    double operator()(const double& logEobs,
                      const double& logEsrc,
                      const double& theta = 0.0,
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0) const;

    // Implemented pure virtual methods
    void                clear(void);
    GCTAEdispPerfTable* clone(void) const;
    void                load(const std::string& filename);
    std::string         filename(void) const;
    GEnergy             mc(GRan&         ran,
                           const double& logE,
                           const double& theta = 0.0,
                           const double& phi = 0.0,
                           const double& zenith = 0.0,
                           const double& azimuth = 0.0) const;
    GEbounds            ebounds_obs(const double& logEsrc,
                                    const double& theta = 0.0,
                                    const double& phi = 0.0,
                                    const double& zenith = 0.0,
                                    const double& azimuth = 0.0) const;
    GEbounds            ebounds_src(const double& logEobs,
                                    const double& theta = 0.0,
                                    const double& phi = 0.0,
                                    const double& zenith = 0.0,
                                    const double& azimuth = 0.0) const;
};


/***********************************************************************//**
 * @brief GCTAEdispPerfTable class extension
 ***************************************************************************/
%extend GCTAEdispPerfTable {
    GCTAEdispPerfTable copy() {
        return (*self);
    }
};
