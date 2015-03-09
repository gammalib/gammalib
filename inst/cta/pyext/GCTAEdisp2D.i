/***************************************************************************
 *             GCTAEdisp2D.i - CTA 2D energy dispersion class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Florent Forest                                   *
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
 * @file GCTAEdisp2D.i
 * @brief CTA 2D energy dispersion class definition
 * @author Florent Forest
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEdisp2D.hpp"
%}


/***********************************************************************//**
 * @class GCTAEdisp2D
 *
 * @brief CTA 2D energy dispersion class
 ***************************************************************************/
class GCTAEdisp2D : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdisp2D(void);
    explicit GCTAEdisp2D(const std::string& filename);
    GCTAEdisp2D(const GCTAEdisp2D& edisp);
    virtual ~GCTAEdisp2D(void);

    // Operators
    double operator()(const double& logEobs, 
                      const double& logEsrc, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0) const;

    // Implemented methods
    void         clear(void);
    GCTAEdisp2D* clone(void) const;
    std::string  classname(void) const;
    void         load(const std::string& filename);
    std::string  filename(void) const;
    GEnergy      mc(GRan&         ran,
                    const double& logE,
                    const double& theta = 0.0,
                    const double& phi = 0.0,
                    const double& zenith = 0.0,
                    const double& azimuth = 0.0) const;
    GEbounds     ebounds_obs(const double& logEsrc,
                             const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;
    GEbounds     ebounds_src(const double& logEobs,
                             const double& theta = 0.0,
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0) const;

    // Methods
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFits& file);
    void                     write(GFitsBinTable& hdu) const;
    void                     save(const std::string& filename,
                                  const bool& clobber = false) const;
};


/***********************************************************************//**
 * @brief GCTAEdisp2D class extension
 ***************************************************************************/
%extend GCTAEdisp2D {
    GCTAEdisp2D copy() {
        return (*self);
    }
};
