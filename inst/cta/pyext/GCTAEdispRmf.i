/***************************************************************************
 *             GCTAEdispRmf.i - CTA RMF energy dispersion class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2018 by Christoph Deil & Ellis Owen                 *
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
 * @file GCTAEdispRmf.i
 * @brief CTA RMF energy dispersion class definition
 * @author Christoph Deil & Ellis Owen
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEdispRmf.hpp"
%}


/***********************************************************************//**
 * @class GCTAEdispRmf
 *
 * @brief CTA Redistribution Matrix File (RMF) energy dispersion class
 ***************************************************************************/
class GCTAEdispRmf : public GCTAEdisp {

public:
    // Constructors and destructors
    GCTAEdispRmf(void);
    explicit GCTAEdispRmf(const GFilename& filename);
    GCTAEdispRmf(const GCTAEdispRmf& psf);
    virtual ~GCTAEdispRmf(void);

    // Operators
    double operator()(const GEnergy& ereco,
                      const GEnergy& etrue,
                      const double&  theta = 0.0,
                      const double&  phi = 0.0,
                      const double&  zenith = 0.0,
                      const double&  azimuth = 0.0) const;

    // Implemented pure virtual methods
    void          clear(void);
    GCTAEdispRmf* clone(void) const;
    std::string   classname(void) const;
    void          load(const GFilename& filename);
    GFilename     filename(void) const;
    GEnergy       mc(GRan&          ran,
                     const GEnergy& etrue,
                     const double&  theta = 0.0,
                     const double&  phi = 0.0,
                     const double&  zenith = 0.0,
                     const double&  azimuth = 0.0) const;
    GEbounds      ebounds_obs(const double& logEsrc,
                              const double& theta = 0.0,
                              const double& phi = 0.0,
                              const double& zenith = 0.0,
                              const double& azimuth = 0.0) const;
    GEbounds      ebounds_src(const double& logEobs,
                              const double& theta = 0.0,
                              const double& phi = 0.0,
                              const double& zenith = 0.0,
                              const double& azimuth = 0.0) const;
    double        prob_erecobin(const GEnergy& ereco_min,
                                const GEnergy& ereco_max,
                                const GEnergy& etrue,
                                const double&  theta) const;

    // Other methods
    const GRmf& rmf(void) const;
};


/***********************************************************************//**
 * @brief Return filename
 ***************************************************************************/
%extend GCTAEdispRmf {
    GCTAEdispRmf copy() {
        return (*self);
    }
};
