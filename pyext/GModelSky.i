/***************************************************************************
 *                      GModelSky.i - Sky model class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSky.i
 * @brief Sky model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSky.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialRadialDisk.hpp"
#include "GModelSpatialRadialGauss.hpp"
#include "GModelSpatialRadialShell.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialDiffuse.hpp"
#include "GModelSpatialDiffuseConst.hpp"
#include "GModelSpatialDiffuseCube.hpp"
#include "GModelSpatialDiffuseMap.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralConst.hpp"
#include "GModelSpectralFunc.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralPlaw2.hpp"
#include "GModelSpectralLogParabola.hpp"
#include "GModelTemporal.hpp"
#include "GModelTemporalConst.hpp"
%}


/***********************************************************************//**
 * @class GModelSky
 *
 * @brief Sky model class
 ***************************************************************************/
class GModelSky : public GModel {
public:
    // Constructors and destructors
    GModelSky(void);
    explicit GModelSky(const std::string& type);
    explicit GModelSky(const GXmlElement& xml);
    explicit GModelSky(const GXmlElement& spatial,
                       const GXmlElement& spectral);
    explicit GModelSky(const GXmlElement& spatial,
                       const GXmlElement& spectral,
                       const GXmlElement& temporal);
    explicit GModelSky(const GModelSpatial& spatial,
                       const GModelSpectral& spectral);
    explicit GModelSky(const GModelSpatial& spatial,
                       const GModelSpectral& spectral,
                       const GModelTemporal& temporal);
    GModelSky(const GModelSky& model);
    virtual ~GModelSky(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GModelSky*  clone(void) const;
    virtual std::string type(void) const;
    virtual bool        is_constant(void) const;
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const;
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const;
    virtual double      npred(const GEnergy& obsEng,
                              const GTime& obsTime,
                              const GObservation& obs) const;
    virtual void        read(const GXmlElement& xml);
    virtual void        write(GXmlElement& xml) const;

    // Other methods
    GModelSpatial*      spatial(void) const;
    GModelSpectral*     spectral(void) const;
    GModelTemporal*     temporal(void) const;
    double              value(const GPhoton& photon);
    GVector             gradients(const GPhoton& photon);
    GPhotons            mc(const double& area,
                           const GSkyDir& dir, const double& radius,
                           const GEnergy& emin, const GEnergy& emax,
                           const GTime& tmin, const GTime& tmax,
                           GRan& ran) const;
};


/***********************************************************************//**
 * @brief GModelSky class extension
 ***************************************************************************/
%extend GModelSky {
};
