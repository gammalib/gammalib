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

/* __ Typemaps ___________________________________________________________ */
%typemap(out) GModelSpatial* {
    if (dynamic_cast<GModelSpatialPointSource*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialPointSource, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpatialRadial*>($1) != NULL) {
        if (dynamic_cast<GModelSpatialRadialDisk*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialRadialDisk, 0 |  0 );
        }
        else if (dynamic_cast<GModelSpatialRadialGauss*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialRadialGauss, 0 |  0 );
        }
        else if (dynamic_cast<GModelSpatialRadialShell*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialRadialShell, 0 |  0 );
        }
        else {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialRadial, 0 |  0 );
        }
    }
    else if (dynamic_cast<GModelSpatialDiffuse*>($1) != NULL) {
        if (dynamic_cast<GModelSpatialDiffuseConst*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialDiffuseConst, 0 |  0 );
        }
        else if (dynamic_cast<GModelSpatialDiffuseCube*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialDiffuseCube, 0 |  0 );
        }
        else if (dynamic_cast<GModelSpatialDiffuseMap*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialDiffuseMap, 0 |  0 );
        }
        else {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialDiffuse, 0 |  0 );
        }
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatial, 0 |  0 );
    }
}
%typemap(out) GModelSpectral* {
    if (dynamic_cast<GModelSpectralConst*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralConst, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralFunc*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralFunc, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralNodes*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralNodes, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralPlaw*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralPlaw, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralPlaw2*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralPlaw2, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralLogParabola*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralLogParabola, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectral, 0 |  0 );
    }
}
%typemap(out) GModelTemporal* {
    if (dynamic_cast<GModelTemporalConst*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelTemporalConst, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelTemporal, 0 |  0 );
    }
}


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
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const;
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const;
    virtual double      npred(const GEnergy& obsEng,
                              const GTime& obsTime,
                              const GObservation& obs) const;
    virtual void        read(const GXmlElement& xml);
    virtual void        write(GXmlElement& xml) const;
    virtual std::string print(void) const;

    // Other methods
    GModelSpatial*      spatial(void) const;
    GModelSpectral*     spectral(void) const;
    GModelTemporal*     temporal(void) const;
    double              value(const GSkyDir& srcDir,
                              const GEnergy& srcEng,
                              const GTime& srcTime);
    GVector             gradients(const GSkyDir& srcDir,
                                  const GEnergy& srcEng,
                                  const GTime& srcTime);
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
