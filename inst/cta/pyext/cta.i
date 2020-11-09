/***************************************************************************
 *                cta.i - Cherenkov Telescope Array module                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2020 by Juergen Knoedlseder                         *
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
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall cta.i                                           *
 ***************************************************************************/
/**
 * @file cat.i
 * @brief Cherenkov Telescope Array module
 * @author Juergen Knoedlseder
 */
%module cta
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialComposite.hpp"
#include "GModelSpatialDiffuse.hpp"
#include "GModelSpatialDiffuseConst.hpp"
#include "GModelSpatialDiffuseCube.hpp"
#include "GModelSpatialDiffuseMap.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialEllipticalDisk.hpp"
#include "GModelSpatialEllipticalGauss.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialRadialDisk.hpp"
#include "GModelSpatialRadialGauss.hpp"
#include "GModelSpatialRadialRing.hpp"
#include "GModelSpatialRadialShell.hpp"
#include "GModelSpatialRadialProfile.hpp"
#include "GModelSpatialRadialProfileDMBurkert.hpp"
#include "GModelSpatialRadialProfileDMEinasto.hpp"
#include "GModelSpatialRadialProfileDMZhao.hpp"
#include "GModelSpatialRadialProfileGauss.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralBrokenPlaw.hpp"
#include "GModelSpectralComposite.hpp"
#include "GModelSpectralConst.hpp"
#include "GModelSpectralExpInvPlaw.hpp"
#include "GModelSpectralExponential.hpp"
#include "GModelSpectralExpPlaw.hpp"
#include "GModelSpectralFunc.hpp"
#include "GModelSpectralGauss.hpp"
#include "GModelSpectralLogParabola.hpp"
#include "GModelSpectralMultiplicative.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralPlawEnergyFlux.hpp"
#include "GModelSpectralPlawPhotonFlux.hpp"
#include "GModelSpectralSmoothBrokenPlaw.hpp"
#include "GModelSpectralSuperExpPlaw.hpp"
#include "GModelSpectralTable.hpp"
#include "GModelTemporal.hpp"
#include "GModelTemporalConst.hpp"
#include "GModelTemporalLightCurve.hpp"
#include "GModelTemporalPhaseCurve.hpp"
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Include GammaLib typemaps __________________________________________ */
%include typemap_GFilename.i
%include typemap_slices.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.base") "GContainer.i";
%import(module="gammalib.base") "GRegistry.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Inform about base classes __________________________________________ */
%import(module="gammalib.obs") "GObservation.i";
%import(module="gammalib.obs") "GEvent.i";
%import(module="gammalib.obs") "GEventBin.i";
%import(module="gammalib.obs") "GEventAtom.i";
%import(module="gammalib.obs") "GEvents.i";
%import(module="gammalib.obs") "GEventCube.i";
%import(module="gammalib.obs") "GEventList.i";
%import(module="gammalib.obs") "GResponse.i";
%import(module="gammalib.obs") "GInstDir.i";
%import(module="gammalib.obs") "GRoi.i";
%import(module="gammalib.model") "GModel.i";
%import(module="gammalib.model") "GModelData.i";
%import(module="gammalib.model") "GModelSpatial.i";
%import(module="gammalib.model") "GModelSpatialComposite.i";
%import(module="gammalib.model") "GModelSpatialDiffuse.i";
%import(module="gammalib.model") "GModelSpatialDiffuseConst.i";
%import(module="gammalib.model") "GModelSpatialDiffuseCube.i";
%import(module="gammalib.model") "GModelSpatialDiffuseMap.i";
%import(module="gammalib.model") "GModelSpatialElliptical.i";
%import(module="gammalib.model") "GModelSpatialEllipticalDisk.i";
%import(module="gammalib.model") "GModelSpatialEllipticalGauss.i";
%import(module="gammalib.model") "GModelSpatialPointSource.i";
%import(module="gammalib.model") "GModelSpatialRadial.i";
%import(module="gammalib.model") "GModelSpatialRadialDisk.i";
%import(module="gammalib.model") "GModelSpatialRadialGauss.i";
%import(module="gammalib.model") "GModelSpatialRadialRing.i";
%import(module="gammalib.model") "GModelSpatialRadialShell.i";
%import(module="gammalib.model") "GModelSpatialRadialProfile.i";
%import(module="gammalib.model") "GModelSpatialRadialProfileDMBurkert.i";
%import(module="gammalib.model") "GModelSpatialRadialProfileDMEinasto.i";
%import(module="gammalib.model") "GModelSpatialRadialProfileDMZhao.i";
%import(module="gammalib.model") "GModelSpatialRadialProfileGauss.i";
%import(module="gammalib.model") "GModelSpectral.i";
%import(module="gammalib.model") "GModelSpectralBrokenPlaw.i";
%import(module="gammalib.model") "GModelSpectralComposite.i";
%import(module="gammalib.model") "GModelSpectralConst.i";
%import(module="gammalib.model") "GModelSpectralExpInvPlaw.i";
%import(module="gammalib.model") "GModelSpectralExponential.i";
%import(module="gammalib.model") "GModelSpectralExpPlaw.i";
%import(module="gammalib.model") "GModelSpectralFunc.i";
%import(module="gammalib.model") "GModelSpectralGauss.i";
%import(module="gammalib.model") "GModelSpectralLogParabola.i";
%import(module="gammalib.model") "GModelSpectralMultiplicative.i";
%import(module="gammalib.model") "GModelSpectralNodes.i";
%import(module="gammalib.model") "GModelSpectralPlaw.i";
%import(module="gammalib.model") "GModelSpectralPlawEnergyFlux.i";
%import(module="gammalib.model") "GModelSpectralPlawPhotonFlux.i";
%import(module="gammalib.model") "GModelSpectralSmoothBrokenPlaw.i";
%import(module="gammalib.model") "GModelSpectralSuperExpPlaw.i";
%import(module="gammalib.model") "GModelSpectralTable.i";
%import(module="gammalib.model") "GModelTemporal.i";
%import(module="gammalib.model") "GModelTemporalConst.i";
%import(module="gammalib.model") "GModelTemporalLightCurve.i";
%import(module="gammalib.model") "GModelTemporalPhaseCurve.i";


/* __ Typemaps ___________________________________________________________ */
%typemap(out) GCTAAeff* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GCTAAeff;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GCTAPsf* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GCTAPsf;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GCTAEdisp* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GCTAEdisp;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GCTABackground* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GCTABackground;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GCTAResponse* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GCTAResponse;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GEvents* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GEvents;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GModelSpatial* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GModelSpatial;
    swig_cast_info *mycast = 0;
    mycast = myinfo->cast;
    while (mycast != 0) {
        if (strcmp(classname, mycast->type->name) == 0) {
            myinfo = mycast->type;
            break;
        }
        mycast = mycast->next;
    }
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GCTAModelSpatial* {
    if ($1 == NULL) {
        $result = Py_None;
        Py_INCREF(Py_None); // Py_None is a singleton so increment its reference
    }
    else {
        char classname[80];
        strcpy(classname, "_p_");
        strcat(classname, result->classname().c_str());
        swig_type_info *myinfo = SWIGTYPE_p_GCTAModelSpatial;
        swig_cast_info *mycast = 0;
        mycast = myinfo->cast;
        while (mycast != 0) {
            if (strcmp(classname, mycast->type->name) == 0) {
                myinfo = mycast->type;
                break;
            }
            mycast = mycast->next;
        }
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 | 0);
    }
}
%typemap(out) GCTAModelRadial* {
    if ($1 == NULL) {
        $result = Py_None;
        Py_INCREF(Py_None); // Py_None is a singleton so increment its reference
    }
    else {
        char classname[80];
        strcpy(classname, "_p_");
        strcat(classname, result->classname().c_str());
        swig_type_info *myinfo = SWIGTYPE_p_GCTAModelRadial;
        swig_cast_info *mycast = 0;
        mycast = myinfo->cast;
        while (mycast != 0) {
            if (strcmp(classname, mycast->type->name) == 0) {
                myinfo = mycast->type;
                break;
            }
            mycast = mycast->next;
        }
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 | 0);
    }
}
%typemap(out) GModelSpectral* {
    if ($1 == NULL) {
        $result = Py_None;
        Py_INCREF(Py_None); // Py_None is a singleton so increment its reference
    }
    else {
        char classname[80];
        strcpy(classname, "_p_");
        strcat(classname, result->classname().c_str());
        swig_type_info *myinfo = SWIGTYPE_p_GModelSpectral;
        swig_cast_info *mycast = 0;
        mycast = myinfo->cast;
        while (mycast != 0) {
            if (strcmp(classname, mycast->type->name) == 0) {
                myinfo = mycast->type;
                break;
            }
            mycast = mycast->next;
        }
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 | 0);
    }
}
%typemap(out) GModelTemporal* {
    if ($1 == NULL) {
        $result = Py_None;
        Py_INCREF(Py_None); // Py_None is a singleton so increment its reference
    }
    else {
        char classname[80];
        strcpy(classname, "_p_");
        strcat(classname, result->classname().c_str());
        swig_type_info *myinfo = SWIGTYPE_p_GModelTemporal;
        swig_cast_info *mycast = 0;
        mycast = myinfo->cast;
        while (mycast != 0) {
            if (strcmp(classname, mycast->type->name) == 0) {
                myinfo = mycast->type;
                break;
            }
            mycast = mycast->next;
        }
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 | 0);
    }
}

/* __ CTA ________________________________________________________________ */
%include "GCTAObservation.i"
%include "GCTAOnOffObservation.i"
%include "GCTAEventCube.i"
%include "GCTAEventList.i"
%include "GCTAEventBin.i"
%include "GCTAEventAtom.i"
%include "GCTAPointing.i"
%include "GCTAInstDir.i"
%include "GCTARoi.i"
%include "GCTAResponse.i"
%include "GCTAResponseIrf.i"
%include "GCTAResponseCube.i"
%include "GCTAResponseTable.i"
%include "GCTAAeff.i"
%include "GCTAAeffPerfTable.i"
%include "GCTAAeffArf.i"
%include "GCTAAeff2D.i"
%include "GCTAPsf.i"
%include "GCTAPsfPerfTable.i"
%include "GCTAPsfVector.i"
%include "GCTAPsf2D.i"
%include "GCTAPsfKing.i"
%include "GCTAPsfTable.i"
%include "GCTAEdisp.i"
%include "GCTAEdispPerfTable.i"
%include "GCTAEdispRmf.i"
%include "GCTAEdisp2D.i"
%include "GCTABackground.i"
%include "GCTABackgroundPerfTable.i"
%include "GCTABackground3D.i"
%include "GCTACubeExposure.i"
%include "GCTACubeBackground.i"
%include "GCTACubePsf.i"
%include "GCTACubeEdisp.i"
%include "GCTAModelBackground.i"
%include "GCTAModelCubeBackground.i"
%include "GCTAModelIrfBackground.i"
%include "GCTAModelAeffBackground.i"
%include "GCTAModelSpatial.i"
%include "GCTAModelSpatialLookup.i"
%include "GCTAModelSpatialGaussSpectrum.i"
%include "GCTAModelSpatialGradient.i"
%include "GCTAModelSpatialMultiplicative.i"
%include "GCTAModelSpatialRegistry.i"
%include "GCTAModelRadial.i"
%include "GCTAModelRadialRegistry.i"
%include "GCTAModelRadialGauss.i"
%include "GCTAModelRadialPolynom.i"
%include "GCTAModelRadialProfile.i"
%include "GCTAModelRadialAcceptance.i"
%include "GCTAModelSkyCube.i"
