/***************************************************************************
 *                          model.i - Model module                         *
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
 * swig -c++ -python -Wall model.i                                         *
 ***************************************************************************/
/**
 * @file model.i
 * @brief Model module
 * @author Juergen Knoedlseder
 */
%module model
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
#include "GException.hpp"
#include "GTools.hpp"
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
%import(module="gammalib.opt") "GOptimizerPar.i";
%import(module="gammalib.opt") "GOptimizerPars.i";

/* __ Inform about other classes (needed for typecasting) ________________ */
%import(module="gammalib.sky") "GSkyRegion.i";

/* __ Typemaps ___________________________________________________________ */
%typemap(out) GModel* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GModel;
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
    if ($1 == NULL) {
        $result = Py_None;
        Py_INCREF(Py_None); // Py_None is a singleton so increment its reference
    }
    else {
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
%typemap(out) GSkyRegion* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GSkyRegion;
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

/* __ Model handling _____________________________________________________ */
%include "GModelPar.i"
%include "GModels.i"
%include "GModel.i"
%include "GModelRegistry.i"
%include "GModelSky.i"
%include "GModelData.i"
%include "GModelSpatial.i"
%include "GModelSpatialRegistry.i"
%include "GModelSpatialPointSource.i"
%include "GModelSpatialRadial.i"
%include "GModelSpatialRadialDisk.i"
%include "GModelSpatialRadialRing.i"
%include "GModelSpatialRadialGauss.i"
%include "GModelSpatialRadialShell.i"
%include "GModelSpatialRadialProfile.i"
%include "GModelSpatialRadialProfileGauss.i"
%include "GModelSpatialRadialProfileDMBurkert.i"
%include "GModelSpatialRadialProfileDMEinasto.i"
%include "GModelSpatialRadialProfileDMZhao.i"
%include "GModelSpatialElliptical.i"
%include "GModelSpatialEllipticalDisk.i"
%include "GModelSpatialEllipticalGauss.i"
%include "GModelSpatialDiffuse.i"
%include "GModelSpatialDiffuseConst.i"
%include "GModelSpatialDiffuseCube.i"
%include "GModelSpatialDiffuseMap.i"
%include "GModelSpatialComposite.i"
%include "GModelSpectral.i"
%include "GModelSpectralRegistry.i"
%include "GModelSpectralBrokenPlaw.i"
%include "GModelSpectralSmoothBrokenPlaw.i"
%include "GModelSpectralConst.i"
%include "GModelSpectralExpInvPlaw.i"
%include "GModelSpectralExponential.i"
%include "GModelSpectralExpPlaw.i"
%include "GModelSpectralSuperExpPlaw.i"
%include "GModelSpectralFunc.i"
%include "GModelSpectralGauss.i"
%include "GModelSpectralLogParabola.i"
%include "GModelSpectralNodes.i"
%include "GModelSpectralPlaw.i"
%include "GModelSpectralPlawPhotonFlux.i"
%include "GModelSpectralPlawEnergyFlux.i"
%include "GModelSpectralTable.i"
%include "GModelSpectralTablePar.i"
%include "GModelSpectralTablePars.i"
%include "GModelSpectralComposite.i"
%include "GModelSpectralMultiplicative.i"
%include "GModelTemporal.i"
%include "GModelTemporalRegistry.i"
%include "GModelTemporalConst.i"
%include "GModelTemporalLightCurve.i"
%include "GModelTemporalPhaseCurve.i"
