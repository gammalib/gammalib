/***************************************************************************
 *                          model.i - Model module                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.base") "GContainer.i";
%import(module="gammalib.base") "GRegistry.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Inform about base classes __________________________________________ */
%import(module="gammalib.opt") "GOptimizerPar.i";
%import(module="gammalib.opt") "GOptimizerPars.i";

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
%typemap(out) GModelSpectral* {
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
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), myinfo, 0 |  0);
}
%typemap(out) GModelTemporal* {
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
%include "GModelSpatialRadialGauss.i"
%include "GModelSpatialRadialShell.i"
%include "GModelSpatialElliptical.i"
%include "GModelSpatialEllipticalDisk.i"
%include "GModelSpatialDiffuse.i"
%include "GModelSpatialDiffuseConst.i"
%include "GModelSpatialDiffuseCube.i"
%include "GModelSpatialDiffuseMap.i"
%include "GModelSpectral.i"
%include "GModelSpectralRegistry.i"
%include "GModelSpectralBrokenPlaw.i"
%include "GModelSpectralConst.i"
%include "GModelSpectralExpPlaw.i"
%include "GModelSpectralSuperExpPlaw.i"
%include "GModelSpectralFunc.i"
%include "GModelSpectralGauss.i"
%include "GModelSpectralLogParabola.i"
%include "GModelSpectralNodes.i"
%include "GModelSpectralPlaw.i"
%include "GModelSpectralPlaw2.i"
%include "GModelTemporal.i"
%include "GModelTemporalRegistry.i"
%include "GModelTemporalConst.i"
