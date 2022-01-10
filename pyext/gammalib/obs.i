/***************************************************************************
 *                       obs.i - Observation module                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2022 by Juergen Knoedlseder                         *
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
 * swig -c++ -python -Wall obs.i                                           *
 ***************************************************************************/
/**
 * @file obs.i
 * @brief Observation module
 * @author Juergen Knoedlseder
 */
%module obs
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
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Include GammaLib typemaps __________________________________________ */
%include typemap_GChatter.i
%include typemap_GFilename.i
%include typemap_slices.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.base") "GContainer.i";
%import(module="gammalib.base") "GRegistry.i";

/* __ Include optimizer class ____________________________________________ */
%import(module="gammalib.opt") "GOptimizerFunction.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Inform about base classes __________________________________________ */
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

/* __ Typemaps ___________________________________________________________ */
%typemap(out) GObservation* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GObservation;
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
%typemap(out) GEvent* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GEvent;
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
%typemap(out) GInstDir* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GInstDir;
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
%typemap(out) GRoi* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GRoi;
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
%typemap(out) GResponse* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GResponse;
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

/* __ Observation handling _______________________________________________ */
%include "GEnergy.i"
%include "GEnergies.i"
%include "GTime.i"
%include "GTimes.i"
%include "GTimeReference.i"
%include "GEbounds.i"
%include "GGti.i"
%include "GCaldb.i"
%include "GObservations.i"
%include "GObservation.i"
%include "GObservationRegistry.i"
%include "GEvents.i"
%include "GEventList.i"
%include "GEventCube.i"
%include "GEvent.i"
%include "GEventAtom.i"
%include "GEventBin.i"
%include "GInstDir.i"
%include "GRoi.i"
%include "GPhases.i"
%include "GResponse.i"
%include "GResponseCache.i"
%include "GResponseVectorCache.i"
%include "GPhoton.i"
%include "GPhotons.i"
%include "GSource.i"
%include "GPulsar.i"
%include "GPulsarEphemeris.i"
