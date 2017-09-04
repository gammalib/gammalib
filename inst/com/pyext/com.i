/***************************************************************************
 *                        com.i - COMPTEL module                           *
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
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall com.i                                           *
 ***************************************************************************/
/**
 * @file com.i
 * @brief COMPTEL module
 * @author Juergen Knoedlseder
 */
%module com
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

/* __ Typemaps ___________________________________________________________ */
%typemap(out) GCOMResponse* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GCOMResponse;
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

/* __ COMPTEL _____________________________________________________________ */
%include "GCOMObservation.i"
%include "GCOMResponse.i"
%include "GCOMEventBin.i"
%include "GCOMEventCube.i"
%include "GCOMEventList.i"
%include "GCOMEventAtom.i"
%include "GCOMRoi.i"
%include "GCOMInstDir.i"
%include "GCOMD1Response.i"
%include "GCOMD2Response.i"
%include "GCOMInstChar.i"
%include "GCOMIaq.i"
%include "GCOMModelDRBFitting.i"
%include "GCOMTim.i"
