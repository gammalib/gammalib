/***************************************************************************
 *                       opt.i - Optimization module                       *
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
 * swig -c++ -python -Wall opt.i                                           *
 ***************************************************************************/
/**
 * @file opt.i
 * @brief Optimization module
 * @author Juergen Knoedlseder
 */
%module opt
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
%include typemap_GChatter.i
%include typemap_GFilename.i
%include typemap_slices.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.base") "GContainer.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Typemaps ___________________________________________________________ */
%typemap(out) GOptimizer* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GOptimizer;
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
%typemap(out) GOptimizerPar* {
    char classname[80];
    strcpy(classname, "_p_");
    strcat(classname, result->classname().c_str());
    swig_type_info *myinfo = SWIGTYPE_p_GOptimizerPar;
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

/* __ Optimizer module ___________________________________________________ */
%include "GOptimizer.i"
%include "GOptimizerLM.i"
%include "GOptimizerPar.i"
%include "GOptimizerPars.i"
%include "GOptimizerFunction.i"
%include "GPythonOptimizerFunction.i"
