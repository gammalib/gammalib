/***************************************************************************
 *                        app.i - Application module                       *
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
 * swig -c++ -python -Wall app.i                                           *
 ***************************************************************************/
/**
 * @file app.i
 * @brief Application module
 * @author Juergen Knoedlseder
 */
%module app
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
#include "GException.hpp"
#include "GTools.hpp"
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i
%include std_string.i
%include std_vector.i
namespace std {
   %template(vectors) vector<std::string>;
};

/* __ Include GammaLib typemaps __________________________________________ */
%include typemap_GChatter.i
%include typemap_GFilename.i
%include typemap_slices.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.base") "GContainer.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Application module _________________________________________________ */
%include "GApplication.i"
%include "GLog.i"
%include "GApplicationPars.i"
%include "GApplicationPar.i"
