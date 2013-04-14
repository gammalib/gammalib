/***************************************************************************
 *                   support module - Python bindings                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Jurgen Knodlseder                           *
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
 * swig -c++ -python -Wall support.i                                      *
 ***************************************************************************/
/**
 * @file support.i
 * @brief Support module
 * @author Juergen Knoedlseder
 */
%module support
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include "GException.hpp"
#include "GTools.hpp"
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
//%import(module="gammalib.base") "GContainer.i";

/* __ Support module _____________________________________________________ */
%include "GException.i"
%include "GTools.i"
%include "GNodeArray.i"
%include "GCsv.i"
%include "GRan.i"
%include "GUrl.i"
%include "GUrlFile.i"
%include "GUrlString.i"
