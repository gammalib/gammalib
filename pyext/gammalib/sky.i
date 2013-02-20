/***************************************************************************
 *                     sky module - Python bindings                        *
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
 * swig -c++ -python -Wall sky.i                                           *
 ***************************************************************************/
/**
 * @file sky.i
 * @brief Sky module
 * @author Juergen Knoedlseder
 */
%module sky
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
%import(module="gammalib.base") "GRegistry.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Skymap handling ____________________________________________________ */
%include "GSkyDir.i"
%include "GSkyPixel.i"
%include "GSkymap.i"
%include "GWcs.i"
%include "GWcsRegistry.i"
%include "GWcslib.i"
%include "GWcsAZP.i"
%include "GWcsCAR.i"
%include "GWcsTAN.i"
%include "GWcsSTG.i"
%include "GWcsHPX.i"
